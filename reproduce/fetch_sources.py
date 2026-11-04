"""Download the corpus's raw sources into data/, byte-identical to the ones the thesis used.

Fetches the inputs build_corpus.py reads.

    uv run reproduce/fetch_sources.py            # fetch missing data
    uv run reproduce/fetch_sources.py --check    # verify data

    data/natural_earth/ne_110m_admin_0_countries.{shp,shx,dbf,prj,cpg}   Natural Earth 1:110m Admin 0, v5.1.1
    data/fonts/*.ttf                                                     ten typefaces, Google Fonts

EVERY FILE IS PINNED BY SHA-256.  The corpus is defined by these exact bytes: a font
revision that moves one control point moves a glyph's vertices, and the stored exact
optima then describe a different polygon.  So a download whose hash does not match is
an error, never a warning, and nothing is written for it.  Upstream changing is exactly
the case this script exists to catch.

SOURCES.
  * Fonts: fonts.gstatic.com, at the VERSIONED URL each file was served from when the
    corpus was pinned (`/s/arbutus/v30/...`).  A new font release gets a new URL, so the
    pinned one keeps naming the old bytes for as long as Google keeps serving it, which
    it has historically done for a long time, but does not promise.  If a pinned URL fails
    or returns other bytes, the script falls back to the Google Fonts download list for the
    ten families, whose manifest names the CURRENT URL of every file.  That succeeds only 
    if the current file is still the pinned one; a new version there is reported as such.
    Each family's licence (OFL.txt, or LICENSE.txt for Fontdiner Swanky) is inlined in the
    same manifest and is written to data/fonts/licenses/<Family>/ beside the fonts.
  * Natural Earth: raw files from the nvkelso/natural-earth-vector repository at the tag
    v5.1.1.  The naciscdn.org zip serves the same bytes today, but it is "latest", not a
    version; the tag cannot move.
  * Zenodo, last for both: record 23938660 (DOI 10.5281/zenodo.23938660) holds every one
    of these files, byte-identical, with the font licences as font-licenses.zip.  It is
    the archival copy that keeps the corpus rebuildable if Google or GitHub stop serving
    the pinned bytes.  The record is pinned by version, not by its concept DOI
    (10.5281/zenodo.23938659), which would follow any later version.

Standard library only, so it runs before anything else is installed.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
import tempfile
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                           # reproduce -> repo
DATA = os.path.join(ROOT, "data")
FONTS_SUBDIR = "fonts"
NATURAL_EARTH_SUBDIR = "natural_earth"

USER_AGENT = "potato-peeling-fetch-sources/1.0"

GOOGLE_FONTS_LIST = (
    "https://fonts.google.com/download/list?family="
    "Arbutus,Audiowide,Fontdiner+Swanky,Matemasie,New+Rocker,Oi,Roboto,Story+Script,"
    "UnifrakturMaguntia,Vampiro+One"
)

GSTATIC = "https://fonts.gstatic.com/s/"

# file name in data/fonts/  ->  (versioned fonts.gstatic.com path, sha256 of the corpus's bytes)
FONTS = {
    "Arbutus-Regular.ttf": (
        GSTATIC + "arbutus/v30/NaPYcZ7dG_5J3poob9JtryO8fMU.ttf",
        "5ec84920a35b08f00d4a18c5d2a0e8a9bfe35041325c2fa9204513f89e50a17c",
    ),
    "Audiowide-Regular.ttf": (
        GSTATIC + "audiowide/v22/l7gdbjpo0cum0ckerWCtkQXPExpQBw.ttf",
        "434a720871336d359378beff5ebff3f9fd654d958693d272c7c6f2e271c7e41c",
    ),
    "FontdinerSwanky-Regular.ttf": (
        GSTATIC + "fontdinerswanky/v24/ijwOs4XgRNsiaI5-hcVb4hQgMvCD4uEfKiGvxts.ttf",
        "62863b9a99c674a8d5cc839f0b49bb30a2b0e8bc2bd86a812b4854c96b65f360",
    ),
    "Matemasie-Regular.ttf": (
        GSTATIC + "matemasie/v4/OD5BuMCN3ne3Gmr7dlL3rEq4DL6w2w.ttf",
        "7980c714c3dbf6f7fb54acb1319a09a16c2f1de1687915021c3eadb8064283e3",
    ),
    "NewRocker-Regular.ttf": (
        GSTATIC + "newrocker/v17/MwQzbhjp3-HImzcCU_cJkGMViblPtXs.ttf",
        "af01192aa3be57df4396f76600e8e38845896060550cee51525194ffcacf3884",
    ),
    "Oi-Regular.ttf": (
        GSTATIC + "oi/v21/w8gXH2EuRqtaut6yjBOG.ttf",
        "d666a21feb156b020209caf81c3e178ced1eda0de58941f3b82ac86e0eac411d",
    ),
    "Roboto-Bold.ttf": (
        GSTATIC + "roboto/v51/KFOMCnqEu92Fr1ME7kSn66aGLdTylUAMQXC89YmC2DPNWuYjammTggvWl0Qn.ttf",
        "4aaf8c5b661a386998c2e70cf2b87e2440f5404e0b8fd81164f0413fb3435ec6",
    ),
    "StoryScript-Regular.ttf": (
        GSTATIC + "storyscript/v3/mem5YaSw02SQ0OlzDuR8Isk-VeJoCqeDjg.ttf",
        "e1073bc62e39e61c1598ef4419e5cfe4d6d4d8d880ddb489be3791a1a797cf40",
    ),
    "UnifrakturMaguntia-Regular.ttf": (
        GSTATIC + "unifrakturmaguntia/v22/WWXPlieVYwiGNomYU-ciRLRvEmK7oaVun2xNNgNa1A.ttf",
        "84488aab45798072c293ceaaaa3e74cf59f109e343f4386df9a1dbd75294efea",
    ),
    "VampiroOne-Regular.ttf": (
        GSTATIC + "vampiroone/v19/gokqH6DoDl5yXvJytFsdLkqnsvhIor3K.ttf",
        "c8616348a896a2559b885576ab3cbaa635e5ea78d3f894e25ec7610c88e0c046",
    ),
}

ZENODO_RECORD = "23938660"                   # version record; concept DOI .23938659
ZENODO_FILES = f"https://zenodo.org/records/{ZENODO_RECORD}/files/"
ZENODO_LICENCES = "font-licenses.zip"


def zenodo_url(name: str) -> str:
    return f"{ZENODO_FILES}{name}?download=1"


NATURAL_EARTH_TAG = "v5.1.1"
NATURAL_EARTH_BASE = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    f"{NATURAL_EARTH_TAG}/110m_cultural/"
)
NATURAL_EARTH = {
    "ne_110m_admin_0_countries.shp": "08e341606e8391e458c3f08deb312de664b56bfae376064c5aa0aee6681a5f55",
    "ne_110m_admin_0_countries.shx": "8b0be2ad97dd484aee5c2ebc98697d5372e832b8ae58a35a661aeef6b985668d",
    "ne_110m_admin_0_countries.dbf": "1fee677cd4e03b367876e03861eb10197e4022a846bf92060e0313432863785b",
    "ne_110m_admin_0_countries.prj": "3259f0e55290a82b1350646f604e8a7ee1e2136c0320a40fad838ab40819fff8",
    "ne_110m_admin_0_countries.cpg": "3ad3031f5503a4404af825262ee8232cc04d4ea6683d42c5dd0a2f2a27ac9824",
}


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read()


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _matches(path: str, sha: str) -> bool:
    if not os.path.isfile(path):
        return False
    with open(path, "rb") as fh:
        return _sha256(fh.read()) == sha


def _write_atomic(path: str, data: bytes) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), prefix=".fetch-")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp, path)
    except BaseException:
        os.unlink(tmp)
        raise


def _attempt(url: str, sha: str) -> tuple[bytes | None, str]:
    """Download url; return (bytes, "") if they hash to sha, else (None, why not)."""
    try:
        data = _get(url)
    except OSError as exc:                      # URLError, HTTPError, timeouts
        return None, f"unreachable ({exc})"
    got = _sha256(data)
    if got != sha:
        return None, f"different bytes (sha256 {got})"
    return data, ""


def _store_first(dest: str, name: str, sha: str, sources) -> None:
    """Try each (label, url) in turn; write the first download whose hash matches.

    `sources` may be lazy: a url of None (e.g. a lookup that failed) is recorded and skipped.
    """
    tried = []
    for label, url in sources:
        if url is None:
            tried.append((label, "-", "not available"))
            continue
        data, why = _attempt(url, sha)
        if data is not None:
            _write_atomic(os.path.join(dest, name), data)
            print(f"  fetched   {name}" + ("" if not tried else f"  (from {label})"))
            return
        tried.append((label, url, why))
        print(f"  {label} for {name}: {why}")
    lines = "\n".join(f"  {lab:<10}{u}\n            {w}" for lab, u, w in tried)
    raise SystemExit(
        f"CANNOT FETCH {name}: no source serves the pinned bytes.\n{lines}\n  expected  {sha}\n"
        "Rebuilding the corpus from any other bytes would give a different corpus from the\n"
        "one the thesis measured."
    )


def font_manifest() -> dict:
    """The Google Fonts download list, with its `)]}'` anti-JSON-hijacking line removed."""
    text = _get(GOOGLE_FONTS_LIST).decode("utf-8")
    if text.startswith(")]}'"):
        text = text.split("\n", 1)[1]
    return json.loads(text)["manifest"]


def fetch_fonts(dest: str, force: bool) -> None:
    manifest = None                             # fetched only if something needs it

    def get_manifest() -> dict | None:
        nonlocal manifest
        if manifest is None:
            try:
                manifest = font_manifest()
            except (OSError, ValueError, KeyError) as exc:
                print(f"  Google Fonts list unavailable ({exc})")
                manifest = {}
        return manifest or None

    def current_url(name: str) -> str | None:
        m = get_manifest()
        for ref in (m or {}).get("fileRefs", []):
            if os.path.basename(ref["filename"]) == name:
                return ref["url"]
        return None

    def sources(name: str, pinned_url: str):
        yield "pinned", pinned_url
        url = current_url(name)
        yield "list", (url if url != pinned_url else None)
        yield "zenodo", zenodo_url(name)

    for name, (pinned_url, sha) in sorted(FONTS.items()):
        if not force and _matches(os.path.join(dest, name), sha):
            print(f"  ok        {name}")
            continue
        _store_first(dest, name, sha, sources(name, pinned_url))

    # Licences travel with the fonts (OFL 1.1 asks for it).  They are not corpus inputs,
    # so a failure to fetch them is a warning, and an existing set is not re-fetched.
    lic_dir = os.path.join(dest, "licenses")
    have = os.listdir(lic_dir) if os.path.isdir(lic_dir) else []
    if not force and len(have) >= len(FONTS):
        return
    m = get_manifest()
    if m:
        for f in m["files"]:
            family, base = f["filename"].split("/", 1)
            if base.upper() in ("OFL.TXT", "LICENSE.TXT"):
                _write_atomic(os.path.join(lic_dir, family, base), f["contents"].encode("utf-8"))
        return
    try:
        blob = _get(zenodo_url(ZENODO_LICENCES))
        with zipfile.ZipFile(io.BytesIO(blob)) as zf:
            for info in zf.infolist():
                parts = info.filename.split("/")
                if not info.is_dir() and ".." not in parts and not info.filename.startswith("/"):
                    _write_atomic(os.path.join(lic_dir, *parts), zf.read(info))
        print(f"  licences  from zenodo ({ZENODO_LICENCES})")
    except (OSError, zipfile.BadZipFile) as exc:
        print(f"  WARNING   could not fetch font licences ({exc})")


def fetch_natural_earth(dest: str, force: bool) -> None:
    for name, sha in NATURAL_EARTH.items():
        if not force and _matches(os.path.join(dest, name), sha):
            print(f"  ok        {name}")
            continue
        _store_first(dest, name, sha, [("github", NATURAL_EARTH_BASE + name), ("zenodo", zenodo_url(name))])


def check(dest: str) -> int:
    pinned = {os.path.join(FONTS_SUBDIR, n): s for n, (_, s) in FONTS.items()}
    pinned |= {os.path.join(NATURAL_EARTH_SUBDIR, n): s for n, s in NATURAL_EARTH.items()}
    bad = [n for n, s in pinned.items() if not _matches(os.path.join(dest, n), s)]
    for n in sorted(pinned):
        print(f"  {'BAD ' if n in bad else 'ok  '}      {n}")
    if bad:
        print(f"{len(bad)} file(s) missing or not the pinned bytes -- run without --check to fetch.")
        return 1
    print("All sources present and byte-identical to the pinned hashes.")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dest", default=DATA, help=f"directory to write into (default {DATA})")
    ap.add_argument("--check", action="store_true", help="verify hashes only; no network")
    ap.add_argument("--force", action="store_true", help="re-download even files that already match")
    args = ap.parse_args(argv)

    if args.check:
        return check(args.dest)
    ne_dir = os.path.join(args.dest, NATURAL_EARTH_SUBDIR)
    fonts_dir = os.path.join(args.dest, FONTS_SUBDIR)
    print(f"Natural Earth {NATURAL_EARTH_TAG} -> {ne_dir}")
    fetch_natural_earth(ne_dir, args.force)
    print(f"Google Fonts -> {fonts_dir}")
    fetch_fonts(fonts_dir, args.force)
    return check(args.dest)


if __name__ == "__main__":
    sys.exit(main())

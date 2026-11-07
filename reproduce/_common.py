"""Shared code for the reproduce/ scripts.

Each script produces one result with one command. Every run writes a new timestamped
database in results/ with an environment sidecar (machine, versions, git commit,
command line, start and end time), and prints a summary.
"""
from __future__ import annotations

import json
import os
import platform
import sqlite3
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                           # reproduce -> repo
DATA = os.path.join(ROOT, "data")
CORPUS = os.path.join(ROOT, "poly_db.sqlite")
RESULTS = os.path.join(HERE, "results")
PUBLISHED = os.path.join(HERE, "published")


def newest_matching(prefix: str) -> str | None:
    """Newest `<prefix>*.sqlite`, preferring a fresh run in results/ over published/."""
    import glob
    for base in (RESULTS, PUBLISHED):
        cands = sorted(glob.glob(os.path.join(base, prefix + "*.sqlite")))
        if cands:
            return cands[-1]
    return None


# --------------------------------------------------------------- capped solving
def _capped_target(fn, args, q):
    try:
        q.put(("ok", fn(*args)))
    except Exception as ex:                                    # noqa: BLE001
        q.put((f"error:{type(ex).__name__}", None))


def carry_over(stem: str, run, table: str, key: str = "pid") -> set:
    """For --resume: copy the newest earlier partial run of `stem` into this run and
    return the keys already done. Call right after entering Run()."""
    import glob
    priors = [p for p in sorted(glob.glob(os.path.join(RESULTS, stem + "_*.sqlite")))
              if os.path.abspath(p) != os.path.abspath(run.path)]
    if not priors:
        print(f"[resume] no prior {stem} partial found; starting fresh", flush=True)
        return set()
    prior = priors[-1]
    rows = sqlite3.connect(prior).execute(f"SELECT * FROM {table}").fetchall()
    for row in rows:
        run.db.execute(f"INSERT OR IGNORE INTO {table} VALUES ({','.join('?' * len(row))})",
                       row)
    run.db.commit()
    done = {r[0] for r in run.db.execute(f"SELECT DISTINCT {key} FROM {table}")}
    print(f"[resume] carried {len(rows)} rows from {os.path.basename(prior)}; "
          f"skipping {len(done)} done {key}s", flush=True)
    return done


def run_capped(fn, args, timeout):
    """Run `fn(*args)` in a forked subprocess with a wall-clock cap (a C++ call cannot
    be interrupted in-process). Returns (status, value, elapsed), status being 'ok',
    'timeout', 'error:<Type>' or 'nooutput'."""
    import multiprocessing as mp
    q = mp.Queue()
    p = mp.Process(target=_capped_target, args=(fn, args, q))
    t0 = time.time()
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate(); p.join()
        return ("timeout", None, float(timeout))
    el = time.time() - t0
    try:
        status, value = q.get_nowait()
    except Exception:                                         # noqa: BLE001
        return ("nooutput", None, el)
    return (status, value, el)


# ---------------------------------------------------------------- corpus access
def require_corpus() -> str:
    """The corpus path, exiting with instructions if it is missing (sqlite3.connect
    would silently create an empty file)."""
    if not os.path.exists(CORPUS):
        sys.exit(f"corpus not found at {CORPUS}\n"
                 f"Build it with reproduce/build_corpus.py (after reproduce/fetch_sources.py).")
    return CORPUS


def polygons(where: str = "", params: tuple = ()) -> list[tuple[int, str, list]]:
    """(pid, name, coords) for the corpus, largest ring only, no repeated last point."""
    from shapely import wkb as shapely_wkb
    conn = sqlite3.connect(require_corpus())
    sql = "SELECT id, name, geometry_wkb FROM polygons"
    if where:
        sql += " WHERE " + where
    sql += " ORDER BY id"
    out = []
    for pid, name, blob in conn.execute(sql, params):
        g = shapely_wkb.loads(blob)
        if g.geom_type == "MultiPolygon":
            g = max(g.geoms, key=lambda x: x.area)
        out.append((pid, name, list(g.exterior.coords)[:-1]))
    conn.close()
    return out


def chang_and_yap_optima(path: str | None = None):
    """(pid, name, nreflex, area_exact, P, Q) per optimum of a stored
    solve_chang_and_yap result (default: newest); P and Q as Fraction pairs."""
    from fractions import Fraction
    path = path or newest_matching("solve_chang_and_yap_")
    if not path:
        sys.exit("no solve_chang_and_yap result; run solve_chang_and_yap.py first")
    coords = {pid: co for pid, _, co in polygons()}
    conn = sqlite3.connect(path)
    rows = conn.execute(
        "SELECT pid, name, nreflex, area_exact, vertices_exact FROM optima "
        "WHERE status='ok' AND vertices_exact IS NOT NULL ORDER BY pid").fetchall()
    conn.close()
    out = []
    for pid, name, nreflex, area_exact, vx in rows:
        P = [(Fraction(x), Fraction(y)) for x, y in coords[pid]]
        Q = [(Fraction(x), Fraction(y)) for x, y in json.loads(vx)]
        out.append((pid, name, nreflex, area_exact, P, Q))
    return out


def polygon_stats(coords):
    """(n, nreflex, area, sha1) for an input polygon, exactly. The hash identifies the
    geometry independently of pid."""
    import hashlib
    from fractions import Fraction as Fr
    pts = [(Fr(str(x)), Fr(str(y))) for x, y in coords]
    n = len(pts)
    two_a = sum(pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1]
                for i in range(n))
    ccw = two_a > 0
    nref = 0
    for i in range(n):
        ax, ay = pts[(i - 1) % n]; bx, by = pts[i]; cx, cy = pts[(i + 1) % n]
        z = (bx - ax) * (cy - by) - (by - ay) * (cx - bx)
        if (z < 0) if ccw else (z > 0):
            nref += 1
    h = hashlib.sha1(repr([(str(x), str(y)) for x, y in pts]).encode()).hexdigest()
    return n, nref, abs(two_a) / 2, h


# ------------------------------------------------------------- result plumbing
def new_result_db(stem: str, schema: str) -> tuple[sqlite3.Connection, str]:
    """A fresh timestamped database.  Refuses to touch an existing one."""
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"{stem}_{datetime.now():%Y%m%d_%H%M%S}.sqlite")
    if os.path.exists(path):
        sys.exit(f"refusing to overwrite {path}")
    db = sqlite3.connect(path)
    db.executescript(schema)
    return db, path


def environment(argv: list[str] | None = None) -> dict:
    """Machine, versions, git commit and command line, recorded beside each result."""
    def run(cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=20).stdout.strip() or None
        except Exception:
            return None

    cpu = None
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip(); break
    except Exception:
        pass
    mem = None
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemTotal"):
                mem = line.split()[1] + " kB"; break
    except Exception:
        pass
    virt = None
    try:
        virt = "hypervisor" in open("/proc/cpuinfo").read()
    except Exception:
        pass

    return {
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cpu": cpu,
        "cpu_count": os.cpu_count(),
        "virtualised_guest": virt,
        "memory_total": mem,
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "compiler": (run(["g++", "--version"]) or "").split("\n")[0] or None,
        "cmake": (run(["cmake", "--version"]) or "").split("\n")[0] or None,
        "git_commit": run(["git", "-C", ROOT, "rev-parse", "HEAD"]),
        "git_dirty": bool(run(["git", "-C", ROOT, "status", "--porcelain"])),
        "argv": argv if argv is not None else sys.argv,
    }


def write_sidecar(result_path: str, env: dict, summary: dict | None = None) -> str:
    side = result_path.rsplit(".", 1)[0] + ".json"
    with open(side, "w") as fh:
        json.dump({"environment": env, "summary": summary or {}}, fh, indent=1)
    return side


def insert_row(db, table, values):
    """INSERT with one placeholder per value."""
    db.execute(f"INSERT INTO {table} VALUES ({','.join('?' * len(values))})", values)


class Run:
    """Context manager: opens the result database, records the environment, writes the
    sidecar on exit, prints where everything went."""

    def __init__(self, stem: str, schema: str):
        self.stem, self.schema = stem, schema
        self.summary: dict = {}

    def __enter__(self):
        self.db, self.path = new_result_db(self.stem, self.schema)
        self.env = environment()
        # Provenance in the database too, in case it is separated from its sidecar.
        self.db.execute("CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT)")
        corpus_sha = None
        try:
            import hashlib
            with open(CORPUS, "rb") as fh:
                corpus_sha = hashlib.sha1(fh.read()).hexdigest()[:16]
        except Exception:
            pass
        for k, v in (("script", self.stem), ("utc", self.env["utc"]),
                     ("git_commit", self.env["git_commit"]),
                     ("git_dirty", str(self.env["git_dirty"])),
                     ("corpus_path", os.path.relpath(CORPUS, ROOT)),
                     ("corpus_sha1_16", corpus_sha),
                     ("argv", " ".join(self.env["argv"])),
                     ("cpu", self.env["cpu"]), ("virtualised_guest", str(self.env["virtualised_guest"]))):
            self.db.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", (k, str(v)))
        self.db.commit()
        self.t0 = time.time()
        print(f"[{self.stem}] -> {os.path.relpath(self.path, ROOT)}", flush=True)
        return self

    def __exit__(self, *exc):
        self.db.commit()
        self.summary["elapsed_s"] = round(time.time() - self.t0, 1)
        side = write_sidecar(self.path, self.env, self.summary)
        self.db.close()
        print(f"\nwrote {os.path.relpath(self.path, ROOT)}"
              f"\n      {os.path.relpath(side, ROOT)}", flush=True)
        return False


def convex_hull_area(coords):
    """Exact area of the convex hull of the vertices (monotone chain)."""
    from fractions import Fraction as Fr
    pts = sorted(set((Fr(str(x)), Fr(str(y))) for x, y in coords))
    if len(pts) < 3:
        return Fr(0)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    hull = lower[:-1] + upper[:-1]
    n = len(hull)
    two_a = sum(hull[i][0] * hull[(i + 1) % n][1] - hull[(i + 1) % n][0] * hull[i][1]
                for i in range(n))
    return abs(two_a) / 2

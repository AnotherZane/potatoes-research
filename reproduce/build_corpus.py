"""Build the polygon corpus, poly_db.sqlite, from the downloaded sources.

    uv run reproduce/fetch_sources.py
    uv run reproduce/build_corpus.py

Reads data/natural_earth/ and data/fonts/, drops convex polygons (exact test), numbers
the rest in (source, source_id) order and writes them with their descriptors.  The
same sources always produce the same file.
"""
import argparse, json, os, shutil, sqlite3, sys, tempfile, time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import DATA, ROOT, convex_hull_area, environment, polygon_stats, write_sidecar

# Import parameters.  Roboto Bold was simplified at 6.0, every other font at 2.0.
FONT_PARAMS = {
    "Roboto-Bold.ttf": {"simplify_tol": 6.0},
}
FONT_DEFAULTS = {"em_size": 1000, "flatten_steps": 12, "simplify_tol": 2.0}

FONTS = os.path.join(DATA, "fonts")
SHAPEFILE = os.path.join(DATA, "natural_earth", "ne_110m_admin_0_countries.shp")
SCHEMA = os.path.join(ROOT, "database", "schema.sql")
CORPUS = os.path.join(ROOT, "poly_db.sqlite")

# What the pinned sources produce.  Any other count means the corpus was not reproduced.
EXPECTED_READ = 978
EXPECTED_KEPT = 824

sys.path.insert(0, os.path.join(ROOT, "importers"))
from font_glyphs import read_font            # noqa: E402
from natural_earth import read_countries     # noqa: E402


def is_convex_exact(coords) -> bool:
    """Exact: convex iff no corner turns the opposite way to the rest (collinear is fine)."""
    pts = [(Fr(str(x)), Fr(str(y))) for x, y in coords]
    n = len(pts)
    if n < 3:
        return True
    pos = neg = False
    for i in range(n):
        ax, ay = pts[(i - 1) % n]
        bx, by = pts[i]
        cx, cy = pts[(i + 1) % n]
        z = (bx - ax) * (cy - by) - (by - ay) * (cx - bx)
        if z > 0: pos = True
        elif z < 0: neg = True
        if pos and neg:
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=CORPUS,
                    help=f"corpus database to write (default {os.path.relpath(CORPUS, ROOT)})")
    a = ap.parse_args()

    if os.path.exists(a.out):
        sys.exit(f"refusing to overwrite {a.out}")
    for p in (SHAPEFILE, SCHEMA):
        if not os.path.exists(p):
            sys.exit(f"source missing: {p}\nRun reproduce/fetch_sources.py first.")
    fonts = sorted(f for f in os.listdir(FONTS) if f.lower().endswith((".ttf", ".otf")))
    if not fonts:
        sys.exit(f"no fonts in {FONTS}\nRun reproduce/fetch_sources.py first.")

    t0 = time.time()
    env = environment()
    tmpdir = tempfile.mkdtemp(prefix="corpus_")

    # Read every source polygon.
    imported = read_countries(SHAPEFILE)
    print(f"  read  {os.path.basename(SHAPEFILE)}: {len(imported)}", flush=True)
    for f in fonts:
        prm = dict(FONT_DEFAULTS); prm.update(FONT_PARAMS.get(f, {}))
        got = read_font(os.path.join(FONTS, f), **prm)
        print(f"  read  {f}: {len(got)}", flush=True)
        imported += got

    keys = [(r["source"], r["source_id"]) for r in imported]
    dupes = sorted({k for k in keys if keys.count(k) > 1})
    if dupes:
        sys.exit(f"duplicate (source, source_id): {dupes[:5]}.  Nothing written.")
    print(f"[build_corpus] read {len(imported)} polygons", flush=True)
    if len(imported) != EXPECTED_READ:
        sys.exit(f"read {len(imported)} polygons, expected {EXPECTED_READ}: the sources did not "
                 "reproduce the corpus.  Nothing written.")

    # Drop degenerate and convex polygons, then number by (source, source_id).
    from shapely import wkb as shapely_wkb
    rows = [(r["name"], r["source"], r["source_id"], r["geometry_wkb"], r["crs"],
             r["coordinate_space"], r["metadata"]) for r in imported]

    kept, dropped = [], []
    for row in rows:
        g = shapely_wkb.loads(row[3])
        gg = max(g.geoms, key=lambda x: x.area) if g.geom_type == "MultiPolygon" else g
        coords = list(gg.exterior.coords)[:-1]
        if len(coords) < 3:
            dropped.append((row[0], row[1], len(coords), "degenerate")); continue
        if is_convex_exact(coords):
            dropped.append((row[0], row[1], len(coords), "convex")); continue
        kept.append((row, coords))

    kept.sort(key=lambda kv: (kv[0][1], kv[0][2]))   # (source, source_id): name is NOT unique

    if len(kept) != EXPECTED_KEPT:
        sys.exit(f"kept {len(kept)} polygons, expected {EXPECTED_KEPT}.  Nothing written.")
    built = os.path.join(tmpdir, "corpus.sqlite")       # moved into place only once complete
    out = sqlite3.connect(built)
    out.executescript(open(SCHEMA).read())
    for new_id, (row, coords) in enumerate(kept, start=1):
        meta = dict(row[6])
        meta["vertices"] = len(coords)
        meta["convex_filtered"] = True
        # polygon_stats() is what the solve scripts record, so poly_sha matches result rows.
        n, nref, area, sha = polygon_stats(coords)
        hull = convex_hull_area(coords)
        xs = [x for x, _ in coords]; ys = [y for _, y in coords]
        q = [(Fr(str(x)), Fr(str(y))) for x, y in coords]
        two_a = sum(q[i][0] * q[(i + 1) % n][1] - q[(i + 1) % n][0] * q[i][1] for i in range(n))
        out.execute("INSERT INTO polygons (id,name,source,source_id,"
                    "geometry_wkb,crs,coordinate_space,metadata_json,"
                    "n_vertices,n_reflex,area_exact,area,hull_area_exact,hull_area,"
                    "orientation,min_x,min_y,max_x,max_y,poly_sha) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (new_id, row[0], row[1], row[2], row[3], row[4], row[5],
                     json.dumps(meta),
                     n, nref, str(area), float(area), str(hull), float(hull),
                     "ccw" if two_a > 0 else "cw", min(xs), min(ys), max(xs), max(ys), sha))
    out.commit()

    by_src = {}
    for row, _ in kept:
        by_src[row[1]] = by_src.get(row[1], 0) + 1
    drop_by = {}
    for _, s, _, why in dropped:
        drop_by[(s, why)] = drop_by.get((s, why), 0) + 1

    print(f"\n=== corpus built: {len(kept)} polygons ===")
    for s, c in sorted(by_src.items()):
        print(f"  kept     {s:16s} {c:5d}")
    for (s, why), c in sorted(drop_by.items()):
        print(f"  dropped  {s:16s} {c:5d}  ({why})")

    out.close()
    shutil.move(built, a.out)
    shutil.rmtree(tmpdir, ignore_errors=True)
    summary = {"kept": len(kept), "dropped": len(dropped),
               "kept_by_source": by_src,
               "dropped_by": {f"{s}:{w}": c for (s, w), c in drop_by.items()},
               "fonts": fonts, "elapsed_s": round(time.time() - t0, 1)}
    write_sidecar(a.out, env, summary)
    print(f"\nwrote {a.out}\n      {a.out.rsplit('.',1)[0]}.json")


if __name__ == "__main__":
    main()

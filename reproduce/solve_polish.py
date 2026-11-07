"""Run refinement (local ascent) from a stored solve result and store the results.

    uv run python reproduce/solve_polish.py --from cabello
    uv run python reproduce/solve_polish.py --from aronov

Refinement was called "polish" when these runs were made; the script, table and result
names keep that name.

The starting polygon is another algorithm's stored answer:

    --from cabello   solve_cabello_baseline_t1_*, its `vertices`
    --from aronov    solve_aronov_*, its `boundary`

The start is passed as rounded doubles, but start_area is the source's exact area and
the refined answer is computed and certified exactly, so only the starting hint is
rounded.

Four move sets, each run on every polygon so the comparison is paired:
    base     push + rotate (the defaults)
    adjoin   + adjoin_rounds=8
    split    + split_steps=3
    both     + both
"""
import argparse, json, os, sqlite3, sys, time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, polygons, polygon_stats, insert_row, newest_matching, carry_over

# --from -> (stored-result prefix, the rounded-vertex column to start from).
# start_area is read from the source's exact area_exact.
SOURCES = {
    "cabello": ("solve_cabello_baseline_t1_", "vertices"),
    "aronov":  ("solve_aronov_", "boundary"),
}
CONFIGS = {
    "base":   {},
    "adjoin": {"adjoin_rounds": 8},
    "split":  {"split_steps": 3},
    "both":   {"adjoin_rounds": 8, "split_steps": 3},
}
SCHEMA = """
CREATE TABLE polish (pid INT, name TEXT, n INT, nreflex INT, poly_sha TEXT, a_of_p TEXT,
                     source TEXT, config TEXT, status TEXT, elapsed REAL, error TEXT,
                     -- start_area: the source's exact area; area_exact: the refined one
                     start_area TEXT, area_exact TEXT,
                     start_vertices TEXT, vertices_exact TEXT,
                     PRIMARY KEY (pid, config));
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", choices=list(SOURCES), required=True)
    ap.add_argument("--configs", default=",".join(CONFIGS))
    ap.add_argument("--seed-result", default=None,
                    help="explicit stored seed result (default: newest matching --from)")
    ap.add_argument("--limit", type=int, help="first N polygons, smallest first (smoke test)")
    ap.add_argument("--resume", action="store_true",
                    help="carry the newest prior partial forward and skip done polygons")
    a = ap.parse_args()

    prefix, seed_col = SOURCES[a.src]
    seed_db = a.seed_result or newest_matching(prefix)
    if not seed_db or not os.path.exists(seed_db):
        sys.exit(f"no {a.src} solve result ({prefix}*.sqlite) in results/ or published/; "
                 f"run solve_{a.src}.py first")
    print(f"[solve_polish] seed: {os.path.basename(seed_db)} ({a.src})")

    from potato_peeling import refinement
    seeds = {}                    # pid -> (start_area_exact, rounded_seed_vertices_json)
    for pid, ae, vx in sqlite3.connect(seed_db).execute(
            f"SELECT pid, area_exact, {seed_col} FROM runs WHERE {seed_col} IS NOT NULL"):
        seeds[pid] = (ae, vx)

    rows = [(pid, nm, co) for pid, nm, co in polygons()]
    rows.sort(key=lambda t: (len(t[2]), t[0]))
    if a.limit:
        rows = rows[:a.limit]
    configs = [c for c in a.configs.split(",") if c in CONFIGS]

    with Run(f"solve_polish_from_{a.src}", SCHEMA) as r:
        print(f"[solve_polish] {len(rows)} polygons x {len(configs)} configs "
              f"({', '.join(configs)}), seeds available for {len(seeds)}", flush=True)
        done = carry_over(f"solve_polish_from_{a.src}", r, "polish") if a.resume else set()
        n_done = n_miss = 0
        gains = {c: [] for c in configs}
        for pid, name, P in rows:
            if pid in done:
                continue
            if pid not in seeds:
                n_miss += 1
                continue
            start_area_s, seed_vx = seeds[pid]
            start = [tuple(pt) for pt in json.loads(seed_vx)]   # rounded doubles
            n_v, nref, a_p, sha = polygon_stats(P)
            start_area = Fr(start_area_s) if start_area_s else None
            first = True
            for cfg in configs:
                t0 = time.time()
                try:
                    res = refinement.refine(P, start, **CONFIGS[cfg])
                    el = time.time() - t0
                    area = Fr(str(res["area_exact"]))
                    if start_area and start_area > 0:
                        gains[cfg].append(float(area / start_area))
                    insert_row(r.db, "polish",
                               (pid, name, len(P), nref, sha, str(a_p), a.src, cfg, "ok",
                                el, None, str(start_area) if start_area else None,
                                str(area),
                                json.dumps(start) if first else None,
                                json.dumps(res.get("vertices_exact"))))
                except Exception as ex:                     # noqa: BLE001
                    insert_row(r.db, "polish",
                               (pid, name, len(P), nref, sha, str(a_p), a.src, cfg,
                                "error", time.time() - t0, f"{type(ex).__name__}: {ex}",
                                str(start_area) if start_area else None, None,
                                None, None))
                    print(f"  ERROR [{pid}] {name} {cfg}: {type(ex).__name__}: {ex}",
                          flush=True)
                first = False
            n_done += 1
            r.db.commit()
            if n_done % 50 == 0:
                print(f"  ... {n_done}/{len(rows)}", flush=True)

        def stat(xs):
            xs = sorted(xs)
            return {} if not xs else {"n": len(xs), "min": xs[0], "median": xs[len(xs)//2],
                                      "mean": sum(xs)/len(xs), "max": xs[-1]}
        print(f"\n=== polish from {a.src}: {n_done} polygons, {n_miss} without a seed ===")
        print("GAIN (polished / seed area), all with a seed:")
        print(f"{'config':10s} {'n':>5s} {'min':>8s} {'median':>8s} {'mean':>8s} {'max':>8s}")
        for cfg in configs:
            g = stat(gains[cfg])
            if g:
                print(f"{cfg:10s} {g['n']:5d} {g['min']:8.4f} {g['median']:8.4f} "
                      f"{g['mean']:8.4f} {g['max']:8.4f}")
        r.summary = {"source": a.src, "seed_result": os.path.basename(seed_db),
                     "polygons": n_done, "without_seed": n_miss,
                     "gain": {c: stat(gains[c]) for c in configs}}


if __name__ == "__main__":
    sys.exit(main() or 0)

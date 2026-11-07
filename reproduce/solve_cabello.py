"""Solve Cabello et al. (2017) over the corpus and store the results.

    uv run python reproduce/solve_cabello.py --variant baseline --threads 1
    uv run python reproduce/solve_cabello.py --variant baseline --threads 8
    uv run python reproduce/solve_cabello.py --variant paper_literal --timeout 60

The answer depends only on the variant and the seed: sampling does not depend on the
thread count, so t1 and t8 runs store identical answers and different timings
(analyse_cabello_threads checks this). analyse_scaling uses the t1 timing, comparable
with single-threaded Chang & Yap. Do not run two solves at once.

Variants:
    baseline       max_edges=64 and reduced sample scales, the configuration that makes
                   the algorithm usable.
    paper_literal  the paper's parameters (max_edges=None, scales 1). Expected to be
                   infeasible, so run under a per-polygon --timeout; a timeout is the
                   result, not a failure.

Every completed answer must be convex and contained; a violation is a bug and the
script exits non-zero.
"""
import argparse, json, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, polygons, polygon_stats, insert_row, run_capped

# variant -> large_potato kwargs.  baseline is the shipped default (nothing overridden).
VARIANTS = {
    "baseline":      {},
    "paper_literal": {"max_edges": None, "sample_scale": 1.0, "s_scale": 1.0},
}

SCHEMA = """
CREATE TABLE runs (pid INT PRIMARY KEY, name TEXT, n INT, nreflex INT, poly_sha TEXT,
                   a_of_p TEXT,                       -- exact (num/den)
                   variant TEXT, threads INT, seed INT, status TEXT,
                   verts INT, area REAL, area_exact TEXT, convex INT, contained INT,
                   elapsed REAL, error TEXT,
                   -- the answer itself, rounded and exact (num/den)
                   vertices TEXT, vertices_exact TEXT,
                   -- the algorithm's counters; solve_seconds excludes binding overhead
                   r INT, r_ab INT, s_ab INT, iterations INT, phi_calls INT,
                   edges_searched INT, graph_edges INT, c2_effective REAL,
                   upper_bound REAL, upper_bound_certified INT, solve_seconds REAL,
                   edge_budget_exceeded INT);
"""


def _solve(coords, seed, threads, kwargs):
    """Top-level so run_capped can fork it for the paper-literal probe."""
    from potato_peeling.cabello import large_potato
    return large_potato(coords, seed=seed, threads=threads, **kwargs)


def _row(pid, name, coords, nref, sha, a_p, variant, threads, seed, status, el, res, error):
    """Build the full runs-table tuple; res is None for timeout/error rows."""
    if res is None:
        return (pid, name, len(coords), nref, sha, str(a_p), variant, threads, seed,
                status, None, None, None, None, None, el, error, None, None,
                *([None] * 12))
    verts = res.get("vertices") or []
    cv, ct = bool(res.get("convex")), bool(res.get("contained"))
    return (pid, name, len(coords), nref, sha, str(a_p), variant, threads, seed, status,
            len(verts), res.get("area"), str(res.get("area_exact")), int(cv), int(ct),
            el, error, json.dumps(verts), json.dumps(res.get("vertices_exact")),
            res.get("r"), res.get("r_ab"), res.get("s_ab"), res.get("iterations"),
            res.get("phi_calls"), res.get("edges_searched"), res.get("graph_edges"),
            res.get("c2_effective"), res.get("upper_bound"),
            int(bool(res.get("upper_bound_certified"))), res.get("seconds"),
            int(bool(res.get("edge_budget_exceeded"))))


def _check_certificate(res, pid, name):
    """A missing certificate is NOT a pass: defaulting an unreturned convex/contained
    flag to True would let a binding regression report a violation as success."""
    if "convex" not in res or "contained" not in res:
        raise KeyError(f"large_potato returned no convex/contained certificate "
                       f"(keys: {sorted(res)[:8]})")
    cv, ct = bool(res["convex"]), bool(res["contained"])
    if not (cv and ct):
        print(f"  NOT CONVEX/CONTAINED [{pid}] {name}: convex={cv} contained={ct}", flush=True)
    return cv and ct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=list(VARIANTS), default="baseline")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--seed", type=int, default=20260911)
    ap.add_argument("--limit", type=int, help="first N polygons, smallest first (smoke test)")
    ap.add_argument("--timeout", type=float, default=60.0,
                    help="per-polygon cap; used by the capped (paper_literal) path")
    ap.add_argument("--budget-hours", type=float, default=2.0,
                    help="total wall-clock for the capped (paper_literal) path")
    ap.add_argument("--nmax", type=int, default=0, help="skip n>nmax (0 = no limit)")
    a = ap.parse_args()

    kwargs = VARIANTS[a.variant]
    capped = (a.variant == "paper_literal")           # the infeasibility probe is capped
    stem = f"solve_cabello_{a.variant}_t{a.threads}"

    polys = polygons()
    polys.sort(key=lambda t: (len(t[2]), t[0]))        # ascending n (matters for the probe)
    if a.nmax:
        polys = [p for p in polys if len(p[2]) <= a.nmax]
    if a.limit:
        polys = polys[:a.limit]

    with Run(stem, SCHEMA) as r:
        print(f"[{stem}] {len(polys)} polygons, variant={a.variant}, threads={a.threads}, "
              f"seed={a.seed}" + (f", timeout={a.timeout}s, budget={a.budget_hours}h"
                                  if capped else ""), flush=True)
        n_ok = n_bad = n_err = n_to = 0
        slowest, total = [], 0.0
        t_start = time.time()
        for pid, name, coords in polys:
            if capped and time.time() - t_start > a.budget_hours * 3600:
                print(f"  budget reached after {n_ok+n_err+n_to} attempted", flush=True); break
            n_v, nref, a_p, sha = polygon_stats(coords)
            if capped:
                status, res, el = run_capped(_solve, (coords, a.seed, a.threads, kwargs),
                                             a.timeout)
                error = status if status.startswith("error") else None
            else:
                t0 = time.time()
                try:
                    res = _solve(coords, a.seed, a.threads, kwargs)
                    status, el, error = "ok", time.time() - t0, None
                except Exception as ex:                 # noqa: BLE001
                    status, res, el, error = "error", None, time.time() - t0, \
                        f"{type(ex).__name__}: {ex}"

            if status == "ok" and res is not None:
                total += el
                if _check_certificate(res, pid, name):
                    n_ok += 1
                else:
                    n_bad += 1
                slowest.append((el, pid, name, len(coords)))
            elif status == "timeout":
                n_to += 1
                print(f"  timeout [{pid}] {name} n={len(coords)} (>{a.timeout}s)", flush=True)
            else:
                n_err += 1
                print(f"  ERROR [{pid}] {name}: {error}", flush=True)
            insert_row(r.db, "runs",
                       _row(pid, name, coords, nref, sha, a_p, a.variant, a.threads,
                            a.seed, status, el, res, error))
            r.db.commit()

        slowest.sort(reverse=True)
        done = n_ok + n_bad
        print(f"\n=== {a.variant} @ {a.threads} thread{'s' if a.threads > 1 else ''}: "
              f"{done} completed, {n_bad} not convex/contained, {n_to} timed out, {n_err} errors ===")
        if done:
            print(f"  completed-solve total {total:8.1f}s   mean {total/done:6.3f}s")
        print("\nslowest completed:")
        for el, pid, name, n in slowest[:5]:
            print(f"  {el:8.2f}s  {pid:<5} {name:<28} n={n}")
        r.summary = {"variant": a.variant, "threads": a.threads, "seed": a.seed,
                     "polygons_attempted": len(polys), "completed": done,
                     "p1_violations": n_bad, "timed_out": n_to, "errors": n_err,
                     "total_solve_s": round(total, 1),
                     "mean_solve_s": round(total / done, 4) if done else None,
                     "slowest_s": round(slowest[0][0], 3) if slowest else None,
                     "capped": capped,
                     "timeout_s": a.timeout if capped else None}
    # Violations always fail; errors fail only when uncapped (paper_literal is capped).
    if n_bad:
        return 1
    if not capped and n_err:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

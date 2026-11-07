"""Compute a certified upper bound on A*(P) for every polygon and store it.

    uv run python reproduce/solve_upper_bound.py --budget-hours 4

U(P) is the minimum of the pivot bound and the width bound, each a valid upper bound on
A* (see upper_bound.py). It is an exact rational: only the search angles are floats,
and the slack is overstated so the bound stays sound. Runs in ascending n under a
per-polygon --timeout and a total --budget-hours; the largest polygons take minutes.

Each bound must satisfy U(P) <= area(P); a violation is a bug and the script exits
non-zero.
"""
import argparse, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, polygons, polygon_stats, insert_row, run_capped, carry_over
import upper_bound as ubmod          # upper_bound.py (validated)

SCHEMA = """
CREATE TABLE bounds (pid INT PRIMARY KEY, name TEXT, n INT, nreflex INT, poly_sha TEXT,
                     a_of_p TEXT, status TEXT, elapsed REAL, error TEXT,
                     -- all exact (num/den): the certified upper bound and its two parts
                     upper_exact TEXT, pivot_exact TEXT, width_exact TEXT);
"""


def _solve(coords, rel_tol, pairs, dirs):
    """Top-level so run_capped can fork it under a hard timeout."""
    from fractions import Fraction as Fr
    from potato_peeling.chang_and_yap import pre_checks
    clean = pre_checks(coords)
    pts = [(Fr(x), Fr(y)) for x, y in clean]
    area_p = abs(ubmod.shoelace(pts))
    bp = ubmod.upper_bound(clean, rel_tol=rel_tol, pairs=pairs)[0]
    bw = ubmod.width_bound(pts, n_uniform=dirs)
    if bw is None:
        bw = area_p
    return {"pivot": str(bp), "width": str(bw), "best": str(min(bp, bw)),
            "area_p": str(area_p)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget-hours", type=float, default=4.0)
    ap.add_argument("--timeout", type=float, default=1800.0, help="per-polygon, seconds")
    ap.add_argument("--rel-tol", type=float, default=1e-3)
    ap.add_argument("--pairs", type=int, default=6, help="two-pivot: top pivots to pair")
    ap.add_argument("--dirs", type=int, default=180, help="uniform slice directions")
    ap.add_argument("--limit", type=int, help="first N polygons, smallest first (smoke)")
    ap.add_argument("--resume", action="store_true",
                    help="carry the newest prior partial forward and skip done polygons")
    a = ap.parse_args()

    from fractions import Fraction as Fr
    polys = polygons()
    polys.sort(key=lambda t: (len(t[2]), t[0]))            # ascending n
    if a.limit:
        polys = polys[:a.limit]

    with Run("solve_upper_bound", SCHEMA) as r:
        print(f"[solve_upper_bound] {len(polys)} polygons, per-polygon timeout "
              f"{a.timeout}s, budget {a.budget_hours}h, pairs={a.pairs}", flush=True)
        done = carry_over("solve_upper_bound", r, "bounds") if a.resume else set()
        n_ok = n_to = n_err = n_bad = 0
        slowest = []
        t_start = time.time()
        for pid, name, coords in polys:
            if pid in done:
                continue
            if time.time() - t_start > a.budget_hours * 3600:
                print(f"  budget reached after {n_ok} bounded", flush=True); break
            n_v, nref, a_p, sha = polygon_stats(coords)
            status, res, el = run_capped(_solve, (coords, a.rel_tol, a.pairs, a.dirs),
                                         a.timeout)
            if status == "ok" and res is not None:
                upper = Fr(res["best"])
                if upper > Fr(res["area_p"]):
                    n_bad += 1
                    print(f"  UNSOUND [{pid}] {name}: U={float(upper):.4f} > "
                          f"area(P)={float(Fr(res['area_p'])):.4f}", flush=True)
                else:
                    n_ok += 1
                slowest.append((el, pid, name, len(coords)))
                insert_row(r.db, "bounds",
                           (pid, name, len(coords), nref, sha, str(a_p), "ok", el, None,
                            res["best"], res["pivot"], res["width"]))
            elif status == "timeout":
                n_to += 1
                insert_row(r.db, "bounds",
                           (pid, name, len(coords), nref, sha, str(a_p), "timeout", el,
                            None, None, None, None))
                print(f"  timeout [{pid}] {name} n={len(coords)} (>{a.timeout}s)", flush=True)
            else:
                n_err += 1
                insert_row(r.db, "bounds",
                           (pid, name, len(coords), nref, sha, str(a_p), "error", el,
                            status, None, None, None))
                print(f"  ERROR [{pid}] {name}: {status}", flush=True)
            r.db.commit()
            if (n_ok + n_to + n_err) % 50 == 0:
                print(f"  ... {n_ok} bounded, {n_to} timed out, n={len(coords)}", flush=True)

        slowest.sort(reverse=True)
        print(f"\n=== upper bound: {n_ok} bounded, {n_to} timed out, {n_err} errors, "
              f"{n_bad} UNSOUND ===")
        print("\nslowest bounded:")
        for el, pid, name, n in slowest[:5]:
            print(f"  {el:8.2f}s  {pid:<5} {name:<28} n={n}")
        r.summary = {"bounded": n_ok, "timed_out": n_to, "errors": n_err,
                     "unsound": n_bad, "pairs": a.pairs, "dirs": a.dirs,
                     "rel_tol": a.rel_tol, "per_polygon_timeout_s": a.timeout}
    return 1 if n_bad else 0            # timeouts expected; an unsound bound is a bug


if __name__ == "__main__":
    sys.exit(main())

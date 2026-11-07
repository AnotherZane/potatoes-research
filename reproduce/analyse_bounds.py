"""Compare the upper bound with A* and with the best lower bound, from stored results.

    uv run python reproduce/analyse_bounds.py

Where A* is known, how tight U is (U/A*, and the share of the area(P) - A* gap it
closes); elsewhere, the bracket [lower, U]. The lower bound is the largest certified
answer per polygon over Cabello, Aronov and both refinement runs. Exact (Fraction); a
violation of lower <= A* <= U is reported as a bug.
"""
import argparse, os, sqlite3, sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, newest_matching, insert_row

SCHEMA = """
CREATE TABLE bracket (pid INT PRIMARY KEY, n INT, a_of_p TEXT,
                      lower TEXT, a_star TEXT, upper TEXT, lower_source TEXT,
                      upper_over_astar REAL, gap_closed REAL, lower_over_upper REAL);
"""


def load(path, tbl, col="area_exact", where="status='ok'"):
    if not path:
        return {}
    return {r[0]: r[1] for r in sqlite3.connect(path).execute(
        f"SELECT pid,{col} FROM {tbl} WHERE {col} IS NOT NULL AND {where}")}


def load_polish_best(path):
    """pid -> largest refined area over the four configs."""
    if not path:
        return {}
    out = {}
    for pid, ae in sqlite3.connect(path).execute(
            "SELECT pid,area_exact FROM polish WHERE status='ok' AND area_exact IS NOT NULL"):
        f = Fr(ae)
        if pid not in out or f > out[pid]:
            out[pid] = f
    return out


def med(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else None


def main():
    argparse.ArgumentParser().parse_args()
    ub = newest_matching("solve_upper_bound_")
    if not ub:
        sys.exit("no solve_upper_bound result; run solve_upper_bound.py first")
    upper = {p: Fr(v) for p, v in load(ub, "bounds", "upper_exact").items()}
    area_p = {p: Fr(v) for p, v in load(ub, "bounds", "a_of_p", where="1=1").items()}
    n_of = {r[0]: r[1] for r in sqlite3.connect(ub).execute("SELECT pid,n FROM bounds")}
    print(f"[analyse_bounds] upper: {os.path.basename(ub)} ({len(upper)} bounded)")

    astar = {p: Fr(v) for p, v in
             load(newest_matching("solve_chang_and_yap_"), "optima").items()}

    # best lower bound per polygon, with its source
    lowers = {}
    for src, path, where in (
            ("cabello", newest_matching("solve_cabello_baseline_t1_"), "status='ok'"),
            ("aronov", newest_matching("solve_aronov_"), "1=1")):
        for p, v in load(path, "runs", where=where).items():
            f = Fr(v)
            if p not in lowers or f > lowers[p][0]:
                lowers[p] = (f, src)
    for src, path in (("polish_cabello", newest_matching("solve_polish_from_cabello_")),
                      ("polish_aronov", newest_matching("solve_polish_from_aronov_"))):
        for p, f in load_polish_best(path).items():
            if p not in lowers or f > lowers[p][0]:
                lowers[p] = (f, src)

    with Run("analyse_bounds", SCHEMA) as r:
        u_over_a, gap_closed, l_over_u = [], [], []
        n_order_viol = n_astar = n_bracket = 0
        for pid in sorted(upper):
            U = upper[pid]
            aP = area_p.get(pid)
            aS = astar.get(pid)
            lo, src = lowers.get(pid, (None, None))
            uoa = float(U / aS) if aS else None
            gc = (float((aP - U) / (aP - aS)) if (aS and aP and aP > aS) else None)
            lou = float(lo / U) if (lo and U) else None
            if aS:
                n_astar += 1
                u_over_a.append(U / aS)
                if gc is not None:
                    gap_closed.append(gc)
                if (lo is not None and lo > aS) or U < aS:      # exact ordering check
                    n_order_viol += 1
                    print(f"  ORDERING VIOLATION [{pid}]: lower={float(lo) if lo else None} "
                          f"A*={float(aS)} upper={float(U)}", flush=True)
            else:
                if lo is not None:
                    n_bracket += 1
                    l_over_u.append(lo / U)
            insert_row(r.db, "bracket",
                       (pid, n_of.get(pid), str(aP) if aP else None,
                        str(lo) if lo else None, str(aS) if aS else None, str(U), src,
                        uoa, gc, lou))
        r.db.commit()

        print(f"\n=== upper bound where A* is known ({n_astar}) ===")
        if u_over_a:
            print(f"  U/A*   : median {float(med(u_over_a)):.4f}  "
                  f"min {float(min(u_over_a)):.4f}  max {float(max(u_over_a)):.4f}")
        if gap_closed:
            print(f"  gap closed (of area(P)->A*): median {med(gap_closed):.3f}")
        print(f"  ordering violations (lower<=A*<=upper): {n_order_viol}")
        print(f"\n=== bracket where A* is UNKNOWN ({n_bracket}) -- the only statement there ===")
        if l_over_u:
            print(f"  lower/upper: median {float(med(l_over_u)):.4f}  "
                  f"min {float(min(l_over_u)):.4f}  (A* lies in [lower, upper])")
        r.summary = {"upper_result": os.path.basename(ub), "bounded": len(upper),
                     "with_astar": n_astar, "bracketed_without_astar": n_bracket,
                     "median_U_over_Astar": float(med(u_over_a)) if u_over_a else None,
                     "median_gap_closed": med(gap_closed) if gap_closed else None,
                     "median_lower_over_upper": float(med(l_over_u)) if l_over_u else None,
                     "ordering_violations": n_order_viol}
    return 1 if n_order_viol else 0


if __name__ == "__main__":
    sys.exit(main())

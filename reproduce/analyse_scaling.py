"""Compare the cost of Chang & Yap and Cabello, from stored solve results.

    uv run python reproduce/analyse_scaling.py

Fits elapsed time against n for the Chang & Yap solve and the 1-thread Cabello baseline
(comparable with the single-threaded exact solver), and reports the cost ratio where
both solved. Fits are reported with range, R^2 and n span: the exponent moves with the
fit window.
"""
import argparse, math, os, sqlite3, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, newest_matching, insert_row

SCHEMA = """
CREATE TABLE fits (algorithm TEXT, cut INT, exponent REAL, r2 REAL, points INT,
                   n_min INT, n_max INT);
CREATE TABLE both (pid INT PRIMARY KEY, n INT, cy_s REAL, cab_s REAL, ratio REAL);
"""


def fit_loglog(pairs):
    """Least squares on (log n, log t). Returns (exponent, R^2, count, n_min, n_max)."""
    pts = [(math.log(n), math.log(t)) for n, t in pairs if n > 0 and t and t > 0]
    k = len(pts)
    if k < 3:
        return None
    mx = sum(p[0] for p in pts) / k
    my = sum(p[1] for p in pts) / k
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    if sxx == 0:
        return None
    b = sxy / sxx
    aa = my - b * mx
    ss_res = sum((p[1] - (aa + b * p[0])) ** 2 for p in pts)
    ss_tot = sum((p[1] - my) ** 2 for p in pts)
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    ns = [n for n, t in pairs if n > 0 and t and t > 0]
    return b, r2, k, min(ns), max(ns)


def report_fits(db, label, algo, pairs, cuts):
    print(f"\n  {label}")
    print(f"    {'range':>10s} {'exponent':>10s} {'R^2':>7s} {'points':>7s} {'n-span':>12s}")
    out = []
    for lo in cuts:
        f = fit_loglog([(n, t) for n, t in pairs if n >= lo])
        if f:
            b, r2, k, n0, n1 = f
            print(f"    {'n >= '+str(lo):>10s} {b:10.2f} {r2:7.3f} {k:7d} {f'{n0}-{n1}':>12s}")
            insert_row(db, "fits", (algo, lo, round(b, 3), round(r2, 4), k, n0, n1))
            out.append({"cut": lo, "exponent": round(b, 3), "r2": round(r2, 4),
                        "points": k, "n_min": n0, "n_max": n1})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chang-and-yap", default=None, help="stored C&Y solve (default: newest)")
    ap.add_argument("--cabello", default=None,
                    help="stored Cabello BASELINE 1-thread solve (default: newest)")
    a = ap.parse_args()

    cy_path = a.chang_and_yap or newest_matching("solve_chang_and_yap_")
    cab_path = a.cabello or newest_matching("solve_cabello_baseline_t1_")
    if not cy_path:
        sys.exit("no solve_chang_and_yap result; run solve_chang_and_yap.py first")
    if not cab_path:
        sys.exit("no solve_cabello_baseline_t1 result; run solve_cabello.py --variant "
                 "baseline --threads 1 first")
    if "_t1_" not in os.path.basename(cab_path):
        sys.exit(f"refusing {os.path.basename(cab_path)}: the exact side is single-threaded, "
                 f"so the approximate side must be the 1-THREAD baseline")
    print(f"[analyse_scaling] exact : {os.path.basename(cy_path)}")
    print(f"[analyse_scaling] approx: {os.path.basename(cab_path)}")

    cy = {pid: (n, el) for pid, n, el in sqlite3.connect(cy_path).execute(
        "SELECT pid,n,elapsed FROM optima WHERE status='ok' AND elapsed IS NOT NULL")}
    cab = {pid: (n, el) for pid, n, el in sqlite3.connect(cab_path).execute(
        "SELECT pid,n,elapsed FROM runs WHERE status='ok' AND elapsed IS NOT NULL")}

    cy_pairs = [(n, el) for n, el in cy.values()]
    cab_pairs = [(n, el) for n, el in cab.values()]

    with Run("analyse_scaling", SCHEMA) as r:
        cuts = [3, 10, 20, 30]
        f_cy = report_fits(r.db, "Chang & Yap (exact, 1 thread)", "chang_and_yap",
                           cy_pairs, cuts)
        f_cab = report_fits(r.db, "Cabello (approximate, baseline, 1 thread)", "cabello",
                            cab_pairs, cuts)

        both = sorted(set(cy) & set(cab))
        for pid in both:
            n, cy_s = cy[pid]
            _, cab_s = cab[pid]
            insert_row(r.db, "both", (pid, n, cy_s, cab_s,
                                      (cy_s / cab_s) if cab_s else None))
        tot_cy = sum(cy[p][1] for p in both)
        tot_cab = sum(cab[p][1] for p in both)
        print(f"\n  on the {len(both)} polygons BOTH solved:")
        if both:
            print(f"    Chang & Yap total {tot_cy:9.1f}s   mean {tot_cy/len(both):7.3f}s")
            print(f"    Cabello     total {tot_cab:9.1f}s   mean {tot_cab/len(both):7.3f}s")
            if tot_cab:
                print(f"    ratio             {tot_cy/tot_cab:9.2f}x")
        r.summary = {"exact_result": os.path.basename(cy_path),
                     "approx_result": os.path.basename(cab_path),
                     "exact_solved": len(cy), "approx_solved": len(cab),
                     "both_solved": len(both),
                     "fits_chang_and_yap": f_cy, "fits_cabello": f_cab,
                     "total_cy_s": round(tot_cy, 1), "total_cab_s": round(tot_cab, 1),
                     "ratio_cy_over_cab": round(tot_cy / tot_cab, 2) if tot_cab else None}
    return 0


if __name__ == "__main__":
    sys.exit(main())

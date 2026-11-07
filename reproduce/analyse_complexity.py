"""Fit cost exponents and matrix densities, from stored solve results.

    uv run python reproduce/analyse_complexity.py

|Xi| against n, the A and M matrix densities (entries / |Xi|^2), the Chang & Yap time
exponent, and Cabello's first sample size r against n, each reported with its fit
range, R^2 and n span. |Xi| and the matrix counts come from solve_chang_and_yap, r from
solve_cabello. Memory is measured separately by solve_memory.py.
"""
import argparse, os, sqlite3, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, newest_matching, insert_row
from analyse_scaling import fit_loglog

SCHEMA = """
CREATE TABLE fits (metric TEXT, cut INT, exponent REAL, r2 REAL, points INT,
                   n_min INT, n_max INT);
CREATE TABLE density (pid INT PRIMARY KEY, n INT, n_xi INT, a_entries INT, m_entries INT,
                      a_density REAL, m_density REAL);
"""


def report(db, metric, pairs, cuts=(3, 10, 20, 30)):
    print(f"\n  {metric}")
    print(f"    {'range':>10s} {'exponent':>10s} {'R^2':>7s} {'points':>7s} {'n-span':>12s}")
    out = []
    for lo in cuts:
        f = fit_loglog([(n, v) for n, v in pairs if n >= lo])
        if f:
            b, r2, k, n0, n1 = f
            print(f"    {'n >= '+str(lo):>10s} {b:10.2f} {r2:7.3f} {k:7d} {f'{n0}-{n1}':>12s}")
            insert_row(db, "fits", (metric, lo, round(b, 3), round(r2, 4), k, n0, n1))
            out.append({"cut": lo, "exponent": round(b, 3), "r2": round(r2, 4),
                        "points": k, "n_min": n0, "n_max": n1})
    return out


def med(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else None


def main():
    argparse.ArgumentParser().parse_args()
    cy = newest_matching("solve_chang_and_yap_")
    cab = newest_matching("solve_cabello_baseline_t1_")
    if not cy:
        sys.exit("no solve_chang_and_yap result; run it first")
    print(f"[analyse_complexity] C&Y: {os.path.basename(cy)}")

    rows = sqlite3.connect(cy).execute(
        "SELECT pid,n,elapsed,n_xi,a_entries,m_entries FROM optima WHERE status='ok'").fetchall()
    time_pairs = [(n, el) for _, n, el, _, _, _ in rows if el]
    xi_pairs = [(n, x) for _, n, _, x, _, _ in rows if x]

    r_pairs = []
    if cab:
        for n, rr in sqlite3.connect(cab).execute(
                "SELECT n,r FROM runs WHERE status='ok' AND r IS NOT NULL"):
            r_pairs.append((n, rr))
        print(f"[analyse_complexity] Cabello: {os.path.basename(cab)}")

    with Run("analyse_complexity", SCHEMA) as r:
        f_time = report(r.db, "chang_and_yap_time", time_pairs)
        f_xi = report(r.db, "xi", xi_pairs)
        f_r = report(r.db, "cabello_r", r_pairs) if r_pairs else []

        adens, mdens = [], []
        for pid, n, el, xi, ae, me in rows:
            if not xi:
                continue
            ad = ae / xi ** 2 if ae is not None else None
            md = me / xi ** 2 if me is not None else None
            if ad is not None:
                adens.append(ad)
            if md is not None:
                mdens.append(md)
            insert_row(r.db, "density", (pid, n, xi, ae, me, ad, md))
        r.db.commit()

        print(f"\n  A/M matrix density (entries / |Xi|^2), over {len(adens)} solves:")
        if adens:
            print(f"    A: median {med(adens)*100:.1f}%  min {min(adens)*100:.1f}%  "
                  f"max {max(adens)*100:.1f}%")
        if mdens:
            print(f"    M: median {med(mdens)*100:.1f}%  min {min(mdens)*100:.1f}%  "
                  f"max {max(mdens)*100:.1f}%")
        r.summary = {"chang_and_yap_result": os.path.basename(cy),
                     "cabello_result": os.path.basename(cab) if cab else None,
                     "fits_time": f_time, "fits_xi": f_xi, "fits_cabello_r": f_r,
                     "a_density_median": med(adens), "m_density_median": med(mdens),
                     "solves": len(rows)}
    return 0


if __name__ == "__main__":
    sys.exit(main())

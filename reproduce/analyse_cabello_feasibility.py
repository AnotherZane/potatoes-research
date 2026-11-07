"""Show that the paper-literal Cabello parameters are infeasible, from stored results.

    uv run python reproduce/analyse_cabello_feasibility.py

Reads the paper_literal solve (run under a per-polygon timeout) and the baseline solve,
and reports how much of the corpus each completes, by n band, and the runtime on the
polygons both completed.

This measures time. Memory is the sharper limit: with max_edges=None and full sample
sizes, a single n=5 polygon was observed using over 36 GB, so no paper-literal run
completes to compare quality against. That observation is not re-run.
"""
import argparse, os, sqlite3, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, newest_matching, insert_row

SCHEMA = """
CREATE TABLE feasibility (pid INT PRIMARY KEY, n INT, pl_status TEXT, pl_elapsed REAL,
                          base_elapsed REAL, blowup REAL);
"""
BANDS = [(0, 10), (11, 20), (21, 30), (31, 50), (51, 999)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper-literal", default=None)
    ap.add_argument("--baseline", default=None,
                    help="baseline at the SAME thread count as the paper_literal run")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()

    pl = a.paper_literal or newest_matching(f"solve_cabello_paper_literal_t{a.threads}_")
    base = a.baseline or newest_matching(f"solve_cabello_baseline_t{a.threads}_")
    if not pl:
        sys.exit("no paper_literal solve; run solve_cabello.py --variant paper_literal "
                 f"--threads {a.threads} --timeout T first")
    if not base:
        sys.exit(f"no baseline t{a.threads} solve to compare against")
    print(f"[analyse_cabello_feasibility] paper-literal: {os.path.basename(pl)}")
    print(f"[analyse_cabello_feasibility] baseline     : {os.path.basename(base)}")

    plrows = {r[0]: (r[1], r[2], r[3]) for r in sqlite3.connect(pl).execute(
        "SELECT pid,n,status,elapsed FROM runs")}
    baseel = {r[0]: r[1] for r in sqlite3.connect(base).execute(
        "SELECT pid,elapsed FROM runs WHERE status='ok'")}

    with Run(f"analyse_cabello_feasibility_t{a.threads}", SCHEMA) as r:
        band_done = {b: 0 for b in BANDS}
        band_tot = {b: 0 for b in BANDS}
        blowups = []
        for pid, (n, status, el) in plrows.items():
            b = next(bd for bd in BANDS if bd[0] <= n <= bd[1])
            band_tot[b] += 1
            if status == "ok":
                band_done[b] += 1
            be = baseel.get(pid)
            bl = (el / be) if (status == "ok" and be and be > 0) else None
            if bl:
                blowups.append(bl)
            insert_row(r.db, "feasibility", (pid, n, status, el, be, bl))
        r.db.commit()

        attempted = len(plrows)
        completed = sum(1 for v in plrows.values() if v[1] == "ok")
        timed_out = sum(1 for v in plrows.values() if v[1] == "timeout")
        print(f"\n=== paper-literal Cabello: {completed}/{attempted} completed, "
              f"{timed_out} timed out (t{a.threads}) ===")
        print(f"{'n band':>10s} {'attempted':>10s} {'completed':>10s} {'rate':>7s}")
        for b in BANDS:
            if band_tot[b]:
                print(f"{f'{b[0]}-{b[1]}':>10s} {band_tot[b]:10d} {band_done[b]:10d} "
                      f"{band_done[b]/band_tot[b]:7.2f}")
        if blowups:
            blowups.sort()
            med = blowups[len(blowups)//2]
            print(f"\n  runtime blow-up (paper-literal / baseline) on {len(blowups)} both "
                  f"completed: median {med:.0f}x, max {blowups[-1]:.0f}x")
        r.summary = {"paper_literal_result": os.path.basename(pl),
                     "baseline_result": os.path.basename(base),
                     "threads": a.threads, "attempted": attempted,
                     "completed": completed, "timed_out": timed_out,
                     "completion_rate": round(completed / attempted, 3) if attempted else None,
                     "blowup_median": round(blowups[len(blowups)//2], 1) if blowups else None,
                     "blowup_max": round(max(blowups), 1) if blowups else None}
    return 0


if __name__ == "__main__":
    sys.exit(main())

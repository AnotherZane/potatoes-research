"""Look for suboptimal Chang & Yap answers with the T2 test, from stored results.

    uv run python reproduce/analyse_t2.py

T2 says the midpoint of every edge of an optimum lies in the convex hull of that edge's
contacts with dP. For each edge that fails it, t2_audit.audit_one builds the rotation T2
suggests and certifies, in exact rationals, that the result is convex, contained in P
and strictly larger. A certified body is a proof of suboptimality whatever T2's status;
T2 only says where to look.

P is the corpus polygon as the solver received it (its doubles as Fractions) and Q the
solver's vertices_exact.
"""
import argparse, json, os, sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, polygons, polygon_stats, newest_matching, insert_row
from t2_audit import audit_one          # the validated oracle

SCHEMA = """
CREATE TABLE t2 (pid INT PRIMARY KEY, name TEXT, n INT, nreflex INT, poly_sha TEXT,
                 status TEXT, edges INT, t2_fail INT, unconfirmed INT, beaten INT,
                 best_rel REAL, area_exact TEXT);
CREATE TABLE beaten (pid INT, name TEXT, edge INT, lo TEXT, hi TEXT,
                     old_area TEXT, new_area TEXT, rel_gain REAL, new_poly TEXT);
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--optima", default=None, help="stored C&Y solve (default: newest)")
    ap.add_argument("--limit", type=int)
    a = ap.parse_args()

    cy = a.optima or newest_matching("solve_chang_and_yap_")
    if not cy:
        sys.exit("no solve_chang_and_yap result; run it first (it supplies the optima)")
    print(f"[analyse_t2] optima: {os.path.basename(cy)}")

    coords = {pid: co for pid, nm, co in polygons()}     # P as C&Y received it
    import sqlite3
    rows = sqlite3.connect(cy).execute(
        "SELECT pid,name,area_exact,vertices_exact FROM optima "
        "WHERE status='ok' AND vertices_exact IS NOT NULL ORDER BY pid").fetchall()
    if a.limit:
        rows = rows[:a.limit]

    with Run("analyse_t2", SCHEMA) as r:
        print(f"[analyse_t2] {len(rows)} exact optima", flush=True)
        npoly = nt2 = nbeat = ntestfail = nunconf = 0
        worst = []
        for pid, name, ae, vx in rows:
            P = [(Fr(x), Fr(y)) for x, y in coords[pid]]     # dyadic: exactly C&Y's input
            Q = [(Fr(x), Fr(y)) for x, y in json.loads(vx)]
            res = audit_one(P, Q, expect_area=ae)
            if res["status"] == "DEGENERATE":
                continue
            npoly += 1
            n_v, nref, a_p, sha = polygon_stats(coords[pid])
            if res["status"] == "TESTFAIL":
                ntestfail += 1
                insert_row(r.db, "t2", (pid, name, n_v, nref, sha, "TESTFAIL", len(Q),
                                        0, 0, 0, 0.0, ae))
                print(f"  TESTFAIL [{pid}] {name} (skull did not re-certify)", flush=True)
                continue
            t2f, beat = res["t2_fail"], res["beaten"]
            nt2 += t2f
            nunconf += res.get("unconfirmed", 0)
            for (ei, lo, hi, g, delta, cand) in res["hits"]:
                a0 = res["area"]
                insert_row(r.db, "beaten",
                           (pid, name, ei, str(lo), str(hi), str(a0), str(a0 + g),
                            float(g / a0),
                            json.dumps([[str(x), str(y)] for x, y in cand])))
            if beat:
                nbeat += 1
                worst.append((res["best_rel"], pid, name, t2f, beat))
            insert_row(r.db, "t2", (pid, name, n_v, nref, sha, "OK", len(Q), t2f,
                                    res.get("unconfirmed", 0), beat, res["best_rel"], ae))
            r.db.commit()
            if npoly % 50 == 0:
                print(f"  ... {npoly}, {nt2} T2-failing edges, {nbeat} beaten", flush=True)

        print(f"\n=== T2 audit over {npoly} exact optima ===")
        print(f"  harness self-test failures (skipped)     {ntestfail}")
        print(f"  edges failing T2                         {nt2}")
        print(f"  of those, an improvement CERTIFIED       {nt2 - nunconf}")
        print(f"  OPTIMA BEATEN (convex+contained+larger)  {nbeat} / {npoly}")
        if worst:
            worst.sort(reverse=True)
            print("\n  largest relative gains (rel, pid, name, t2fail, beaten):")
            for w in worst[:12]:
                print(f"    {w[0]:.3e}  id {w[1]:<5} {w[2]:<22} t2fail {w[3]} beaten {w[4]}")
        r.summary = {"optima_result": os.path.basename(cy), "audited": npoly,
                     "testfail": ntestfail, "t2_failing_edges": nt2,
                     "certified_improvements": nt2 - nunconf, "optima_beaten": nbeat,
                     "beaten_ids": [w[1] for w in worst]}
    return 0


if __name__ == "__main__":
    sys.exit(main())

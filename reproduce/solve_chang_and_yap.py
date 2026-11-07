"""Solve Chang & Yap (1986), the exact algorithm, and store the results.

    uv run python reproduce/solve_chang_and_yap.py --budget-hours 6

The stored area_exact is A*(P), which the analyses divide by. The solver cannot reach
the whole corpus (it is memory-bound and its cost grows steeply with n), so it runs in
ascending n under a per-polygon --timeout and a total --budget-hours, and records
timeouts as results. Each completed answer is certified convex and contained in exact
arithmetic; a violation is a bug and the script exits non-zero.
"""
import argparse, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, polygons, polygon_stats, insert_row, run_capped, carry_over

SCHEMA = """
CREATE TABLE optima (pid INT PRIMARY KEY, name TEXT, n INT, nreflex INT, poly_sha TEXT,
                     a_of_p TEXT,                         -- exact (num/den)
                     status TEXT, elapsed REAL, error TEXT,
                     -- A*(P) and the optimum, exact (num/den)
                     area_exact TEXT, vertices_exact TEXT, convex INT, contained INT,
                     -- elapsed is the solve alone, child_wall the whole child process;
                     -- n_xi = |Xi|, a_entries/m_entries = entries in A and M
                     child_wall REAL, n_xi INT, a_entries INT, m_entries INT);
"""


def _solve(coords):
    """Top-level so run_capped can fork it under a hard timeout.

    convex_skull_exact is timed on its own, so `elapsed` excludes fork, import and
    telemetry overhead."""
    import json, time
    from potato_peeling.chang_and_yap import convex_skull_exact
    t0 = time.time()
    r = convex_skull_exact(coords)
    solve_s = time.time() - t0
    return {"area_exact": str(r.get("area_exact")),
            "vertices_exact": json.dumps(r.get("vertices_exact")),
            "convex": bool(r.get("convex")), "contained": bool(r.get("contained")),
            "solve_s": solve_s, "n_xi": r.get("n_xi"), "a_entries": r.get("a_entries"),
            "m_entries": r.get("m_entries")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget-hours", type=float, default=6.0)
    ap.add_argument("--timeout", type=float, default=600.0, help="per-polygon, seconds")
    ap.add_argument("--nmax", type=int, default=70, help="skip n>nmax (memory ceiling)")
    ap.add_argument("--limit", type=int, help="first N polygons, smallest first (smoke test)")
    ap.add_argument("--resume", action="store_true",
                    help="carry the newest prior partial forward and skip done polygons")
    a = ap.parse_args()

    polys = [(p, nm, c) for p, nm, c in polygons() if len(c) <= a.nmax]
    polys.sort(key=lambda t: (len(t[2]), t[0]))            # ascending n
    if a.limit:
        polys = polys[:a.limit]

    with Run("solve_chang_and_yap", SCHEMA) as r:
        print(f"[solve_chang_and_yap] {len(polys)} polygons with n<={a.nmax}, "
              f"per-polygon timeout {a.timeout}s, budget {a.budget_hours}h", flush=True)
        done = carry_over("solve_chang_and_yap", r, "optima") if a.resume else set()
        n_ok = n_to = n_err = n_bad = 0
        slowest = []
        t_start = time.time()
        for pid, name, coords in polys:
            if pid in done:
                continue
            if time.time() - t_start > a.budget_hours * 3600:
                print(f"  budget reached after {n_ok} solved", flush=True); break
            n_v, nref, a_p, sha = polygon_stats(coords)
            status, res, el = run_capped(_solve, (coords,), a.timeout)
            error = status if status.startswith("error") or status == "nooutput" else None
            cv = ct = None
            if status == "ok" and res is not None:
                cv, ct = int(bool(res["convex"])), int(bool(res["contained"]))
                if cv and ct:
                    n_ok += 1
                else:
                    n_bad += 1
                    print(f"  NOT CERTIFIED [{pid}] {name}: convex={cv} contained={ct}",
                          flush=True)
                slowest.append((res["solve_s"], pid, name, len(coords)))
                insert_row(r.db, "optima",
                           (pid, name, len(coords), nref, sha, str(a_p), "ok",
                            res["solve_s"], None, res["area_exact"], res["vertices_exact"],
                            cv, ct, el, res["n_xi"], res["a_entries"], res["m_entries"]))
            elif status == "timeout":
                n_to += 1
                insert_row(r.db, "optima",
                           (pid, name, len(coords), nref, sha, str(a_p), "timeout", None,
                            None, None, None, None, None, el, None, None, None))
            else:
                n_err += 1
                print(f"  ERROR [{pid}] {name}: {error}", flush=True)
                insert_row(r.db, "optima",
                           (pid, name, len(coords), nref, sha, str(a_p), "error", None,
                            error, None, None, None, None, el, None, None, None))
            r.db.commit()
            if (n_ok + n_to + n_err) % 50 == 0:
                print(f"  ... {n_ok} solved, {n_to} timed out, n={len(coords)}", flush=True)

        slowest.sort(reverse=True)
        ns = [n for _, _, _, n in slowest]
        print(f"\n=== exact optima: {n_ok} solved, {n_to} timed out, {n_err} errors, "
              f"{n_bad} NOT certified ===")
        if ns:
            print(f"  n range solved: {min(ns)}-{max(ns)}")
        print("\nslowest solved:")
        for el, pid, name, n in slowest[:5]:
            print(f"  {el:8.2f}s  {pid:<5} {name:<28} n={n}")
        tot = len(polygons())
        r.summary = {"corpus": tot, "solved": n_ok, "with_exact_optimum": n_ok,
                     "coverage_pct": round(100 * n_ok / tot, 1),
                     "timed_out": n_to, "errors": n_err, "not_certified": n_bad,
                     "n_solved": [min(ns), max(ns)] if ns else None,
                     "per_polygon_timeout_s": a.timeout, "nmax": a.nmax}
    return 1 if n_bad else 0            # timeouts/errors are expected; miscertification is a bug


if __name__ == "__main__":
    sys.exit(main())

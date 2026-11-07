"""Solve Aronov et al. (2011) over the corpus and store the results.

    uv run python reproduce/solve_aronov.py

Stores per polygon the exact answer, its exact area, the timing and the mesh statistics.
The solver checks each answer exactly before it reaches Python:

    O1  DP area == reconstructed area == boundary area
    O2  boundary simple, CCW and convex
    O3  the region is a non-empty set of mesh triangles

A failure is a bug and the script exits non-zero. Takes about a minute.
"""
import argparse, json, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, polygons, polygon_stats, insert_row

SCHEMA = """
CREATE TABLE runs (pid INT PRIMARY KEY, name TEXT, n INT, nreflex INT, poly_sha TEXT,
                   a_of_p TEXT, verts INT,           -- a_of_p exact (num/den)
                   area REAL, skull_valid INT, elapsed REAL, error TEXT,
                   -- the answer itself, rounded and exact
                   boundary TEXT, area_exact TEXT, boundary_exact TEXT,
                   -- mesh size, triangles used, and whether the answer is the whole mesh
                   faces INT, edges INT, used_faces INT, mesh_area REAL,
                   is_whole_mesh INT, mesh_valid INT);
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, help="first N polygons (smoke test)")
    a = ap.parse_args()

    from potato_peeling.aronov import skull_summary

    polys = polygons()
    polys.sort(key=lambda t: (len(t[2]), t[0]))     # ascending n, like the other solves
    if a.limit:
        polys = polys[:a.limit]

    with Run("solve_aronov", SCHEMA) as r:
        print(f"[solve_aronov] {len(polys)} polygons", flush=True)
        n_pass = n_fail = n_err = 0
        slowest = []
        for pid, name, coords in polys:
            n_v, nref, a_p, sha = polygon_stats(coords)
            t0 = time.time()
            try:
                s = skull_summary(coords)
                el = time.time() - t0
                ok = bool(s.get("skull_valid"))
                n_pass += ok
                n_fail += (not ok)
                insert_row(r.db, "runs",
                    (pid, name, len(coords), nref, sha, str(a_p),
                              s.get("vertices"), s.get("area"), int(ok), el, None,
                              json.dumps(s.get("boundary")),
                              str(s.get("area_exact")),
                              json.dumps(s.get("boundary_exact")),
                              s.get("faces"), s.get("edges"), s.get("used_faces"),
                              s.get("mesh_area"), int(bool(s.get("is_whole_mesh"))),
                              int(bool(s.get("mesh_valid")))))
                slowest.append((el, pid, name, len(coords)))
                if not ok:
                    print(f"  FAIL [{pid}] {name} n={len(coords)}", flush=True)
            except Exception as ex:
                n_err += 1
                insert_row(r.db, "runs",
                    (pid, name, len(coords), nref, sha, str(a_p),
                              None, None, 0, time.time() - t0,
                              f"{type(ex).__name__}: {ex}", None, None, None, *([None] * 6)))
                print(f"  ERROR [{pid}] {name}: {type(ex).__name__}: {ex}", flush=True)
            r.db.commit()

        slowest.sort(reverse=True)
        print(f"\n=== {n_pass}/{len(polys)} passed, {n_fail} failed, {n_err} errored ===")
        print("\nslowest:")
        for el, pid, name, n in slowest[:5]:
            print(f"  {el:8.2f}s  {pid:<5} {name:<28} n={n}")
        r.summary = {"polygons": len(polys), "passed": n_pass,
                     "failed": n_fail, "errored": n_err,
                     "slowest_s": round(slowest[0][0], 3) if slowest else None,
                     "slowest_n": slowest[0][3] if slowest else None}

    return 1 if (n_fail or n_err) else 0


if __name__ == "__main__":
    sys.exit(main())

"""Check that thesis_code/ reproduces the published results exactly.

    uv run python reproduce/verify_thesis_code.py aronov
    uv run python reproduce/verify_thesis_code.py cabello
    uv run python reproduce/verify_thesis_code.py chang_and_yap --budget 900
    uv run python reproduce/verify_thesis_code.py refinement --budget 900
    uv run python reproduce/verify_thesis_code.py all

Re-solves rows of the databases in reproduce/published/ and compares the exact answers:
area and vertices (Aronov: boundary), plus |Xi| and the A/M entry counts for Chang &
Yap. Cabello uses each row's published seed (baseline, 1 thread). Refinement starts
from the published Cabello and Aronov answers, as solve_polish.py does.

--budget caps the published solve time re-run per solver (per source and config for
refinement), cheapest rows first; without it everything is re-run (about 12 h for
refinement). Exits non-zero on any mismatch. Writes nothing.
"""
import argparse
import glob
import json
import os
import sqlite3
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import PUBLISHED, polygons  # noqa: E402


def published(prefix):
    found = sorted(glob.glob(os.path.join(PUBLISHED, prefix + "*.sqlite")))
    if not found:
        sys.exit(f"no published {prefix}*.sqlite")
    return found[-1]


def within(rows, budget, elapsed_col=1):
    """Rows (sorted cheapest first) whose cumulative published time fits the budget."""
    if budget is None:
        return rows
    out, spent = [], 0.0
    for r in rows:
        if spent + r[elapsed_col] > budget:
            break
        spent += r[elapsed_col]
        out.append(r)
    return out


def report(label, n, total, bad, t0):
    print(f"{label:28s} {n:4d}/{total} re-solved, {len(bad)} mismatches {bad[:8]}  "
          f"({time.time() - t0:.0f}s)", flush=True)
    return len(bad)


def check_aronov(coords, budget):
    from potato_peeling.aronov import skull_summary
    rows = sqlite3.connect(published("solve_aronov_")).execute(
        "SELECT pid, elapsed, area_exact, boundary_exact FROM runs "
        "WHERE area_exact IS NOT NULL ORDER BY elapsed").fetchall()
    run, t0, bad = within(rows, budget), time.time(), []
    for pid, _, ae, be in run:
        s = skull_summary(coords[pid])
        if (str(s["area_exact"]), json.dumps(s["boundary_exact"])) != (ae, be):
            bad.append(pid)
    return report("aronov", len(run), len(rows), bad, t0)


def check_cabello(coords, budget):
    from potato_peeling.cabello import large_potato
    rows = sqlite3.connect(published("solve_cabello_baseline_t1_")).execute(
        "SELECT pid, elapsed, area_exact, vertices_exact, seed FROM runs "
        "WHERE status='ok' ORDER BY elapsed").fetchall()
    run, t0, bad = within(rows, budget), time.time(), []
    for pid, _, ae, ve, seed in run:
        r = large_potato(coords[pid], seed=seed, threads=1)
        if (str(r["area_exact"]), json.dumps(r["vertices_exact"])) != (ae, ve):
            bad.append(pid)
    return report("cabello baseline", len(run), len(rows), bad, t0)


def check_chang_and_yap(coords, budget):
    from potato_peeling.chang_and_yap import convex_skull_exact
    rows = sqlite3.connect(published("solve_chang_and_yap_")).execute(
        "SELECT pid, elapsed, area_exact, vertices_exact, n_xi, a_entries, m_entries "
        "FROM optima WHERE status='ok' ORDER BY elapsed").fetchall()
    run, t0, bad = within(rows, budget), time.time(), []
    for pid, _, ae, ve, nxi, na, nm in run:
        r = convex_skull_exact(coords[pid])
        got = (str(r["area_exact"]), json.dumps(r["vertices_exact"]),
               r["n_xi"], r["a_entries"], r["m_entries"])
        if got != (ae, ve, nxi, na, nm):
            bad.append(pid)
    return report("chang_and_yap", len(run), len(rows), bad, t0)


def check_refinement(coords, budget):
    from potato_peeling import refinement
    from solve_polish import CONFIGS, SOURCES
    failures = 0
    for src, (seed_prefix, seed_col) in SOURCES.items():
        seeds = dict(sqlite3.connect(published(seed_prefix)).execute(
            f"SELECT pid, {seed_col} FROM runs WHERE {seed_col} IS NOT NULL"))
        pub = sqlite3.connect(published(f"solve_polish_from_{src}_"))
        starts = dict(pub.execute("SELECT pid, start_vertices FROM polish WHERE config='base'"))
        for cfg in CONFIGS:
            rows = pub.execute(
                "SELECT pid, elapsed, area_exact, vertices_exact FROM polish "
                "WHERE config=? AND status='ok' ORDER BY elapsed", (cfg,)).fetchall()
            run = within(rows, budget)
            t0, bad = time.time(), []
            for pid, _, ae, ve in run:
                start = [tuple(pt) for pt in json.loads(seeds[pid])]
                if json.dumps(start) != starts[pid]:
                    bad.append(pid)          # not the published starting polygon
                    continue
                res = refinement.refine(coords[pid], start, **CONFIGS[cfg])
                got = (str(Fr(str(res["area_exact"]))), json.dumps(res.get("vertices_exact")))
                if got != (ae, ve):
                    bad.append(pid)
            failures += report(f"refinement {src}/{cfg}", len(run), len(rows), bad, t0)
    return failures


CHECKS = {"aronov": check_aronov, "cabello": check_cabello,
          "chang_and_yap": check_chang_and_yap, "refinement": check_refinement}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("solver", choices=[*CHECKS, "all"])
    ap.add_argument("--budget", type=float, default=None,
                    help="seconds of published solve time to re-run (default: everything)")
    a = ap.parse_args()

    coords = {pid: co for pid, _, co in polygons()}
    names = list(CHECKS) if a.solver == "all" else [a.solver]
    failures = sum(CHECKS[name](coords, a.budget) for name in names)
    print("\nPASS: identical to the published results" if not failures
          else f"\nFAIL: {failures} mismatches")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

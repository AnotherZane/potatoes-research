"""Measure answer quality, from stored solve results.

    uv run python reproduce/analyse_quality.py

The ratio to A* of the Cabello and Aronov answers and of every refinement run, and
refinement's gain over its starting answer. A ratio needs A*, so it exists only where
Chang & Yap completed; a gain (refined / start) is reported for every polygon. Ratios
are computed in exact arithmetic (Fraction).
"""
import argparse, os, sqlite3, sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, polygons, newest_matching, insert_row

SCHEMA = """
CREATE TABLE quality (pid INT, n INT, source TEXT, config TEXT,
                      a_star TEXT, start_area TEXT, area_exact TEXT,
                      ratio REAL, start_ratio REAL, gain REAL,
                      PRIMARY KEY (pid, source, config));
"""


def stat(xs):
    xs = sorted(xs)
    if not xs:
        return {}
    n = len(xs)
    return {"n": n, "min": xs[0], "median": xs[n // 2], "mean": sum(xs) / n,
            "max": xs[-1], "below_0.95": sum(1 for x in xs if x < 0.95)}


def show(label, s):
    if s:
        print(f"{label:22s} {s['n']:5d} {s['min']:8.4f} {s['median']:8.4f} "
              f"{s['mean']:8.4f} {s['max']:8.4f} {s['below_0.95']:6d}")


def load_areas(path, tbl, col="area_exact", key="pid"):
    """pid -> exact area string, from a stored solve result."""
    if not path:
        return {}
    return {r[0]: r[1] for r in sqlite3.connect(path).execute(
        f"SELECT {key},{col} FROM {tbl} WHERE {col} IS NOT NULL")}


def main():
    argparse.ArgumentParser().parse_args()
    n_of = {pid: len(co) for pid, nm, co in polygons()}

    cy = newest_matching("solve_chang_and_yap_")
    if not cy:
        sys.exit("no solve_chang_and_yap result; run it first (A* is the denominator)")
    astar = {p: Fr(a) for p, a in load_areas(cy, "optima").items()}
    print(f"[analyse_quality] A*: {os.path.basename(cy)} ({len(astar)} optima)")

    cab = load_areas(newest_matching("solve_cabello_baseline_t1_"), "runs")
    aro = load_areas(newest_matching("solve_aronov_"), "runs")

    def load_polish(src):
        path = newest_matching(f"solve_polish_from_{src}_")
        if not path:
            return {}, None
        out = {}
        for pid, cfg, sa, ae in sqlite3.connect(path).execute(
                "SELECT pid,config,start_area,area_exact FROM polish "
                "WHERE status='ok' AND area_exact IS NOT NULL"):
            out.setdefault(pid, {})[cfg] = (sa, ae)
        return out, os.path.basename(path)

    pol_cab, pol_cab_f = load_polish("cabello")
    pol_aro, pol_aro_f = load_polish("aronov")

    with Run("analyse_quality", SCHEMA) as r:
        by = {}                       # (source,config) -> list of ratios
        gain_by = {}                  # (source,config) -> list of gains

        def record(pid, source, config, start_area, area):
            a_star = astar.get(pid)
            ratio = float(area / a_star) if a_star else None
            s_ratio = float(start_area / a_star) if (a_star and start_area) else None
            gain = float(area / start_area) if (start_area and start_area > 0) else None
            if ratio is not None:
                by.setdefault((source, config), []).append(ratio)
            if gain is not None:
                gain_by.setdefault((source, config), []).append(gain)
            insert_row(r.db, "quality",
                       (pid, n_of.get(pid), source, config,
                        str(a_star) if a_star else None,
                        str(start_area) if start_area is not None else None,
                        str(area), ratio, s_ratio, gain))

        # bodies with no separate "start": ratio to A* only
        for pid, ae in cab.items():
            record(pid, "cabello", "-", None, Fr(ae))
        for pid, ae in aro.items():
            record(pid, "aronov", "-", None, Fr(ae))
        # refinement runs: the start is the source's answer
        for src, pol in (("polish_from_cabello", pol_cab), ("polish_from_aronov", pol_aro)):
            for pid, cfgs in pol.items():
                for cfg, (sa, ae) in cfgs.items():
                    record(pid, src, cfg, Fr(sa) if sa else None, Fr(ae))
        r.db.commit()

        hdr = f"{'body / config':22s} {'n':>5s} {'min':>8s} {'median':>8s} {'mean':>8s} {'max':>8s} {'<0.95':>6s}"
        print(f"\n=== RATIO to A* (only where the exact optimum is known) ===\n{hdr}")
        show("cabello", stat(by.get(("cabello", "-"), [])))
        show("aronov", stat(by.get(("aronov", "-"), [])))
        for src in ("polish_from_cabello", "polish_from_aronov"):
            for cfg in ("base", "adjoin", "split", "both"):
                show(f"{src.split('_')[-1]}+{cfg}", stat(by.get((src, cfg), [])))

        print("\n=== GAIN (polished / seed area), corpus-wide ===")
        print(f"{'source / config':22s} {'n':>5s} {'min':>8s} {'median':>8s} {'mean':>8s} {'max':>8s}")
        for src in ("polish_from_cabello", "polish_from_aronov"):
            for cfg in ("base", "adjoin", "split", "both"):
                g = stat(gain_by.get((src, cfg), []))
                if g:
                    print(f"{src.split('_')[-1]+'/'+cfg:22s} {g['n']:5d} {g['min']:8.4f} "
                          f"{g['median']:8.4f} {g['mean']:8.4f} {g['max']:8.4f}")

        r.summary = {
            "a_star_result": os.path.basename(cy), "a_star_count": len(astar),
            "polish_cabello_result": pol_cab_f, "polish_aronov_result": pol_aro_f,
            "ratio": {f"{s}|{c}": stat(v) for (s, c), v in by.items()},
            "gain": {f"{s}|{c}": stat(v) for (s, c), v in gain_by.items()},
        }
    return 0


if __name__ == "__main__":
    sys.exit(main())

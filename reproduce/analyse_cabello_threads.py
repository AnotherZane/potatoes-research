"""Check that 1-thread and 8-thread Cabello runs give identical answers.

    uv run python reproduce/analyse_cabello_threads.py

Reads the baseline t1 and t8 solves. Sampling does not depend on the thread count, so
every area must be the same rational; a difference is a bug and the script exits
non-zero.
"""
import argparse, os, sqlite3, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import newest_matching


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t1", default=None)
    ap.add_argument("--t8", default=None)
    a = ap.parse_args()
    t1 = a.t1 or newest_matching("solve_cabello_baseline_t1_")
    t8 = a.t8 or newest_matching("solve_cabello_baseline_t8_")
    if not t1 or not t8:
        sys.exit("need both baseline t1 and t8 solve results "
                 "(run solve_cabello.py --variant baseline at both thread counts)")
    print(f"[analyse_cabello_threads] t1: {os.path.basename(t1)}")
    print(f"[analyse_cabello_threads] t8: {os.path.basename(t8)}")

    def areas(p):
        return {r[0]: r[1] for r in sqlite3.connect(p).execute(
            "SELECT pid,area_exact FROM runs WHERE status='ok' AND area_exact IS NOT NULL")}
    A, B = areas(t1), areas(t8)
    common = sorted(set(A) & set(B))
    diff = [p for p in common if A[p] != B[p]]
    print(f"\n{len(common)} polygons completed in both")
    print(f"  areas identical : {len(common) - len(diff)}")
    print(f"  areas DIFFERING : {len(diff)}   <-- must be 0: rng.hpp makes sampling "
          f"independent of thread count")
    for p in diff[:10]:
        print(f"    [{p}]: {str(A[p])[:28]} vs {str(B[p])[:28]}")
    return 1 if diff else 0


if __name__ == "__main__":
    sys.exit(main())

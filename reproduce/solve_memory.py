"""Measure the memory of Chang & Yap solves and store it.

    uv run python reproduce/solve_memory.py --budget-hours 2

Memory, not time, limits Chang & Yap in practice. For each polygon this records peak
RSS (VmHWM, what the OOM killer acts on) and, per phase, the stored-entry and bucket
counts of the A and M matrices, separating the per-entry lazy-exact DAGs from the
reserved bucket arrays.

Each solve runs in a fresh subprocess: peak RSS never falls, so two solves in one
process would both report the larger, and a fork would inherit the parent's RSS. An OOM
kill is recorded as a row with no peak. By default it walks a ladder of one polygon per
n-band up to --max-n; --all does every polygon.
"""
import argparse, json, os, subprocess, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import Run, polygons, polygon_stats, insert_row, carry_over

GB = 1024 ** 3
SCHEMA = """
CREATE TABLE mem (pid INT PRIMARY KEY, name TEXT, n INT, nreflex INT, poly_sha TEXT,
                  a_of_p TEXT, status TEXT, peak_rss INT, baseline_rss INT, xi INT,
                  a_stored INT, a_buckets INT, m_stored INT, m_buckets INT,
                  seconds REAL, area_exact TEXT, error TEXT);
"""
SENTINEL = "@@MEM@@"


def run_child(pid):
    """One clean solve with the memory probe on; print the report as a sentinel line."""
    from potato_peeling import chang_and_yap as cy
    co = polygons(where="id=?", params=(pid,))[0][2]
    cy.set_mem_probe(True)
    t0 = time.time()
    res = cy.convex_skull_exact(co)
    rep = cy.mem_report()
    rep["seconds"] = time.time() - t0
    rep["area_exact"] = str(res.get("area_exact"))
    print(SENTINEL + json.dumps(rep))


def ladder(polys, max_n):
    """One polygon per n-decade, ascending -- the space trend without the whole corpus."""
    picked, seen = [], set()
    for pid, name, co in polys:
        n = len(co)
        if n > max_n:
            break
        b = n // 10
        if b in seen:
            continue
        seen.add(b)
        picked.append((pid, name, co))
    return picked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--child", type=int, default=None, help=argparse.SUPPRESS)
    ap.add_argument("--all", action="store_true", help="every polygon, not just the ladder")
    ap.add_argument("--max-n", type=int, default=140)
    ap.add_argument("--timeout", type=float, default=900.0, help="per-polygon, seconds")
    ap.add_argument("--budget-hours", type=float, default=2.0)
    ap.add_argument("--limit", type=int, help="first N (smoke test)")
    ap.add_argument("--resume", action="store_true",
                    help="carry the newest prior partial forward and skip done polygons")
    a = ap.parse_args()

    if a.child is not None:
        run_child(a.child)
        return 0

    polys = polygons()
    polys.sort(key=lambda t: (len(t[2]), t[0]))
    polys = [p for p in polys if len(p[2]) <= a.max_n]
    if not a.all:
        polys = ladder(polys, a.max_n)
    if a.limit:
        polys = polys[:a.limit]

    with Run("solve_memory", SCHEMA) as r:
        print(f"[solve_memory] {len(polys)} polygons ({'all' if a.all else 'ladder'}, "
              f"n<={a.max_n}), per-polygon timeout {a.timeout}s, budget {a.budget_hours}h",
              flush=True)
        done = carry_over("solve_memory", r, "mem") if a.resume else set()
        n_ok = n_oom = n_to = n_err = 0
        t_start = time.time()
        for pid, name, co in polys:
            if pid in done:
                continue
            if time.time() - t_start > a.budget_hours * 3600:
                print(f"  budget reached after {n_ok} measured", flush=True); break
            n_v, nref, a_p, sha = polygon_stats(co)
            try:
                proc = subprocess.run(
                    [sys.executable, os.path.abspath(__file__), "--child", str(pid)],
                    capture_output=True, text=True, timeout=a.timeout)
                line = next((l for l in proc.stdout.splitlines()
                             if l.startswith(SENTINEL)), None)
            except subprocess.TimeoutExpired:
                proc, line = None, None
            if line:
                rep = json.loads(line[len(SENTINEL):])
                ph = rep.get("phases") or [{}]
                last = ph[-1] if ph else {}
                peak = rep.get("peak_rss")
                insert_row(r.db, "mem",
                           (pid, name, len(co), nref, sha, str(a_p), "ok", peak,
                            rep.get("baseline_rss"), last.get("xi"), last.get("a_stored"),
                            last.get("a_buckets"), last.get("m_stored"),
                            last.get("m_buckets"), rep.get("seconds"),
                            rep.get("area_exact"), None))
                n_ok += 1
                print(f"  [{pid}] {name} n={len(co)}  peak {peak/GB:.3f} GiB  "
                      f"{rep.get('seconds', 0):.1f}s", flush=True)
            else:
                # No report: a timeout, or an OOM kill (rc=-9).  The ceiling IS the datum.
                if proc is None:
                    status = "timeout"; n_to += 1
                elif proc.returncode == -9:
                    status = "oom"; n_oom += 1
                else:
                    status = "error"; n_err += 1
                err = None if proc is None else (proc.stderr or "")[-160:]
                insert_row(r.db, "mem",
                           (pid, name, len(co), nref, sha, str(a_p), status, None, None,
                            None, None, None, None, None, None, None, err))
                print(f"  [{pid}] {name} n={len(co)}  {status.upper()}"
                      + (f" rc={proc.returncode}" if proc else ""), flush=True)
            r.db.commit()

        print(f"\n=== memory ladder: {n_ok} measured, {n_oom} OOM, {n_to} timed out, "
              f"{n_err} errors ===")
        rows = list(r.db.execute("SELECT n,peak_rss FROM mem WHERE status='ok' "
                                 "AND peak_rss IS NOT NULL ORDER BY n"))
        if rows:
            print(f"  peak RSS: n={rows[0][0]} -> {rows[0][1]/GB:.3f} GiB,  "
                  f"n={rows[-1][0]} -> {rows[-1][1]/GB:.3f} GiB")
        r.summary = {"measured": n_ok, "oom": n_oom, "timed_out": n_to, "errors": n_err,
                     "mode": "all" if a.all else "ladder", "max_n": a.max_n,
                     "largest_measured_n": rows[-1][0] if rows else None,
                     "peak_gib_at_max": round(rows[-1][1] / GB, 3) if rows else None}
    return 0


if __name__ == "__main__":
    sys.exit(main())

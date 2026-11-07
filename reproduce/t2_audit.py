"""Look for suboptimal Chang & Yap optima, using T2 to say where.

    uv run python reproduce/t2_audit.py [--optima SOLVE.sqlite] [--out FRESH.sqlite]

T2 says the midpoint of every edge of a maximum-area Q lies in the convex hull of that
edge's contacts with dP. Where it fails, build the rotation it suggests and certify in
exact rationals that the result is convex, contained in P and strictly larger: a proof
of suboptimality that does not depend on T2.

Containment splits each candidate edge at its exact intersections with dP and tests
each piece's midpoint. The original Q is tested the same way first; if it fails, the
polygon is reported TESTFAIL and no claim is made.

Reads a stored solve_chang_and_yap result (default: the newest). analyse_t2.py runs
audit_one with provenance and stores the results.
"""

import argparse
import json
import os
import sqlite3
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyse_lemma4 import (cross, edge_data, merge_collinear,   # noqa: E402
                            reflex_flags, seg_overlap)

DELTAS = [Fr(1, 10 ** k) for k in range(2, 13)]


def clip(poly, a, b):
    """Clip a convex polygon by the closed half-plane LEFT of the directed line a->b."""
    out = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        sp, sq = cross(a, b, p), cross(a, b, q)
        if sp >= 0:
            out.append(p)
        if (sp > 0 > sq) or (sp < 0 < sq):
            t = sp / (sp - sq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def area(poly):
    n = len(poly)
    return abs(sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1]
                   for i in range(n))) / 2


def is_convex_ccw(poly):
    n = len(poly)
    return all(cross(poly[i - 1], poly[i], poly[(i + 1) % n]) >= 0 for i in range(n))


def point_in_poly(P, pt):
    """Exact membership in the CLOSED simple polygon P (boundary counts as inside)."""
    n = len(P)
    for j in range(n):
        u, w = P[j], P[(j + 1) % n]
        if cross(u, w, pt) == 0 and \
                min(u[0], w[0]) <= pt[0] <= max(u[0], w[0]) and \
                min(u[1], w[1]) <= pt[1] <= max(u[1], w[1]):
            return True
    inside = False
    for j in range(n):
        u, w = P[j], P[(j + 1) % n]
        if (u[1] > pt[1]) != (w[1] > pt[1]):
            if u[0] + (pt[1] - u[1]) * (w[0] - u[0]) / (w[1] - u[1]) > pt[0]:
                inside = not inside
    return inside


def seg_in_poly(P, A, B):
    """Exact: is the segment AB contained in the closed simple polygon P?"""
    ts = {Fr(0), Fr(1)}
    for j in range(len(P)):
        for lo, hi in seg_overlap(A, B, P[j], P[(j + 1) % len(P)]):
            ts.add(lo)
            ts.add(hi)
    ts = sorted(ts)
    for k in range(len(ts) - 1):
        t = (ts[k] + ts[k + 1]) / 2
        if not point_in_poly(P, (A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1]))):
            return False
    return True


def contained(P, K):
    return all(seg_in_poly(P, K[i], K[(i + 1) % len(K)]) for i in range(len(K)))


def rebuild(Q, ei, pivot, newdir):
    """Q with edge ei's line replaced by the line through `pivot` along `newdir`."""
    m = len(Q)
    big = Fr(10) ** 7
    cur = [(-big, -big), (big, -big), (big, big), (-big, big)]
    q2 = (pivot[0] + newdir[0], pivot[1] + newdir[1])
    for i in range(m):
        if i == ei:
            a, b, ref = pivot, q2, Q[(ei + 2) % m]
        else:
            a, b, ref = Q[i], Q[(i + 1) % m], Q[(i + 2) % m]
        cur = clip(cur, a, b) if cross(a, b, ref) > 0 else clip(cur, b, a)
        if len(cur) < 3:
            return None
    return cur


def improve(P, Q, ei, lo, hi):
    """Build the rotation T2 names for a T2-violating edge; return (gain, delta, poly).

    The midpoint t=1/2 escapes the contact hull [lo,hi].  Pivot strictly between the
    hull and the midpoint, and rotate so the midpoint side swings outward: every
    contact then moves inward (strictly, as the pivot separates them), and the area moves at
    |e| * v(midpoint) > 0.
    """
    m = len(Q)
    A, B = Q[ei], Q[(ei + 1) % m]
    d = (B[0] - A[0], B[1] - A[1])
    nout = (d[1], -d[0])
    if cross(A, B, Q[(ei + 2) % m]) < 0:
        nout = (-nout[0], -nout[1])
    half = Fr(1, 2)
    if hi < half:
        z, sgn = (hi + half) / 2, 1        # far side (t > z) swings out
    elif lo > half:
        z, sgn = (half + lo) / 2, -1       # near side (t < z) swings out
    else:
        return None
    pivot = (A[0] + z * d[0], A[1] + z * d[1])
    a0 = area(Q)
    best = None
    for delta in DELTAS:
        nd = (d[0] + sgn * delta * nout[0], d[1] + sgn * delta * nout[1])
        cand = rebuild(Q, ei, pivot, nd)
        if cand is None or len(cand) < 3:
            continue
        g = area(cand) - a0
        if g <= 0 or not is_convex_ccw(cand):
            continue
        if contained(P, cand):
            if best is None or g > best[0]:
                best = (g, delta, cand)
    return best


def audit_one(P, Q, expect_area=None):
    """Audit one (polygon, skull) pair. Returns a dict:
      status   "OK" | "TESTFAIL" | "DEGENERATE"
      t2_fail  edges whose midpoint escapes the hull of that edge's contacts
      beaten   of those, how many yielded a CERTIFIED larger convex contained body
      best_rel largest relative gain found
      hits     list of (edge_index, lo, hi, gain, delta, candidate_polygon)

    `expect_area`, when given, is the independently known exact area; a mismatch is
    a harness self-test failure and the pair is reported TESTFAIL with no claim made.
    """
    Q = merge_collinear(Q)
    if len(P) < 3 or len(Q) < 3:
        return {"status": "DEGENERATE", "t2_fail": 0, "beaten": 0,
                "best_rel": 0.0, "hits": [], "area": None}
    if not is_convex_ccw(Q):
        Q = Q[::-1]
    a0 = area(Q)
    if ((expect_area is not None and a0 != Fr(expect_area))
            or not is_convex_ccw(Q) or not contained(P, Q)):
        return {"status": "TESTFAIL", "t2_fail": 0, "beaten": 0,
                "best_rel": 0.0, "hits": [], "area": a0}
    rf = reflex_flags(P)
    m = len(Q)
    t2f = beat = unconfirmed = 0
    best_rel = 0.0
    hits = []
    for i in range(m):
        pts, spans, _c = edge_data(P, rf, Q[i], Q[(i + 1) % m])
        allc = pts + spans
        if not allc:
            continue
        lo = min(c[0] for c in allc)
        hi = max(c[1] for c in allc)
        if lo <= Fr(1, 2) <= hi:
            continue
        t2f += 1
        r = improve(P, Q, i, lo, hi)
        if r is None:
            unconfirmed += 1
            continue
        g, delta, cand = r
        beat += 1
        rel = float(g / a0)
        best_rel = max(best_rel, rel)
        hits.append((i, lo, hi, g, delta, cand))
    return {"status": "OK", "t2_fail": t2f, "beaten": beat, "best_rel": best_rel,
            "hits": hits, "area": a0, "unconfirmed": unconfirmed}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--optima", default=None,
                    help="stored solve_chang_and_yap result (default: newest)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--limit", type=int, default=10 ** 9)
    args = ap.parse_args()

    from _common import chang_and_yap_optima
    rows = chang_and_yap_optima(args.optima)
    out = None
    if args.out:
        if os.path.exists(args.out):
            print(f"refusing to overwrite {args.out}")
            return 2
        out = sqlite3.connect(args.out)
        out.execute("CREATE TABLE beaten (pid INT, name TEXT, nreflex INT, edge INT,"
                    " lo TEXT, hi TEXT, old_area TEXT, new_area TEXT, rel_gain REAL,"
                    " delta TEXT, new_poly TEXT)")
        out.execute("CREATE TABLE summary (pid INT PRIMARY KEY, name TEXT, nreflex INT,"
                    " edges INT, t2_fail INT, beaten INT, best_rel REAL, status TEXT)")

    npoly = nbeat = ntestfail = nt2 = nunconfirmed = 0
    worst = []
    for pid, name, nr, ae, P, Q in rows[:args.limit]:
        res = audit_one(P, Q, expect_area=ae)
        if res["status"] == "DEGENERATE":
            continue
        npoly += 1
        status = res["status"]
        if status == "TESTFAIL":
            ntestfail += 1
            if out:
                out.execute("INSERT INTO summary VALUES (?,?,?,?,?,?,?,?)",
                            (pid, name, nr, len(Q), 0, 0, 0.0, status))
            continue
        a0 = res["area"]
        t2f, beat, best_rel = res["t2_fail"], res["beaten"], res["best_rel"]
        nt2 += t2f
        nunconfirmed += res["unconfirmed"]
        m = len(merge_collinear(Q))
        for (i, lo, hi, g, delta, cand) in res["hits"]:
            if out:
                out.execute("INSERT INTO beaten VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                            (pid, name, nr, i, str(lo), str(hi), str(a0),
                             str(a0 + g), float(g / a0), str(delta),
                             json.dumps([[str(x), str(y)] for x, y in cand])))
        if beat:
            nbeat += 1
            worst.append((best_rel, pid, name, t2f, beat))
        if out:
            out.execute("INSERT INTO summary VALUES (?,?,?,?,?,?,?,?)",
                        (pid, name, nr, m, t2f, beat, best_rel, status))
            out.commit()
        if npoly % 25 == 0:
            print(f"  ... {npoly} polygons, {nt2} T2-failing edges, "
                  f"{nbeat} optima BEATEN", flush=True)

    print(f"\n=== T2 audit over {npoly} exact optima ===")
    print(f"  harness self-test failures (skipped)   {ntestfail}")
    print(f"  edges failing T2                       {nt2}")
    print(f"  of those, an improvement CERTIFIED     {nt2 - nunconfirmed}")
    print(f"  of those, no improvement found         {nunconfirmed}")
    print(f"  OPTIMA BEATEN (convex+contained+larger) {nbeat} / {npoly}")
    if worst:
        worst.sort(reverse=True)
        print("\n  largest relative gains (rel_gain, pid, name, t2_fail_edges, beaten):")
        for w in worst[:12]:
            print(f"    {w[0]:.3e}  id {w[1]:<5} {w[2]:<22} t2fail {w[3]}  beaten {w[4]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

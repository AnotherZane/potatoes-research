"""Exact structural tests on the edges of stored Chang & Yap optima.

    uv run python reproduce/analyse_lemma4.py [--optima SOLVE.sqlite]

Reads the optima of a stored solve_chang_and_yap result (default: the newest) as exact
rationals. Edges of Q are maximal flat pieces (collinear vertices merged). For each
edge e, with contact set C = closure(e) n dP and midpoint g:

  T1        C != {}
  T2        g in conv(C)
  T2r       g in conv(the reflex corners of P on e, plus collinear overlaps):
            Chang & Yap's "balanced" condition, 1986 p161 (2)
  sharpT2   g in conv(B), B = C minus endpoint contacts protected by an adjacent edge
            of Q that is flush with dP there

  S1        no convex corner of P lies in the relative interior of e
  S2        every isolated contact in the relative interior of e is a reflex corner

  Lemma 4   alphaS  e is contained in an edge of P
            alphaW  e overlaps an edge of P in positive length
            beta    the midpoint of e is a reflex corner of P
            gamma   e contains >= 2 corners of P, at least one reflex

Every predicate is an exact rational comparison.
"""

import argparse
import sys
from fractions import Fraction as Fr


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def reflex_flags(P):
    n = len(P)
    a2 = sum(P[i][0] * P[(i + 1) % n][1] - P[(i + 1) % n][0] * P[i][1] for i in range(n))
    ccw = a2 > 0
    return [(cross(P[i - 1], P[i], P[(i + 1) % n]) < 0) if ccw else
            (cross(P[i - 1], P[i], P[(i + 1) % n]) > 0) for i in range(n)]


def merge_collinear(Q):
    V = [p for i, p in enumerate(Q) if p != Q[i - 1]]
    out = [V[i] for i in range(len(V)) if cross(V[i - 1], V[i], V[(i + 1) % len(V)]) != 0]
    return out if len(out) >= 3 else V


def seg_overlap(A, B, C, D):
    """Exact [A,B] n [C,D] as t-intervals along AB; a point is (t, t)."""
    d1 = (B[0] - A[0], B[1] - A[1])
    d2 = (D[0] - C[0], D[1] - C[1])
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if den != 0:
        t = ((C[0] - A[0]) * d2[1] - (C[1] - A[1]) * d2[0]) / den
        u = ((C[0] - A[0]) * d1[1] - (C[1] - A[1]) * d1[0]) / den
        return [(t, t)] if 0 <= t <= 1 and 0 <= u <= 1 else []
    if cross(A, B, C) != 0:
        return []
    L = d1[0] * d1[0] + d1[1] * d1[1]
    if L == 0:
        return []
    tc = ((C[0] - A[0]) * d1[0] + (C[1] - A[1]) * d1[1]) / L
    td = ((D[0] - A[0]) * d1[0] + (D[1] - A[1]) * d1[1]) / L
    lo, hi = (tc, td) if tc <= td else (td, tc)
    lo, hi = max(lo, Fr(0)), min(hi, Fr(1))
    return [(lo, hi)] if lo <= hi else []


def edge_data(P, rf, A, B):
    """Contacts of segment AB with dP: isolated points and collinear spans."""
    n = len(P)
    pts, spans = [], []
    for j in range(n):
        for iv in seg_overlap(A, B, P[j], P[(j + 1) % n]):
            (spans if iv[0] != iv[1] else pts).append(iv)
    corners = []                       # corners of P lying on e
    for k, v in enumerate(P):
        if cross(A, B, v) == 0:
            d = (B[0] - A[0], B[1] - A[1])
            t = ((v[0] - A[0]) * d[0] + (v[1] - A[1]) * d[1]) / (d[0] ** 2 + d[1] ** 2)
            if 0 <= t <= 1:
                corners.append((t, k, rf[k]))
    return pts, spans, corners


def in_span(t, spans):
    return any(lo <= t <= hi for lo, hi in spans)


def halfplane_inside_cone(P, rf, y, n_out):
    """Is the half-plane {x : <n_out, x-y> <= 0} contained in the tangent cone T_y P?

    Chang & Yap's `D n C+ subset P` (1986 p158) as a cone condition. It decides
    protection: a contact y on the boundary of edge F cannot block F if another edge at
    y already satisfies it.

    A half-plane fits only where:
      * y in the relative interior of an edge of P -- the half-plane must BE that
        edge's inward half-plane;
      * y a REFLEX corner -- the exterior wedge cone(r1,r2) is convex, so the
        half-plane fits iff <n_out, r1> >= 0 and <n_out, r2> >= 0.
    A convex corner admits no half-plane at all.
    """
    n = len(P)
    for k, v in enumerate(P):
        if v == y:
            if not rf[k]:
                return False                      # convex corner: no half-plane fits
            r1 = (P[k - 1][0] - v[0], P[k - 1][1] - v[1])
            r2 = (P[(k + 1) % n][0] - v[0], P[(k + 1) % n][1] - v[1])
            return (n_out[0] * r1[0] + n_out[1] * r1[1] >= 0
                    and n_out[0] * r2[0] + n_out[1] * r2[1] >= 0)
    for j in range(n):
        u, w = P[j], P[(j + 1) % n]
        if cross(u, w, y) == 0:
            d = (w[0] - u[0], w[1] - u[1])
            L = d[0] * d[0] + d[1] * d[1]
            t = ((y[0] - u[0]) * d[0] + (y[1] - u[1]) * d[1]) / L
            if 0 < t < 1:
                # T_y P is this edge's inward half-plane; the only half-plane it
                # contains is itself, i.e. n_out must be its outward normal.
                return n_out[0] * d[0] + n_out[1] * d[1] == 0 and (
                    n_out[0] * (P[j - 1][0] - y[0])
                    + n_out[1] * (P[j - 1][1] - y[1])) <= 0
    return False


def outward_normal(A, B, interior_pt):
    """Outward normal of the half-plane bounded by line AB that holds interior_pt."""
    d = (B[0] - A[0], B[1] - A[1])
    n_out = (d[1], -d[0])
    s = n_out[0] * (interior_pt[0] - A[0]) + n_out[1] * (interior_pt[1] - A[1])
    return n_out if s <= 0 else (-n_out[0], -n_out[1])


def analyse(P, Q):
    rf = reflex_flags(P)
    m = len(Q)
    # any interior point of the convex Q, used only to orient half-planes
    ctr = (sum(p[0] for p in Q) / m, sum(p[1] for p in Q) / m)
    per = []
    raw = [edge_data(P, rf, Q[i], Q[(i + 1) % m]) for i in range(m)]
    for i in range(m):
        A, B = Q[i], Q[(i + 1) % m]
        pts, spans, corners = raw[i]
        allc = pts + spans
        half = Fr(1, 2)
        t1 = bool(allc)
        lo = min(c[0] for c in allc) if allc else None
        hi = max(c[1] for c in allc) if allc else None
        t2 = t1 and lo <= half <= hi

        rc = [t for t, k, r in corners if r]
        rlo = min(rc + [s[0] for s in spans]) if (rc or spans) else None
        rhi = max(rc + [s[1] for s in spans]) if (rc or spans) else None
        t2r = rlo is not None and rlo <= half <= rhi

        # --- blocking set -------------------------------------------------
        # An endpoint contact is NOT blocking when the adjacent edge of Q is
        # flush with dP there: the neighbour already lies in that face of P, so
        # material added past e at that point stays on P's interior side.
        prev_spans = raw[(i - 1) % m][1]        # edge ending at A
        next_spans = raw[(i + 1) % m][1]        # edge starting at B
        # conservative rule: the neighbour is flush with an edge of P at the shared vertex.
        prot_A = in_span(Fr(1), prev_spans)     # A is the far end of the previous edge
        prot_B = in_span(Fr(0), next_spans)     # B is the near end of the next edge
        # sharp rule: the neighbour's half-plane fits inside T_y P at the shared
        # vertex.  Strictly weaker blocking set, so a strictly stronger condition.
        prot_A2 = halfplane_inside_cone(P, rf, A, outward_normal(Q[(i - 1) % m], A, ctr))
        prot_B2 = halfplane_inside_cone(P, rf, B, outward_normal(B, Q[(i + 2) % m], ctr))
        def blocking(pA, pB):
            out = []
            for a, b in allc:
                # a SPAN touching an endpoint still blocks: e itself lies in dP there
                if a == b and a == 0 and pA:
                    continue
                if a == b and a == 1 and pB:
                    continue
                out.append((a, b))
            return out

        B_set = blocking(prot_A, prot_B)
        blo = min(c[0] for c in B_set) if B_set else None
        bhi = max(c[1] for c in B_set) if B_set else None
        sharp = blo is not None and blo <= half <= bhi

        B2 = blocking(prot_A2, prot_B2)
        b2lo = min(c[0] for c in B2) if B2 else None
        b2hi = max(c[1] for c in B2) if B2 else None
        sharp2 = b2lo is not None and b2lo <= half <= b2hi

        # --- structural claims --------------------------------------------
        # S1: a corner of P strictly inside e, not swallowed by a collinear
        # overlap, must be REFLEX (a line through a convex corner leaves P).
        s1 = all(r for t, k, r in corners if 0 < t < 1 and not in_span(t, spans))
        # S2: every isolated contact in relint(e) must BE such a reflex corner --
        # there is no third way to touch a facet's interior.
        isolated = [t for t, t_ in pts if 0 < t < 1 and not in_span(t, spans)]
        s2 = all(any(t == tc and r for tc, k, r in corners) for t in isolated)
        # The configuration Lemma 4 would forbid and T2 alone permits: the START
        # vertex of e lies on dP, is NOT a corner of P, and NEITHER incident edge
        # of Q is flush with dP there -- a "bare node" pinned at a point in the
        # relative interior of an edge of P.  Counted, not asserted.
        bare_node = (any(a <= 0 <= b for a, b in allc)
                     and not any(tc == 0 for tc, k, r in corners)
                     and not in_span(Fr(0), spans)
                     and not in_span(Fr(1), prev_spans))

        # --- Lemma 4 -------------------------------------------------------
        alphaS = any(lo_ == 0 and hi_ == 1 for lo_, hi_ in spans)
        alphaW = bool(spans)
        mid = ((A[0] + B[0]) / 2, (A[1] + B[1]) / 2)
        beta = any(v == mid and rf[k] for k, v in enumerate(P))
        gamma = len(corners) >= 2 and any(r for _t, _k, r in corners)

        per.append(dict(t1=t1, t2=t2, t2r=t2r, sharp=sharp, sharp2=sharp2, nblock2=len(B2), s1=s1, s2=s2,
                        alphaS=alphaS, alphaW=alphaW, beta=beta, gamma=gamma,
                        ncorner=len(corners), nreflex=len(rc), nspan=len(spans),
                        nblock=len(B_set), lo=lo, hi=hi,
                        bare_node=bare_node, nisolated=len(isolated),
                        pure_beta=(len(isolated) == 1 and not spans
                                   and len(B_set) == 1)))
    return per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--optima", default=None,
                    help="stored solve_chang_and_yap result (default: newest)")
    ap.add_argument("--show", type=int, default=6)
    args = ap.parse_args()

    from _common import chang_and_yap_optima
    rows = chang_and_yap_optima(args.optima)
    tot = 0
    fails = {k: 0 for k in ('t1', 't2', 't2r', 'sharp', 'sharp2', 's1', 's2', 'L4', 'L4S')}
    census = {}
    examples = {k: [] for k in fails}
    npoly = 0
    extra = {'bare_node': 0, 'pure_beta': 0, 'relint_only': 0}
    for pid, name, nr, _, P, Q in rows:
        Q = merge_collinear(Q)
        if len(P) < 3 or len(Q) < 3:
            continue
        npoly += 1
        for i, r in enumerate(analyse(P, Q)):
            tot += 1
            l4 = r['alphaW'] or r['beta'] or r['gamma']
            l4s = r['alphaS'] or r['beta'] or r['gamma']
            for k, v in (('t1', r['t1']), ('t2', r['t2']), ('t2r', r['t2r']),
                         ('sharp', r['sharp']), ('sharp2', r['sharp2']),
                         ('s1', r['s1']), ('s2', r['s2']),
                         ('L4', l4), ('L4S', l4s)):
                if not v:
                    fails[k] += 1
                    if len(examples[k]) < args.show:
                        examples[k].append((pid, name, i, r['ncorner'],
                                            r['nreflex'], r['nspan'],
                                            str(r['lo']), str(r['hi'])))
            key = ''.join(c for c, v in (('A', r['alphaS']), ('a', r['alphaW'] and not r['alphaS']),
                                         ('b', r['beta']), ('g', r['gamma'])) if v) or 'NONE'
            census[key] = census.get(key, 0) + 1
            extra['bare_node'] += r['bare_node']
            extra['pure_beta'] += r['pure_beta']
            extra['relint_only'] += (r['lo'] is not None and r['lo'] > 0 and r['hi'] < 1)

    print(f"\n=== EXACT structural analysis: {npoly} optima, {tot} edges ===\n")
    label = {
        't1': "T1     e meets dP",
        't2': "T2     midpoint in conv(all contacts)",
        't2r': "T2r    midpoint in conv(reflex contacts)   [C&Y 'balanced', p161(2)]",
        'sharp': "sharpT2 midpoint in conv(BLOCKING), flush-protection rule",
        'sharp2': "sharpT2* midpoint in conv(BLOCKING), tangent-cone protection rule",
        's1': "S1     no convex corner of P inside e",
        's2': "S2     every isolated interior contact is a reflex corner",
        'L4': "L4     alpha(overlap) or beta or gamma",
        'L4S': "L4S    alpha(e inside an edge of P) or beta or gamma",
    }
    for k in ('t1', 't2', 't2r', 'sharp', 'sharp2', 's1', 's2', 'L4', 'L4S'):
        print(f"  {label[k]:<66} {fails[k]:5d} failures")
    print(f"\n  census (A=e in edge of P, a=partial overlap, b=beta, g=gamma): "
          f"{dict(sorted(census.items()))}")
    print(f"  edges whose contacts are all in the relative INTERIOR   "
          f"{extra['relint_only']}")
    print(f"  pure (beta) edges: one isolated reflex contact, blocking  "
          f"{extra['pure_beta']}")
    print(f"  BARE NODES: vertex of Q on dP, not a corner, neither edge flush "
          f"{extra['bare_node']}   (Lemma 4 forbids; T2 alone permits)")
    for k in fails:
        if examples[k]:
            print(f"\n  {k} examples (pid,name,edge,ncorner,nreflex,nspan,lo,hi):")
            for e in examples[k]:
                print("     ", e)
    return 0


if __name__ == "__main__":
    sys.exit(main())

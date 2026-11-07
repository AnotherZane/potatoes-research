"""A certified upper bound on A*(P).

PIVOT BOUND. Let Q* be an optimal convex subset and u a point of dP (in practice a
reflex corner). u is not interior to Q* (it is not interior to P), so a line through u
has Q* in one closed half-plane H+, and

    A*  <=  U(u) := max over theta of area(P ∩ H+(u, theta)),
    A*  <=  min over u of U(u).                                              (*)

COMPUTING U(u). As the line rotates, the area of P ∩ H+ is Lipschitz in theta:

    d/dtheta area(P ∩ H+)  =  (1/2) [ m-(theta) - m+(theta) ]

where m+ is the sum of (r2^2 - r1^2) over the maximal intervals [r1,r2] in which
the FORWARD ray from u lies inside P, and m- the same for the backward ray.
Both are bounded by R_u^2 with R_u = max distance from u to P, so

    |d/dtheta area|  <=  L := R_u^2 / 2.

So on [t1, t2] with sampled areas a1, a2 the maximum is at most
(a1 + a2 + L*(t2-t1)) / 2, and branch and bound on theta converges to U(u) from
above; every intermediate value is already a valid bound. A pivot whose best sample
exceeds the best bound so far cannot be the minimiser in (*) and is abandoned.

Areas are exact (Fraction, Sutherland-Hodgman); only angles are floats, and the
angular gap in the slack is an exact difference plus a margin, so it is never
understated.

WIDTH BOUND: see width_bound_dir. The result is the smaller of the two.

`--check-clip` tests the clipper against shapely. `--verify` checks the bound
independently: each pivot's bound must dominate a dense shapely sweep, and the
width bound a numerical integral.

    uv run python reproduce/upper_bound.py --suite
    uv run python reproduce/upper_bound.py --sawtooth
    uv run python reproduce/upper_bound.py --suite --verify
    uv run python reproduce/upper_bound.py --ids 499 210 314
    uv run python reproduce/upper_bound.py --optima --limit 60
"""

import argparse
import heapq
import math
import sqlite3
import sys
import time
from fractions import Fraction

from shapely import wkb as shapely_wkb

from potato_peeling.chang_and_yap import pre_checks

ANG_EPS = Fraction(1, 10 ** 12)   # see u_pivot: the slack must never be understated


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #
def shoelace(pts):
    """Signed area of a closed path, exact."""
    n = len(pts)
    if n < 3:
        return Fraction(0)
    s = Fraction(0)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return s / 2


def clip_halfplane(pts, ux, uy, dx, dy):
    """Sutherland-Hodgman clip of a closed path against the closed half-plane left of
    the directed line (u, d). For a concave subject the extra edges lie on the clip
    line and add no area, so the shoelace of the output is the clipped area."""
    def side(p):
        # cross(d, p - u) > 0  <=>  p strictly left of the directed line
        return dx * (p[1] - uy) - dy * (p[0] - ux)

    out = []
    n = len(pts)
    for i in range(n):
        a = pts[i]
        b = pts[(i + 1) % n]
        sa = side(a)
        sb = side(b)
        if sa >= 0:
            out.append(a)
        if (sa > 0 and sb < 0) or (sa < 0 and sb > 0):
            t = sa / (sa - sb)
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return out


def halfplane_area(pts, u, d):
    """Exact area of P ∩ H+(u, d)."""
    return abs(shoelace(clip_halfplane(pts, u[0], u[1], d[0], d[1])))


def reflex_indices(pts):
    """Reflex corners of a CCW polygon: a right turn."""
    n = len(pts)
    out = []
    for i in range(n):
        ax, ay = pts[(i - 1) % n]
        bx, by = pts[i]
        cx, cy = pts[(i + 1) % n]
        if (bx - ax) * (cy - by) - (by - ay) * (cx - bx) < 0:
            out.append(i)
    return out


def max_radius_sq(pts, u):
    """max over x in P of |x - u|^2, attained at a vertex."""
    return max((x - u[0]) ** 2 + (y - u[1]) ** 2 for x, y in pts)


def convex_hull(pts):
    """Monotone chain, exact."""
    p = sorted(set(pts))
    if len(p) < 3:
        return p

    def half(seq):
        out = []
        for q in seq:
            while len(out) >= 2:
                (x1, y1), (x2, y2) = out[-2], out[-1]
                if (x2 - x1) * (q[1] - y1) - (y2 - y1) * (q[0] - x1) > 0:
                    break
                out.pop()
            out.append(q)
        return out[:-1]

    return half(p) + half(reversed(p))


class RadialBound:
    """Per-sector upper bound on |x - u|^2 over P, from P's convex hull: the farthest
    hull point in a sector is a hull vertex in it or on one of its bounding rays. A
    global radius is far too loose for elongated polygons."""

    def __init__(self, pts, u, nsec=256):
        self.nsec = nsec
        self.u = u
        self.glob = float(max_radius_sq(pts, u))
        hull = convex_hull(pts)
        # hull half-planes n.x <= c, outward normals (hull is CCW)
        hp = []
        for i in range(len(hull)):
            a, b = hull[i], hull[(i + 1) % len(hull)]
            nx, ny = b[1] - a[1], a[0] - b[0]
            hp.append((nx, ny, nx * a[0] + ny * a[1]))

        def ray_exit_sq(th):
            dx, dy = math.cos(th), math.sin(th)
            best = None
            for (nx, ny, c) in hp:
                den = float(nx) * dx + float(ny) * dy
                if den <= 0:
                    continue
                t = (float(c) - float(nx) * float(u[0]) - float(ny) * float(u[1])) / den
                if t >= 0 and (best is None or t < best):
                    best = t
            # Cannot happen for a bounded hull; fall back to the global radius,
            # never 0, which would understate the bound.
            return self.glob if best is None else best * best

        # per-sector: hull vertices whose direction falls in it, plus both edges
        sec = [0.0] * nsec
        for v in hull:
            dx, dy = float(v[0] - u[0]), float(v[1] - u[1])
            if dx == 0 and dy == 0:
                continue
            k = int((math.atan2(dy, dx) % (2 * math.pi)) / (2 * math.pi) * nsec) % nsec
            sec[k] = max(sec[k], dx * dx + dy * dy)
        for k in range(nsec):
            lo = 2 * math.pi * k / nsec
            hi = 2 * math.pi * (k + 1) / nsec
            sec[k] = max(sec[k], ray_exit_sq(lo), ray_exit_sq(hi))
        self.sec = sec

    def over(self, a1, a2):
        """Upper bound on the squared radius for directions in [a1,a2] and its
        opposite sector (the backward ray is what m- measures)."""
        span = a2 - a1
        if span >= math.pi:
            return max(self.sec)
        best = 0.0
        for base in (a1, a1 + math.pi):
            k0 = int((base % (2 * math.pi)) / (2 * math.pi) * self.nsec)
            k1 = int(((base + span) % (2 * math.pi)) / (2 * math.pi) * self.nsec)
            k = k0
            while True:
                best = max(best, self.sec[k % self.nsec])
                if k % self.nsec == k1 % self.nsec:
                    break
                k += 1
                if k - k0 > self.nsec:
                    break
        return best


# --------------------------------------------------------------------------- #
# U(u): branch and bound over the rotation angle
# --------------------------------------------------------------------------- #
def u_pivot(pts, u, tol_abs, cutoff=None, max_evals=200000):
    """Certified max over theta of area(P ∩ H+(u, theta)).

    Returns (bound, best_sample, evals). `bound` is always a valid upper bound on
    the true maximum, whatever the search did.

    `cutoff`: abandon as soon as best_sample exceeds it -- this pivot cannot be
    the minimiser in (*), so its exact value is irrelevant.
    """
    L = Fraction(max_radius_sq(pts, u)) / 2      # Lipschitz constant, exact

    def area_at(theta):
        return halfplane_area(pts, u, (math.cos(theta), math.sin(theta)))

    # Seed with the directions that matter combinatorially: the O(n) angles at
    # which the line passes through a vertex of P, both orientations, plus a
    # uniform floor so no wide gap is left uncovered.
    seeds = set()
    for v in pts:
        if v[0] == u[0] and v[1] == u[1]:
            continue
        a = math.atan2(float(v[1] - u[1]), float(v[0] - u[0])) % (2 * math.pi)
        seeds.add(a)
        seeds.add((a + math.pi) % (2 * math.pi))
    for k in range(32):
        seeds.add(2 * math.pi * k / 32)
    grid = sorted(seeds)
    grid.append(grid[0] + 2 * math.pi)           # close the circle

    vals = [area_at(t) for t in grid[:-1]]
    vals.append(vals[0])
    evals = len(grid) - 1
    best = max(vals)

    # The margin covers the sample at nominal angle t being taken at the angle of
    # the rounded (cos t, sin t), about 1e-16 away.
    def sawtooth(alo, ahi, lo, hi):
        gap = Fraction(hi) - Fraction(lo) + ANG_EPS
        return (alo + ahi + L * gap) / 2

    # priority queue of intervals, keyed by the sawtooth bound (negated for heapq)
    heap = []
    for i in range(len(grid) - 1):
        b = sawtooth(vals[i], vals[i + 1], grid[i], grid[i + 1])
        heapq.heappush(heap, (-b, grid[i], grid[i + 1], vals[i], vals[i + 1]))

    while heap:
        nb, t1, t2, a1, a2 = heap[0]
        bound = -nb
        if bound <= best + tol_abs:
            break
        if cutoff is not None and best > cutoff:
            break                                 # cannot be the minimiser
        if evals >= max_evals:
            break
        heapq.heappop(heap)
        tm = (t1 + t2) / 2
        am = area_at(tm)
        evals += 1
        if am > best:
            best = am
        for (lo, hi, alo, ahi) in ((t1, tm, a1, am), (tm, t2, am, a2)):
            heapq.heappush(heap, (-sawtooth(alo, ahi, lo, hi), lo, hi, alo, ahi))

    bound = max(best, -heap[0][0]) if heap else best
    return bound, best, evals


# --------------------------------------------------------------------------- #
# U(u,v): the two-pivot bound
# --------------------------------------------------------------------------- #
# The pivot argument applies to u and v at once, so
#
#     A*  <=  U(u,v) := max over (a,b) of area(P ∩ H+(u,a) ∩ H+(v,b))  <=  min(U(u), U(v)).
#
# f is Lipschitz in each angle with the single-pivot constants, so on a box of side
# (h1,h2) with opposite corner values f_c and f_c',
#
#     f(x) <= f_c  + Lu*(x1-c1)      + Lv*(x2-c2)
#     f(x) <= f_c' + Lu*(h1-(x1-c1)) + Lv*(h2-(x2-c2))
#  => f(x) <= ( f_c + f_c' + Lu*h1 + Lv*h2 ) / 2
#
# taken as the min over both diagonals. The worst box is split along the axis with
# more slack (2 new evaluations per split).
def halfplane_area2(pts, u, du, v, dv):
    """Exact area of P ∩ H+(u,du) ∩ H+(v,dv), by clipping against each in turn."""
    path = clip_halfplane(pts, u[0], u[1], du[0], du[1])
    if len(path) < 3:
        return Fraction(0)
    path = clip_halfplane(path, v[0], v[1], dv[0], dv[1])
    if len(path) < 3:
        return Fraction(0)
    return abs(shoelace(path))


def uv_pair(pts, u, v, tol_abs, cutoff=None, init=12, max_evals=4000):
    """Certified max over (a,b) of area(P ∩ H+(u,a) ∩ H+(v,b)), as
    (bound, best_sample, evals). Stops early once the best sample exceeds `cutoff`."""
    # Per-box Lipschitz constants from the radius over the box's angle range.
    RBu, RBv = RadialBound(pts, u), RadialBound(pts, v)
    Lu_glob = Fraction(max_radius_sq(pts, u)) / 2
    Lv_glob = Fraction(max_radius_sq(pts, v)) / 2

    def Lu_of(a1, a2):
        return min(Lu_glob, Fraction(RBu.over(a1, a2) * (1 + 1e-9)) / 2 + ANG_EPS)

    def Lv_of(b1, b2):
        return min(Lv_glob, Fraction(RBv.over(b1, b2) * (1 + 1e-9)) / 2 + ANG_EPS)

    TAU = 2 * math.pi

    cache = {}

    def f(a, b):
        key = (a, b)
        hit = cache.get(key)
        if hit is None:
            hit = halfplane_area2(pts, u, (math.cos(a), math.sin(a)),
                                  v, (math.cos(b), math.sin(b)))
            cache[key] = hit
        return hit

    def box_bound(a1, a2, b1, b2, f11, f12, f21, f22):
        slack = Lu_of(a1, a2) * (Fraction(a2) - Fraction(a1) + ANG_EPS) \
            + Lv_of(b1, b2) * (Fraction(b2) - Fraction(b1) + ANG_EPS)
        return min(f11 + f22, f12 + f21) / 2 + slack / 2

    grid = [TAU * k / init for k in range(init + 1)]
    heap = []
    best = Fraction(0)
    evals = 0
    for i in range(init):
        for j in range(init):
            a1, a2, b1, b2 = grid[i], grid[i + 1], grid[j], grid[j + 1]
            vals = (f(a1, b1), f(a1, b2), f(a2, b1), f(a2, b2))
            best = max(best, *vals)
            heapq.heappush(heap, (-box_bound(a1, a2, b1, b2, *vals),
                                  a1, a2, b1, b2, vals))
    evals = len(cache)

    while heap:
        nb, a1, a2, b1, b2, vals = heap[0]
        if -nb <= best + tol_abs:
            break
        if cutoff is not None and best > cutoff:
            break
        if evals >= max_evals:
            break
        heapq.heappop(heap)
        f11, f12, f21, f22 = vals
        # split along the axis carrying more slack
        if (Lu_of(a1, a2) * (Fraction(a2) - Fraction(a1))
                >= Lv_of(b1, b2) * (Fraction(b2) - Fraction(b1))):
            am = (a1 + a2) / 2
            m1, m2 = f(am, b1), f(am, b2)
            best = max(best, m1, m2)
            heapq.heappush(heap, (-box_bound(a1, am, b1, b2, f11, f12, m1, m2),
                                  a1, am, b1, b2, (f11, f12, m1, m2)))
            heapq.heappush(heap, (-box_bound(am, a2, b1, b2, m1, m2, f21, f22),
                                  am, a2, b1, b2, (m1, m2, f21, f22)))
        else:
            bm = (b1 + b2) / 2
            m1, m2 = f(a1, bm), f(a2, bm)
            best = max(best, m1, m2)
            heapq.heappush(heap, (-box_bound(a1, a2, b1, bm, f11, m1, f21, m2),
                                  a1, a2, b1, bm, (f11, m1, f21, m2)))
            heapq.heappush(heap, (-box_bound(a1, a2, bm, b2, m1, f12, m2, f22),
                                  a1, a2, bm, b2, (m1, f12, m2, f22)))
        evals = len(cache)

    bound = max(best, -heap[0][0]) if heap else best
    return bound, best, evals


def two_pivot_bound(pts, single_bound, single_detail, pivots, top_m=6,
                    rel_tol=1e-3, verbose=False):
    """min of U(u,v) over pairs of the `top_m` most promising pivots (any subset of
    pairs gives a valid bound)."""
    if len(pivots) < 2:
        return single_bound
    area_P = abs(shoelace(pts))
    tol_abs = area_P * Fraction(rel_tol).limit_denominator(10 ** 9)

    # rank by best sample: with the cutoff, most pivots' bounds are loose
    order = sorted(range(len(pivots)), key=lambda i: single_detail[i][2])
    cand = order[:max(2, top_m)]

    best_bound = single_bound
    for x in range(len(cand)):
        for y in range(x + 1, len(cand)):
            u, v = pivots[cand[x]], pivots[cand[y]]
            b, s, ev = uv_pair(pts, u, v, tol_abs, cutoff=best_bound)
            if b < best_bound:
                best_bound = b
                if verbose:
                    print(f"    pair ({cand[x]},{cand[y]}) -> {float(b):.6f} "
                          f"(evals {ev})", flush=True)
    return best_bound


# --------------------------------------------------------------------------- #
# W(theta): the slice / width bound
# --------------------------------------------------------------------------- #
# Slice P by lines perpendicular to a direction w. Q*'s slice at each level is one
# interval inside P, hence inside one maximal sub-segment of P's slice:
#
#     len(Q* ∩ line_t)  <=  maxseg(t) := longest maximal sub-segment of P ∩ line_t
#
# and by Cavalieri, for EVERY direction w,
#
#     A*  <=  W(w) := integral of maxseg(t) dt.                              (**)
#
# (*) is tight where a reflex corner pinches P, (**) where a direction splits it.
#
# Between consecutive vertex projections the crossed edges are fixed, so each
# sub-segment length is linear in t and maxseg is convex on the cell, hence below
# its chord:
#
#     integral over the cell  <=  (maxseg(a) + maxseg(b)) * (b - a) / 2
#
# with the endpoint values from the cell's own linear pieces. Exact throughout.
def width_bound_dir(pts, w):
    """Exact upper bound on the integral of maxseg for a rational direction w. Works
    in the frame (x . v, x . w), v = w rotated 90 degrees, dividing by |w|^2 at the
    end, so nothing is irrational."""
    p, q = w
    vx, vy = -q, p
    W2 = p * p + q * q
    n = len(pts)
    T = [x * p + y * q for (x, y) in pts]
    events = sorted(set(T))
    if len(events) < 2:
        return Fraction(0)

    total = Fraction(0)
    for a, b in zip(events, events[1:]):
        tm = (a + b) / 2
        lin = []                       # (slope, intercept) of S(t) per crossed edge
        for i in range(n):
            j = (i + 1) % n
            if (T[i] < tm) == (T[j] < tm):
                continue               # edge does not straddle this level
            ex, ey = pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]
            ew = ex * p + ey * q       # nonzero: the edge straddles
            ev = ex * vx + ey * vy
            slope = ev / ew
            lin.append((slope, pts[i][0] * vx + pts[i][1] * vy - slope * T[i]))
        if len(lin) % 2:
            # The count must be even (tm avoids vertices); skipping the cell would
            # under-count, so fail instead.
            raise AssertionError(
                f"odd crossing count {len(lin)} at t={float(tm)}; the slice logic is wrong")
        if len(lin) == 0:
            continue                   # slice genuinely empty: contributes nothing
        lin.sort(key=lambda sc: sc[0] * tm + sc[1])
        # inside intervals alternate: (0,1), (2,3), ...
        segs = [(lin[2 * k + 1][0] - lin[2 * k][0], lin[2 * k + 1][1] - lin[2 * k][1])
                for k in range(len(lin) // 2)]
        fa = max(s * a + c for s, c in segs)
        fb = max(s * b + c for s, c in segs)
        total += (fa + fb) * (b - a) / 2
    return total / W2


def _candidate_directions(pts, n_uniform=180):
    """Edge directions, edge normals, and a uniform sweep -- all as exact rationals."""
    dirs = set()
    n = len(pts)
    for i in range(n):
        j = (i + 1) % n
        ex, ey = pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]
        if ex or ey:
            dirs.add((ex, ey))         # slice lines perpendicular to the edge
            dirs.add((-ey, ex))        # slice lines parallel to the edge
    for k in range(n_uniform):
        th = math.pi * k / n_uniform
        dirs.add((Fraction(round(math.cos(th) * 10 ** 7), 10 ** 7),
                  Fraction(round(math.sin(th) * 10 ** 7), 10 ** 7)))
    return [d for d in dirs if d[0] or d[1]]


def width_bound(pts, n_uniform=180, exact_top=8):
    """min over directions of W(w): rank in floats, evaluate the best few exactly
    (any direction is valid)."""
    cands = _candidate_directions(pts, n_uniform)
    fpts = [(float(x), float(y)) for x, y in pts]
    scored = []
    for w in cands:
        fw = (float(w[0]), float(w[1]))
        scored.append((_width_float(fpts, fw), w))
    scored.sort(key=lambda r: r[0])
    best = None
    for _, w in scored[:exact_top]:
        val = width_bound_dir(pts, w)
        if best is None or val < best:
            best = val
    return best


def _width_float(fpts, w):
    """Float twin of width_bound_dir, used only to RANK directions."""
    p, q = w
    vx, vy = -q, p
    W2 = p * p + q * q
    n = len(fpts)
    T = [x * p + y * q for (x, y) in fpts]
    events = sorted(set(T))
    if len(events) < 2:
        return float("inf")
    total = 0.0
    for a, b in zip(events, events[1:]):
        tm = (a + b) / 2
        lin = []
        for i in range(n):
            j = (i + 1) % n
            if (T[i] < tm) == (T[j] < tm):
                continue
            ex, ey = fpts[j][0] - fpts[i][0], fpts[j][1] - fpts[i][1]
            ew = ex * p + ey * q
            if ew == 0.0:
                continue
            slope = (ex * vx + ey * vy) / ew
            lin.append((slope, fpts[i][0] * vx + fpts[i][1] * vy - slope * T[i]))
        if len(lin) < 2 or len(lin) % 2:
            continue
        lin.sort(key=lambda sc: sc[0] * tm + sc[1])
        segs = [(lin[2 * k + 1][0] - lin[2 * k][0], lin[2 * k + 1][1] - lin[2 * k][1])
                for k in range(len(lin) // 2)]
        total += (max(s * a + c for s, c in segs)
                  + max(s * b + c for s, c in segs)) * (b - a) / 2
    return total / W2


def pivot_points(pts, edge_samples=0):
    """Reflex corners, plus `edge_samples` points per edge. Any point of dP is a
    valid pivot; convex corners give only area(P)."""
    out = [pts[i] for i in reflex_indices(pts)]
    if edge_samples:
        n = len(pts)
        for i in range(n):
            a, b = pts[i], pts[(i + 1) % n]
            for k in range(1, edge_samples + 1):
                lam = Fraction(k, edge_samples + 1)
                out.append((a[0] + lam * (b[0] - a[0]), a[1] + lam * (b[1] - a[1])))
    return out


def upper_bound(coords, rel_tol=1e-3, verbose=False, edge_samples=0, pairs=0):
    """(min over pivots of U(p), per-pivot detail); `pairs` > 0 adds the two-pivot
    bound over that many pivots."""
    pts = [(Fraction(x), Fraction(y)) for x, y in coords]
    area_P = abs(shoelace(pts))
    if not reflex_indices(pts):
        return area_P, []                        # convex: A* = area(P), exactly

    pivots = pivot_points(pts, edge_samples)
    tol_abs = area_P * Fraction(rel_tol).limit_denominator(10 ** 9)
    best_bound = area_P                          # trivially valid
    detail = []
    for idx, u in enumerate(pivots):
        b, s, ev = u_pivot(pts, u, tol_abs, cutoff=best_bound)
        detail.append((idx, b, s, ev))
        if b < best_bound:
            best_bound = b
        if verbose:
            print(f"    pivot {idx:<4} U<= {float(b):14.4f}  best sample {float(s):14.4f}"
                  f"  evals {ev}", flush=True)

    if pairs:
        best_bound = two_pivot_bound(pts, best_bound, detail, pivots, top_m=pairs,
                                     rel_tol=rel_tol, verbose=verbose)
    return best_bound, detail


# --------------------------------------------------------------------------- #
# inputs
# --------------------------------------------------------------------------- #
def sawtooth(t):
    n = 2 * t
    top = [(60.0 * i / n, 1.0 if i % 2 == 0 else 0.75) for i in range(n + 1)]
    return [(0.0, 0.0), (60.0, 0.0)] + top[::-1]


def load_corpus(ids=None, limit=None):
    from _common import require_corpus
    conn = sqlite3.connect(require_corpus())
    if ids:
        q = ("SELECT id,name,geometry_wkb FROM polygons WHERE id IN (%s)"
             % ",".join("?" * len(ids)))
        rows = conn.execute(q, ids).fetchall()
    else:
        rows = conn.execute("SELECT id,name,geometry_wkb FROM polygons").fetchall()
    conn.close()
    out = []
    for pid, name, blob in rows:
        g = shapely_wkb.loads(blob)
        if g.geom_type == "MultiPolygon":
            g = max(g.geoms, key=lambda x: x.area)
        c = list(g.exterior.coords)[:-1]
        if len(c) >= 3:
            out.append((pid, name, c))
    out.sort(key=lambda r: len(r[2]))
    return out[:limit] if limit else out


# --------------------------------------------------------------------------- #
# independent checks
# --------------------------------------------------------------------------- #
def verify(coords, dense=3000, edge_samples=0):
    """Check the bound without this file's geometry: each pivot's bound (and a few
    pair bounds) must dominate a dense shapely sweep, and the width bound a numerical
    integral. Returns (ok, messages)."""
    from shapely.geometry import Polygon as ShPoly
    pts = [(Fraction(x), Fraction(y)) for x, y in coords]
    P = ShPoly([(float(x), float(y)) for x, y in pts]).buffer(0)
    big = 10 * max(abs(float(v)) for p in pts for v in p) + 1000
    area_P = abs(shoelace(pts))
    msgs, ok = [], True

    def sh_quad(u, th):
        """A shapely half-plane, built independently of anything in this file."""
        d = (math.cos(th), math.sin(th))
        nx, ny = -d[1], d[0]
        ux, uy = float(u[0]), float(u[1])
        return ShPoly([(ux + d[0] * big, uy + d[1] * big),
                       (ux - d[0] * big, uy - d[1] * big),
                       (ux - d[0] * big + nx * big, uy - d[1] * big + ny * big),
                       (ux + d[0] * big + nx * big, uy + d[1] * big + ny * big)]).buffer(0)

    def sh_halfplane_area(u, th):
        return P.intersection(sh_quad(u, th)).area

    # (1) per-pivot: certified bound must dominate a dense independent sweep
    tol_abs = area_P * Fraction(1, 1000)
    worst = 0.0
    for u in pivot_points(pts, edge_samples):
        b, _, _ = u_pivot(pts, u, tol_abs)          # no cutoff: we want the real bound
        bf = float(b)
        hi = max(sh_halfplane_area(u, 2 * math.pi * k / dense) for k in range(dense))
        if hi > bf * (1 + 1e-9):
            ok = False
            msgs.append(f"  PIVOT MISS at {(float(u[0]), float(u[1]))}: "
                        f"certified {bf:.6f} but shapely found {hi:.6f}")
        worst = max(worst, hi / bf if bf else 0.0)
    msgs.append(f"  (1) branch-and-bound vs dense shapely sweep, {dense} angles/pivot: "
                f"worst sample/bound = {worst:.6f} (must be <= 1)")

    # (1b) the two-pivot bound against a dense independent sweep
    pv = pivot_points(pts, edge_samples)
    if len(pv) >= 2:
        ranked = sorted(pv, key=lambda u: float(u_pivot(pts, u, tol_abs)[1]))[:3]
        m = 90
        worst2 = 0.0
        for i in range(len(ranked)):
            for j in range(i + 1, len(ranked)):
                u, v = ranked[i], ranked[j]
                b, _, _ = uv_pair(pts, u, v, tol_abs)
                bf = float(b)
                # dense 2-D sweep with shapely, from two independently built quads
                hi = 0.0
                for ka in range(m):
                    ra = P.intersection(sh_quad(u, 2 * math.pi * ka / m))
                    if ra.is_empty or ra.area <= hi:
                        continue          # cannot beat hi once intersected further
                    for kb in range(m):
                        av = ra.intersection(sh_quad(v, 2 * math.pi * kb / m)).area
                        if av > hi:
                            hi = av
                if hi > bf * (1 + 1e-9):
                    ok = False
                    msgs.append(f"  PAIR MISS: certified {bf:.6f} but shapely found {hi:.6f}")
                worst2 = max(worst2, hi / bf if bf else 0.0)
        msgs.append(f"  (1b) two-pivot 2-D branch and bound vs {m}x{m} shapely sweep: "
                    f"worst sample/bound = {worst2:.6f} (must be <= 1)")

    # (2) width bound vs a brute-force numerical integral of maxseg
    w = _candidate_directions(pts, 24)[0]
    exact = float(width_bound_dir(pts, w))
    p, q = float(w[0]), float(w[1])
    vx, vy = -q, p
    W2 = p * p + q * q
    fp = [(float(x), float(y)) for x, y in pts]
    T = [x * p + y * q for x, y in fp]
    lo, hi_t = min(T), max(T)
    steps = 4000
    total = 0.0
    for k in range(steps):
        t = lo + (hi_t - lo) * (k + 0.5) / steps
        cs = []
        for i in range(len(fp)):
            j = (i + 1) % len(fp)
            if (T[i] < t) == (T[j] < t):
                continue
            lam = (t - T[i]) / (T[j] - T[i])
            x = fp[i][0] + lam * (fp[j][0] - fp[i][0])
            y = fp[i][1] + lam * (fp[j][1] - fp[i][1])
            cs.append(x * vx + y * vy)
        cs.sort()
        if len(cs) >= 2 and len(cs) % 2 == 0:
            total += max(cs[2 * m + 1] - cs[2 * m] for m in range(len(cs) // 2))
    approx = total * (hi_t - lo) / steps / W2
    if exact < approx * (1 - 1e-6):
        ok = False
        msgs.append(f"  WIDTH MISMATCH: exact {exact:.6f} < numeric {approx:.6f}")
    msgs.append(f"  (2) width bound vs {steps}-step numeric integral: "
                f"exact {exact:.6f}, numeric {approx:.6f}, ratio {exact/approx:.6f} (must be >= 1)")
    return ok, msgs


def check_clip():
    from shapely.geometry import Polygon as ShPoly
    import random
    random.seed(7)
    shapes = [
        sawtooth(4),
        [(0, 0), (10, 0), (10, 10), (6, 10), (6, 7), (4, 7), (4, 10), (0, 10)],
        [(0, 0), (4, 0), (4, 4), (2, 1), (0, 4)],
        [(0, 0), (10, 0), (10, 6), (8, 6), (8, 1), (7, 1), (7, 6), (0, 6)],
    ]
    worst = 0.0
    for sh in shapes:
        P = ShPoly(sh).buffer(0)
        pts = [(Fraction(x), Fraction(y)) for x, y in sh]
        for _ in range(200):
            ux = random.uniform(min(p[0] for p in sh), max(p[0] for p in sh))
            uy = random.uniform(min(p[1] for p in sh), max(p[1] for p in sh))
            th = random.uniform(0, 2 * math.pi)
            d = (math.cos(th), math.sin(th))
            mine = float(halfplane_area(pts, (Fraction(ux), Fraction(uy)), d))
            # shapely reference: clip by a big rectangle on the left of the line
            big = 10 * max(abs(v) for p in sh for v in p) + 100
            nx, ny = -d[1], d[0]                  # left normal
            quad = [(ux + d[0] * big, uy + d[1] * big),
                    (ux - d[0] * big, uy - d[1] * big),
                    (ux - d[0] * big + nx * big, uy - d[1] * big + ny * big),
                    (ux + d[0] * big + nx * big, uy + d[1] * big + ny * big)]
            ref = P.intersection(ShPoly(quad).buffer(0)).area
            worst = max(worst, abs(mine - ref) / max(ref, 1e-9))

            # and the SEQUENTIAL two-half-plane clip used by the two-pivot bound
            th2 = random.uniform(0, 2 * math.pi)
            d2 = (math.cos(th2), math.sin(th2))
            ux2 = random.uniform(min(p[0] for p in sh), max(p[0] for p in sh))
            uy2 = random.uniform(min(p[1] for p in sh), max(p[1] for p in sh))
            mine2 = float(halfplane_area2(pts, (Fraction(ux), Fraction(uy)), d,
                                          (Fraction(ux2), Fraction(uy2)), d2))
            nx2, ny2 = -d2[1], d2[0]
            quad2 = [(ux2 + d2[0] * big, uy2 + d2[1] * big),
                     (ux2 - d2[0] * big, uy2 - d2[1] * big),
                     (ux2 - d2[0] * big + nx2 * big, uy2 - d2[1] * big + ny2 * big),
                     (ux2 + d2[0] * big + nx2 * big, uy2 + d2[1] * big + ny2 * big)]
            ref2 = P.intersection(ShPoly(quad).buffer(0)).intersection(
                ShPoly(quad2).buffer(0)).area
            worst = max(worst, abs(mine2 - ref2) / max(ref2, 1e-9))
    print(f"clipper self-test: worst relative disagreement with shapely = {worst:.3e}")
    return worst < 1e-9


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", nargs="*", type=int)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--suite", action="store_true", help="the 5-polygon C&Y regression set")
    ap.add_argument("--sawtooth", action="store_true", help="the sawtooth-bar ladder, A* = 45 exactly")
    ap.add_argument("--optima", nargs="?", const="",
                    help="compare against A* from a stored solve_chang_and_yap result "
                         "(default: newest); alone, runs every polygon with a known A*")
    ap.add_argument("--rel-tol", type=float, default=1e-3)
    ap.add_argument("--which", choices=("pivot", "width", "both"), default="both")
    ap.add_argument("--dirs", type=int, default=180, help="uniform slice directions to rank")
    ap.add_argument("--pairs", type=int, default=0,
                    help="two-pivot bound: number of top pivots to draw pairs from")
    ap.add_argument("--edge-samples", type=int, default=0,
                    help="extra pivots per edge (every boundary point is a legal pivot)")
    ap.add_argument("--check-clip", action="store_true")
    ap.add_argument("--verify", action="store_true",
                    help="independent audit of the bound (see verify())")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    if args.check_clip:
        return 0 if check_clip() else 1

    jobs = []          # (label, coords, known_A_star or None)

    if args.sawtooth:
        for t in (1, 2, 3, 4, 5, 6, 8, 12):
            jobs.append((f"sawtooth t={t}", sawtooth(t), 45.0))

    star = {}
    if args.optima is not None:
        from _common import newest_matching
        path = args.optima or newest_matching("solve_chang_and_yap_")
        if not path:
            print("no solve_chang_and_yap result; run solve_chang_and_yap.py first")
            return 2
        conn = sqlite3.connect(path)
        star = {pid: float(Fraction(ae)) for pid, ae in conn.execute(
            "SELECT pid, area_exact FROM optima WHERE status='ok' AND area_exact IS NOT NULL")}
        conn.close()

    if args.suite:
        names = ["Ireland", "Vietnam", "Yemen", "New Zealand_1", "New Zealand_2"]
        known = dict(zip(names, [6.138621, 8.810944, 29.292857, 12.802720, 7.286795]))
        from _common import require_corpus
        conn = sqlite3.connect(require_corpus())
        for nm in names:
            row = conn.execute("SELECT id,name,geometry_wkb FROM polygons WHERE name=?",
                               (nm,)).fetchone()
            g = shapely_wkb.loads(row[2])
            if g.geom_type == "MultiPolygon":
                g = max(g.geoms, key=lambda x: x.area)
            jobs.append((f"[{row[0]}] {nm}", list(g.exterior.coords)[:-1], known[nm]))
        conn.close()

    if args.ids or (star and not args.suite and not args.sawtooth):
        for pid, name, coords in load_corpus(args.ids or sorted(star), args.limit):
            jobs.append((f"[{pid}] {name}", coords, star.get(pid)))

    if not jobs:
        print("nothing to do; pass --suite, --sawtooth, --ids or --optima")
        return 2

    if args.verify:
        all_ok = True
        for label, coords, _ in jobs:
            clean = pre_checks(coords)
            print(f"{label}  (n={len(clean)})", flush=True)
            ok, msgs = verify(clean, edge_samples=args.edge_samples)
            for m in msgs:
                print(m, flush=True)
            all_ok &= ok
        print(f"\nVERIFY: {'PASS' if all_ok else '*** FAIL ***'}")
        return 0 if all_ok else 1

    want_p = args.which in ("pivot", "both")
    want_w = args.which in ("width", "both")
    print(f"{'polygon':<24} {'n':>4} {'refl':>5} {'area(P)':>12} {'PIVOT':>12} {'WIDTH':>12} "
          f"{'BEST':>12} {'known A*':>12} {'p/A*':>7} {'w/A*':>7} {'best/A*':>8} {'ok':>4} {'s':>6}")
    rows, bad = [], 0
    for label, coords, known in jobs:
        clean = pre_checks(coords)
        pts = [(Fraction(x), Fraction(y)) for x, y in clean]
        areaP = abs(shoelace(pts))
        t0 = time.time()
        bp = (upper_bound(clean, rel_tol=args.rel_tol, verbose=args.verbose,
                          edge_samples=args.edge_samples,
                          pairs=args.pairs)[0] if want_p else areaP)
        bw = width_bound(pts, n_uniform=args.dirs) if want_w else areaP
        if bw is None:
            bw = areaP
        dt = time.time() - t0
        best = min(bp, bw)
        bpf, bwf, bf = float(bp), float(bw), float(best)
        if known:
            ok = bf >= known - 1e-9 * max(1.0, known)
            bad += (not ok)
            rows.append((bpf / known, bwf / known, bf / known))
            print(f"{label:<24} {len(clean):>4} {len(reflex_indices(pts)):>5} {float(areaP):>12.4f} "
                  f"{bpf:>12.4f} {bwf:>12.4f} {bf:>12.4f} {known:>12.6f} "
                  f"{bpf/known:>7.3f} {bwf/known:>7.3f} {bf/known:>8.4f} "
                  f"{'OK' if ok else 'FAIL':>4} {dt:>6.1f}")
        else:
            print(f"{label:<24} {len(clean):>4} {len(reflex_indices(pts)):>5} {float(areaP):>12.4f} "
                  f"{bpf:>12.4f} {bwf:>12.4f} {bf:>12.4f} {'-':>12} "
                  f"{'-':>7} {'-':>7} {'-':>8} {'-':>4} {dt:>6.1f}")

    if rows:
        def stats(k):
            v = sorted(r[k] for r in rows)
            return v[0], v[len(v) // 2], v[-1]
        for k, nm in ((0, "pivot"), (1, "width"), (2, "BEST ")):
            lo, md, hi = stats(k)
            print(f"\n  {nm}: min {lo:.4f}  median {md:.4f}  max {hi:.4f}", end="")
        print(f"\n\n{len(rows)} with a known A*.")
        print(f"VALIDITY: {bad} bound(s) below a known A*  "
              f"({'SOUND' if bad == 0 else '*** UNSOUND -- the bound is wrong ***'})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

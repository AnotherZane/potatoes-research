#include "potato_peeling_core/algorithms/chang_and_yap/chang_and_yap.hpp"
#include <CGAL/Boolean_set_operations_2.h>
#include "potato_peeling_core/geometry.hpp"
#include <iterator>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <variant>

#include "potato_peeling_core/algorithms/chang_and_yap/memprobe.hpp"

namespace potato_peeling::chang_and_yap {

    // pre_checks() orients P counter-clockwise, and everything here assumes it:
    // reflex corners are right turns, C+ is left of its directed chord, and
    // signed areas are CCW-positive as in the paper (p174). The paper's clockwise
    // corner order R(C,C') (p173) is the increasing-index arc here.

#pragma region Helpers
    // Shared core of the containment test. Assumes no self intersections.
    //
    // `endpoints_known_inside` lets callers that already know both endpoints lie
    // in closed P skip two O(n) bounded_side scans; generate_extremal_chords
    // passes polygon vertices, which are on the boundary.

    static bool segment_in_closed_polygon(const Polygon &polygon, const Segment &segment,
                                          const bool endpoints_known_inside) {
        const Point &src = segment.source();
        const Point &tgt = segment.target();

        // If any of the endpoints are outside, it's not contained
        if (!endpoints_known_inside
            && (polygon.bounded_side(src) == CGAL::ON_UNBOUNDED_SIDE
                || polygon.bounded_side(tgt) == CGAL::ON_UNBOUNDED_SIDE)) {
            return false;
        }

        // Each crossing carries its sort key, computed once per point, and its
        // provenance, which the sweep below uses to tell a transversal crossing
        // from a degenerate one.
        struct Cross {
            FT key;  // squared distance from src — the sort key
            Point p;
            int producers;    // how many polygon edges produced this point
            bool degenerate;  // vertex hit, collinear overlap, or a segment endpoint
            Point ea, eb;     // the producing edge, when producers == 1
        };
        std::vector<Cross> intersections;
        intersections.reserve(polygon.size() + 2);
        const auto add = [&](const Point &p, const int producers, const bool degen, const Point &ea,
                             const Point &eb) {
            intersections.push_back({CGAL::squared_distance(src, p), p, producers, degen, ea, eb});
        };
        // The segment's own endpoints are never toggle sites: nothing is known about
        // how dP behaves at them, so they are marked degenerate and the cell they
        // open is anchored with an exact test.
        add(src, 0, true, src, src);
        add(tgt, 0, true, tgt, tgt);

        for (auto edge = polygon.edges_begin(); edge != polygon.edges_end(); ++edge) {
            auto result = CGAL::intersection(*edge, segment);
            if (!result) continue;

            if (const Point *p = std::get_if<Point>(&*result)) {
                // A hit AT an edge endpoint is a CORNER of P. The boundary may touch
                // the segment's line there and turn back without crossing it, so a
                // corner is never a toggle site.
                const bool at_corner = (*p == edge->source() || *p == edge->target());
                add(*p, 1, at_corner, edge->source(), edge->target());
            } else if (const Segment *s = std::get_if<Segment>(&*result)) {
                // The segment runs ALONG this edge. Both ends of the overlap bound a
                // stretch of dP rather than crossing it.
                add(s->source(), 1, true, edge->source(), edge->target());
                add(s->target(), 1, true, edge->source(), edge->target());
            }
        }

        // Sort the intersection points by distance from src. They all lie on
        // the (collinear) segment line, so distance order equals order along the
        // segment — the sequencing the sub-segment containment test below needs.
        std::sort(intersections.begin(), intersections.end(),
                  [](const Cross &a, const Cross &b) { return a.key < b.key; });

        // Merge duplicates (adjacent after the sort), accumulating provenance: a
        // point two edges produced is a corner of P, and a point degenerate for
        // any producer stays degenerate.
        std::vector<Cross> merged;
        merged.reserve(intersections.size());
        for (const auto &c : intersections) {
            if (!merged.empty() && merged.back().p == c.p) {
                merged.back().producers += c.producers;
                merged.back().degenerate = merged.back().degenerate || c.degenerate;
            } else {
                merged.push_back(c);
            }
        }

        // Sweep the sub-segments. At a transversal crossing (one edge, interior to
        // it, segment ends strictly on opposite sides of its line) the side
        // toggles in O(1); anything else re-anchors with an exact O(n)
        // bounded_side(). O(n log n) in general position.
        CGAL::Bounded_side side = CGAL::ON_BOUNDARY;
        bool side_known = false;  // false => the next cell must be anchored exactly

        for (std::size_t i = 0; i + 1 < merged.size(); ++i) {
            if (side_known) {
                const Cross &x = merged[i];  // the crossing separating cell i-1 from cell i
                bool transversal = false;
                if (!x.degenerate && x.producers == 1) {
                    const CGAL::Orientation o_src = CGAL::orientation(x.ea, x.eb, src);
                    const CGAL::Orientation o_tgt = CGAL::orientation(x.ea, x.eb, tgt);
                    transversal
                        = o_src != CGAL::COLLINEAR && o_tgt != CGAL::COLLINEAR && o_src != o_tgt;
                }
                if (transversal)
                    side = (side == CGAL::ON_BOUNDED_SIDE) ? CGAL::ON_UNBOUNDED_SIDE
                                                           : CGAL::ON_BOUNDED_SIDE;
                else
                    side_known = false;
            }

            if (!side_known) {
                // Since we collect all intersections, it is enough to check just the midpoint of a segment
                side = polygon.bounded_side(CGAL::midpoint(merged[i].p, merged[i + 1].p));
                side_known = side != CGAL::ON_BOUNDARY;
            }

            // If even one sub-segment midpoint is outside, the whole segment is not contained
            if (side == CGAL::ON_UNBOUNDED_SIDE) {
                return false;
            }
        }

        return true;
    }

    // Assumes no self intersections
    bool is_contained(const Polygon &polygon, const Segment &segment) {
        return segment_in_closed_polygon(polygon, segment, false);
    }

    // Membership in the chord (p158/p173), which §2 defines as the maximal
    // segment of the line through v_i, v_j contained in P. The vertex segment
    // [v_i, v_j] is narrower and would reject bounding nodes the optimum needs.
    // Chord::start/end hold the maximal endpoints.
    bool on_chord_segment(const Point &x, const Chord &c) {
        return CGAL::are_ordered_along_line(c.start, x, c.end);
    }

    // For any chord C_{i,j} ∈ Ξ, we call the point midway between
    // the double-pivot v_i, and v_j the reference point of C_{i,j}
    Point reference_point(const Chord &C, const Polygon &poly) {
        return CGAL::midpoint(poly[C.i], poly[C.j]);
    }

    // Signed area of triangle (a, b, c): positive when CCW. This is Section 6's
    // signed area 2*AREA = sum x_i (y_{i+1} - y_i) specialized to a triangle.
    // Callers rely on the SIGN (the convex-corner guard in compute_M_stage_1, the
    // signed sum in length_1_chain_area), so this must stay signed.
    FT triangle_area_signed(const Point &a, const Point &b, const Point &c) {
        return CGAL::area(a, b, c);
    }

#pragma endregion

#pragma region Step 1

    // Does the chord extend past corner v_k in direction d (away from the other
    // corner)? True iff d lies in the closed interior cone at v_k. A chord through
    // a reflex corner does not always extend past it: the line can graze the
    // corner and leave P there.
    bool extends_past_corner(const Polygon &polygon, const std::size_t k, const Vector &d) {
        const std::size_t n = polygon.size();
        const Vector e_out = polygon[(k + 1) % n] - polygon[k];
        const Vector e_in = polygon[(k + n - 1) % n] - polygon[k];
        const auto cross2
            = [](const Vector &u, const Vector &v) -> FT { return u.x() * v.y() - u.y() * v.x(); };
        // Interior directions sweep CCW from e_out to e_in. Convex corner
        // (cone < pi): d between both bounds; reflex corner (cone > pi):
        // d outside the exterior cone, i.e. past either bound. Boundary
        // (collinear with an edge) counts as inside — closed P.
        if (cross2(e_out, e_in) > FT(0))
            return cross2(e_out, d) >= FT(0) && cross2(d, e_in) >= FT(0);
        return cross2(e_out, d) >= FT(0) || cross2(d, e_in) >= FT(0);
    }

    // Endpoints of the maximal chord on the line through v_i, v_j: it extends past
    // any corner it merely passes through and stops where it leaves P. It does not
    // decide which chords exist (generate_extremal_chords does); on_chord_segment
    // tests membership in it. Returns nullopt when the maximal piece cannot be
    // recovered cleanly; the caller then uses the vertex endpoints.
    std::optional<std::pair<Point, Point>> maximal_chord_endpoints(const Polygon &polygon,
                                                                   std::size_t i, std::size_t j) {
        const Line geometric_line(polygon[i], polygon[j]);

        const auto &bbox = polygon.bbox();
        Point bottom_left(bbox.xmin(), bbox.ymin()), bottom_right(bbox.xmax(), bbox.ymin());
        Point top_right(bbox.xmax(), bbox.ymax()), top_left(bbox.xmin(), bbox.ymax());

        std::vector<Segment> box_edges
            = {Segment(bottom_left, bottom_right), Segment(bottom_right, top_right),
               Segment(top_right, top_left), Segment(top_left, bottom_left)};

        std::vector<Point> box_intersections;
        for (const auto &edge : box_edges) {
            if (auto res = CGAL::intersection(geometric_line, edge)) {
                if (const Point *p = std::get_if<Point>(&*res)) {
                    box_intersections.push_back(*p);
                } else if (const Segment *s = std::get_if<Segment>(&*res)) {
                    box_intersections.push_back(s->source());
                    box_intersections.push_back(s->target());
                }
            }
        }

        // The line must cross the box at least twice to intersect the polygon interior.
        if (box_intersections.size() < 2) return std::nullopt;

        std::sort(box_intersections.begin(), box_intersections.end(),
                  SortAlongLine(geometric_line.direction()));
        Segment clipped_segment(box_intersections.front(), box_intersections.back());

        // Crossings carry provenance, as in segment_in_closed_polygon.
        struct LineCross {
            Point p;
            int producers;
            bool degenerate;
            Point ea, eb;
        };
        std::vector<LineCross> crossings;
        crossings.reserve(polygon.size() + 2);
        for (auto edge = polygon.edges_begin(); edge != polygon.edges_end(); ++edge) {
            auto result = CGAL::intersection(*edge, clipped_segment);
            if (!result) continue;

            if (const Point *p = std::get_if<Point>(&*result)) {
                const bool at_corner = (*p == edge->source() || *p == edge->target());
                crossings.push_back({*p, 1, at_corner, edge->source(), edge->target()});
            } else if (const Segment *s = std::get_if<Segment>(&*result)) {
                crossings.push_back({s->source(), 1, true, edge->source(), edge->target()});
                crossings.push_back({s->target(), 1, true, edge->source(), edge->target()});
            }
        }

        // Sort ∂P-crossings along the line so consecutive pairs bound sub-segments.
        const SortAlongLine along(geometric_line.direction());
        std::sort(crossings.begin(), crossings.end(),
                  [&](const LineCross &a, const LineCross &b) { return along(a.p, b.p); });

        // Merge duplicates, accumulating provenance.
        std::vector<LineCross> intersections;
        intersections.reserve(crossings.size());
        for (const auto &x : crossings) {
            if (!intersections.empty() && intersections.back().p == x.p) {
                intersections.back().producers += x.producers;
                intersections.back().degenerate = intersections.back().degenerate || x.degenerate;
            } else {
                intersections.push_back(x);
            }
        }

        if (intersections.size() < 2) return std::nullopt;

        // p158: "A chord of P is a maximal line segment fully contained in P."
        // Each piece between consecutive crossings is wholly inside or outside
        // closed P (on the boundary counts as inside: a chord can run along an
        // edge). Merging the inside pieces into runs makes the chord maximal. Sides
        // are found by the same sweep as in segment_in_closed_polygon.
        const Point &csrc = clipped_segment.source();
        const Point &ctgt = clipped_segment.target();
        CGAL::Bounded_side side = CGAL::ON_BOUNDARY;
        bool side_known = false;

        std::vector<Segment> segments;
        for (std::size_t k = 0; k + 1 < intersections.size(); ++k) {
            if (side_known) {
                const LineCross &x = intersections[k];
                bool transversal = false;
                if (!x.degenerate && x.producers == 1) {
                    const CGAL::Orientation o_src = CGAL::orientation(x.ea, x.eb, csrc);
                    const CGAL::Orientation o_tgt = CGAL::orientation(x.ea, x.eb, ctgt);
                    transversal
                        = o_src != CGAL::COLLINEAR && o_tgt != CGAL::COLLINEAR && o_src != o_tgt;
                }
                if (transversal)
                    side = (side == CGAL::ON_BOUNDED_SIDE) ? CGAL::ON_UNBOUNDED_SIDE
                                                           : CGAL::ON_BOUNDED_SIDE;
                else
                    side_known = false;
            }
            if (!side_known) {
                side = polygon.bounded_side(
                    CGAL::midpoint(intersections[k].p, intersections[k + 1].p));
                side_known = side != CGAL::ON_BOUNDARY;
            }

            if (side == CGAL::ON_UNBOUNDED_SIDE) continue;

            if (!segments.empty() && segments.back().end() == intersections[k].p)
                segments.back() = Segment(segments.back().start(), intersections[k + 1].p);
            else
                segments.emplace_back(intersections[k].p, intersections[k + 1].p);
        }

        std::vector<Segment> chordz;
        for (const auto &segment : segments) {
            if (CGAL::are_ordered_along_line(segment.start(), polygon[i], polygon[j])
                && CGAL::are_ordered_along_line(polygon[i], polygon[j], segment.end())) {
                chordz.push_back(segment);
            }
        }

        if (chordz.size() == 1) {
            return std::make_pair(chordz.front().start(), chordz.front().end());
        }
        // 0 or >1 spanning pieces: cannot cleanly recover the maximal chord.
        return std::nullopt;
    }

    std::vector<Chord> generate_extremal_chords(const Polygon &polygon,
                                                const std::vector<bool> &reflex_flags) {
        const std::size_t n = polygon.size();

        // Chords are oriented (p173), but building one does not depend on its
        // orientation, so each unordered pair is solved once (oriented i -> j for
        // i < j) and emitted both ways.
        std::vector<std::optional<std::pair<Point, Point>>> solved(n * n);

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                // O(1) necessary condition, before the O(n) test below:
                // [v_i,v_j] cannot leave P at its own endpoints, so v_j - v_i
                // must lie in the closed interior cone at v_i, and v_i - v_j in
                // the cone at v_j.
                if (!extends_past_corner(polygon, i, polygon[j] - polygon[i])
                    || !extends_past_corner(polygon, j, polygon[i] - polygon[j]))
                    continue;

                // Existence test: a chord C_{i,j} exists iff the VERTEX segment
                // [v_i, v_j] is contained in P (p173's definedness condition —
                // maximality never changes which chords exist, only their extent).
                // The endpoints are corners of P, so their test is skipped.
                if (const Segment vseg(polygon[i], polygon[j]);
                    !segment_in_closed_polygon(polygon, vseg, true)) {
                    continue;
                }

                // Maximal endpoints (vertex endpoints on failure). A chord can only
                // run past an endpoint whose interior angle is at least pi, so two
                // convex endpoints already bound the maximal chord.
                Point start = polygon[i], end = polygon[j];
                if (reflex_flags[i] || reflex_flags[j]) {
                    if (const auto ext = maximal_chord_endpoints(polygon, i, j)) {
                        start = ext->first;
                        end = ext->second;
                    }
                }

                solved[i * n + j] = std::make_pair(start, end);
            }
        }

        // Emit in (i ascending, then j ascending) order. Xi's order matters:
        // composition commits to the first valid decomposition it finds.
        std::vector<Chord> chords;
        chords.reserve(n * (n - 1));
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (i == j) continue;

                if (i < j) {
                    if (const auto &e = solved[i * n + j])
                        chords.emplace_back(e->first, e->second, i, j);
                } else {
                    if (const auto &e = solved[j * n + i])
                        chords.emplace_back(e->second, e->first, i, j);
                }
            }
        }
        return chords;
    }

#pragma endregion

#pragma region Step 4

    // Admissibility (Section 6): the pair (C_{i,j}, C_{k,l}) is admissible iff
    // the corners v_i, v_j, v_k, v_l occur in this cyclic order on the boundary
    // of P: with positions cyclic from v_i, pos(j) <= pos(k) < pos(l), where
    // l == i (C' ends where C starts, Step III) counts as maximal. j == k (C ends
    // where C' starts) is allowed.
    bool is_admissible_pair(const Chord &c, const Chord &c_prime, const std::size_t &n) {
        if (c == c_prime) return false;

        const std::size_t i = c.i, j = c.j;
        const std::size_t k = c_prime.i, l = c_prime.j;

        // j == k && l == i is C_{i,j} paired with its reverse C_{j,i}: only two
        // distinct corners, and Step III needs m > 2 boundary chords.
        if (j == k && l == i) return false;

        // positions measured cyclically from v_i; pos(v_i) = 0
        auto cyclic_less = [n](const std::size_t a, const std::size_t b, const std::size_t start) {
            return ((a - start + n) % n) < ((b - start + n) % n);
        };

        // pos(j) <= pos(k); equality is the j == k adjacency boundary case
        const bool j_leq_k = (j == k) || cyclic_less(j, k, i);
        // pos(k) < pos(l); l == i is the full-wrap closing case, handled
        // explicitly because cyclic_less(k, i, i) is always false. k != i holds
        // for any real chord reaching here, so the closing case is admitted.
        const bool k_less_l = (l == i) ? (k != i) : cyclic_less(k, l, i);

        return j_leq_k && k_less_l;
    }

    // R(C_{i,j}, C_{k,l}) — the corners of P between the two chords, on the
    // side where the chain lives (Section 6: "the corners of P clockwise from
    // v_j to v_k inclusive").
    //
    // The paper's clockwise arc is the increasing-index arc on our CCW polygon.
    std::vector<std::size_t> compute_R(const Chord &c, const Chord &c_prime, const std::size_t &n) {
        const auto &j = c.j;
        const auto &k = c_prime.i;

        std::vector<std::size_t> R;
        for (std::size_t idx = j; idx != k; idx = (idx + 1) % n) R.push_back(idx);
        R.push_back(k);
        return R;
    }

#pragma endregion

#pragma region Step 5

    // Chord c directed as (A,B) with C+ on the left: the half-plane whose small
    // disc at the reflex corner lies inside P (§2 p158). Oriented by whichever
    // endpoint is reflex.
    std::pair<Point, Point> directed_C_plus(const Chord &c, const Polygon &P,
                                            const std::vector<bool> &reflex_flags) {
        Point A = P[c.i];
        Point B = P[c.j];

        std::size_t pivot_idx;
        if (reflex_flags[c.i])
            pivot_idx = c.i;
        else if (reflex_flags[c.j])
            pivot_idx = c.j;
        else {
            // Edge chord: the interior is left of the CCW boundary edge.
            // A convex–convex diagonal has no C+ in the paper (Section 2) and
            // gets this arbitrary orientation. Such chords occur in Ξ but are
            // filtered out of the optimum elsewhere.
            return {A, B};
        }

        const std::size_t prev = (pivot_idx + P.size() - 1) % P.size();
        const std::size_t next = (pivot_idx + 1) % P.size();

        // Interior direction at the reflex corner.
        const Vector v_in(-(P[prev].x() + P[next].x() - 2 * P[pivot_idx].x()),
                          -(P[prev].y() + P[next].y() - 2 * P[pivot_idx].y()));
        const Point test_pt(P[pivot_idx].x() + v_in.x(), P[pivot_idx].y() + v_in.y());

        // Direct (A,B) so that test_pt (interior side) is LEFT.
        if (CGAL::orientation(A, B, test_pt) == CGAL::RIGHT_TURN) std::swap(A, B);
        return {A, B};
    }

    // r strictly on the interior side: the open half-plane. Emptiness is defined
    // in §4 (p165) via the open region ρ(𝒞), so a point on a bounding line is not
    // an obstruction. This is why R may include the on-line endpoints v_j, v_k.
    // (The reconstruction uses §2's closed C+.)
    bool in_C_plus(const Point &r, const Chord &c, const Polygon &P,
                   const std::vector<bool> &reflex_flags) {
        auto [A, B] = directed_C_plus(c, P, reflex_flags);
        return Line(A, B).oriented_side(r) == CGAL::ON_POSITIVE_SIDE;
    }

    // The length-0 test (p174: "R(C,C') ∩ C+ ∩ C'+ is empty"), by §4's definition
    // (p165): no corner of R strictly inside both half-planes. It must be open,
    // since R includes v_j on L0 and v_k on L1.
    bool R_intersect_halfplanes_empty(const std::vector<std::size_t> &R, const Chord &c,
                                      const Chord &c_prime, const Polygon &P,
                                      const std::vector<bool> &reflex_flags) {
        auto [A, B] = directed_C_plus(c, P, reflex_flags);
        const auto c_line = Line(A, B);
        auto [C, D] = directed_C_plus(c_prime, P, reflex_flags);
        const auto c_prime_line = Line(C, D);

        for (const std::size_t idx : R) {
            const Point &r = P[idx];
            if (c_line.oriented_side(r) == CGAL::ON_POSITIVE_SIDE
                && c_prime_line.oriented_side(r) == CGAL::ON_POSITIVE_SIDE)
                return false;
        }
        return true;
    }

    // Each node x_i is affine in x0 by the lever rule x_i = 2 p_i - x_{i-1}
    // (Section 3); iterating it gives the per-node offsets A_i. The paper leaves
    // the closed form implicit.
    std::vector<Vector> compute_lever_offsets(const std::vector<Point> &pivots,
                                              const Point &origin) {
        // A_i as vector from origin
        // A_0 = 0 (x0 = origin + t*dir → offset = 0)
        std::vector<Vector> A;
        A.emplace_back(0, 0);  // A_0

        Vector prev(0, 0);
        for (const auto &pivot : pivots) {
            // Ai = 2*pi - A_{i-1} follows from
            // A_i = 2*p_i - A_{i-1} (accumulated point offset from origin)
            Vector Ai(2 * (pivot.x() - origin.x()) - prev.x(),
                      2 * (pivot.y() - origin.y()) - prev.y());
            A.push_back(Ai);
            prev = Ai;
        }
        return A;
    }

    // The first node x0 of the unique balanced (L0,L1)-chain (Theorem 5). Every
    // node is affine in x0 = origin + t*dir:
    //     x_m = origin + A[m] + (-1)^m * t * dir,
    // and x_m ∈ L1 fixes t. nullopt iff L0 ∥ L1.
    std::optional<Point> solve_balanced_chain_x0(const Line &L0, const Line &L1,
                                                 const std::vector<Point> &pivots) {
        const Point origin = L0.point(0);
        const Vector dir = L0.to_vector();

        const auto A = compute_lever_offsets(pivots, origin);
        const std::size_t m = pivots.size();
        const int sign_m = (m % 2 == 0) ? 1 : -1;

        // x_m evaluated at t = 0 (x0 = origin):
        const Point xm_base(origin.x() + A[m].x(), origin.y() + A[m].y());

        // L1 in implicit form a*x + b*y + c = 0.
        const FT f0 = L1.a() * xm_base.x() + L1.b() * xm_base.y() + L1.c();
        const FT df = FT(sign_m) * (L1.a() * dir.x() + L1.b() * dir.y());

        if (df == FT(0)) return std::nullopt;  // L0 ∥ L1

        const FT t = -f0 / df;
        return Point(origin.x() + t * dir.x(), origin.y() + t * dir.y());
    }

    // type (a) critical point — [x_{i-1}, xi]
    // passes through corner r — collinearity condition
    // (r - x_{i-1}) × (xi - x_{i-1}) = 0 — Section 3
    std::optional<FT>
    compute_critical_param_type_a(const Point &r,
                                  const std::size_t i,  // strip/node index (1-based)
                                  const std::vector<Vector> &A, const Point &origin,
                                  const Vector &dir,
                                  const int sign_i)  // (-1)^i
    {
        auto cross_2d
            = [](const Vector &u, const Vector &v) -> FT { return u.x() * v.y() - u.y() * v.x(); };

        // P = r - x_{i-1} constant part,
        // Q = x_i - x_{i-1} constant part — Section 3/5
        const Vector P(r.x() - origin.x() - A[i - 1].x(), r.y() - origin.y() - A[i - 1].y());
        const Vector Q(A[i].x() - A[i - 1].x(), A[i].y() - A[i - 1].y());

        // δ = (-1)^{i-1} = -(-1)^i = -sign_i
        const int delta = -sign_i;

        const FT numerator = cross_2d(P, Q);
        const FT denominator = FT(delta) * (FT(2) * cross_2d(P, dir) + cross_2d(dir, Q));

        // Degenerate: no critical event from this corner and strip. (FT cannot
        // hold a non-finite value.)
        if (denominator == FT(0)) return std::nullopt;

        return numerator / denominator;
    }

    // type (b) critical point — x_{k-1}, xk, x_{k+1}
    // are collinear — Section 3
    // dk × dk1 = 0 where dk = xk - x_{k-1}, dk1 = x_{k+1} - xk
    // origin cancels in differences x_k - x_{k-1},
    // so Bk = A[k] - A[k-1] with no origin term — Section 3/5
    std::optional<FT> compute_critical_param_type_b(const std::size_t k,
                                                    const std::vector<Vector> &A,
                                                    const Point &origin, const Vector &dir,
                                                    const int sign_k)  // (-1)^k
    {
        auto cross_2d
            = [](const Vector &u, const Vector &v) -> FT { return u.x() * v.y() - u.y() * v.x(); };

        // Bk = A[k] - A[k-1] — constant part of
        // dk = xk - x_{k-1}; origin cancels in difference
        const Vector Bk(A[k].x() - A[k - 1].x(), A[k].y() - A[k - 1].y());

        // Bk1 = A[k+1] - A[k] — constant part of
        // dk1 = x_{k+1} - xk; origin cancels in difference
        const Vector Bk1(A[k + 1].x() - A[k].x(), A[k + 1].y() - A[k].y());

        // dk × dk1 = 0 expands to
        // Bk×Bk1 - 2σt·((Bk+Bk1)×dir) = 0
        // where σ = (-1)^k — Section 3
        const FT numerator = cross_2d(Bk, Bk1);
        const FT denom_base = cross_2d(Vector(Bk.x() + Bk1.x(), Bk.y() + Bk1.y()), dir);
        const FT denominator = FT(2) * FT(sign_k) * denom_base;

        // Degenerate.
        if (denominator == FT(0)) return std::nullopt;

        return numerator / denominator;
    }

    // Sorted, de-duplicated critical parameters t of the (L0, p1..pm)-paths,
    // measured along L0 as offsets from L0.point(0) in units of L0.to_vector().
    // A critical parameter is where the path changes combinatorial structure:
    // type (a) a node segment [x_{i-1}, x_i] passes through a corner r in R, or
    // type (b) three consecutive nodes become collinear (Section 3, Lemma 7).
    // The critical intervals are the maximal sub-intervals between consecutive
    // critical parameters (plus the two unbounded end rays).
    std::vector<FT> compute_critical_params(const Line &L0, const std::vector<Point> &pivots,
                                            const std::vector<std::size_t> &R,
                                            const Polygon &poly) {
        const Point origin = L0.point(0);
        const Vector dir = L0.to_vector();

        // "pick arbitrary x0 on L0, form unique path" — Lemma 7
        const auto A = compute_lever_offsets(pivots, origin);

        // (-1)^i per node index — sign alternation from the lever rule
        std::vector<int> signs(pivots.size() + 1);
        for (std::size_t i = 0; i <= pivots.size(); i++) signs[i] = (i % 2 == 0) ? 1 : -1;

        std::vector<FT> critical_params;

        // Strips determined by the parallel lines L0(x_i) — Lemma 7
        const Vector normal(-dir.y(), dir.x());
        std::vector<Point> nodes = {origin};
        for (const auto &pivot : pivots) {
            Point xj(2 * pivot.x() - nodes.back().x(), 2 * pivot.y() - nodes.back().y());
            nodes.push_back(xj);
        }

        // type (a): node segment through a corner r in R - {pivots} — Section 3.
        // Node heights above L are not monotone in general, so strips OVERLAP
        // and a corner can be swept by several segments, so every containing
        // strip is checked. The band test is only a necessary condition, so each candidate
        // is verified exactly: at its solved t the corner must lie ON the
        // segment [x_i, x_{i+1}], not merely on its line (this also filters
        // the degenerate equal-height strip, where collinearity alone does not
        // imply betweenness).
        for (const std::size_t idx : R) {
            const Point &r = poly[idx];
            bool is_pivot = false;
            for (const auto &p : pivots)
                if (r == p) {
                    is_pivot = true;
                    break;
                }
            if (is_pivot) continue;

            FT r_proj = normal * (r - origin);
            for (std::size_t i = 0; i + 1 < nodes.size(); i++) {
                FT lo = normal * (nodes[i] - origin);
                FT hi = normal * (nodes[i + 1] - origin);
                if (lo > hi) std::swap(lo, hi);
                if (r_proj < lo || r_proj > hi) continue;

                // nullopt => degenerate => skip
                const auto t
                    = compute_critical_param_type_a(r, i + 1, A, origin, dir, signs[i + 1]);
                if (!t) continue;

                const Point xi(origin.x() + A[i].x() + FT(signs[i]) * *t * dir.x(),
                               origin.y() + A[i].y() + FT(signs[i]) * *t * dir.y());
                const Point xi1(origin.x() + A[i + 1].x() + FT(signs[i + 1]) * *t * dir.x(),
                                origin.y() + A[i + 1].y() + FT(signs[i + 1]) * *t * dir.y());
                if (CGAL::are_ordered_along_line(xi, r, xi1)) critical_params.push_back(*t);
            }
        }

        // type (b): three consecutive nodes collinear — Section 3
        for (std::size_t k = 1; k < pivots.size(); k++)
            if (auto t = compute_critical_param_type_b(k, A, origin, dir, signs[k]))
                critical_params.push_back(*t);

        std::sort(critical_params.begin(), critical_params.end());
        critical_params.erase(std::unique(critical_params.begin(), critical_params.end()),
                              critical_params.end());
        return critical_params;
    }

    // (L, p1..pm)-equivalence (Section 3 / Section 5). Two points are
    // equivalent iff the paths they determine share a critical interval, i.e.
    // no critical parameter lies strictly between their x0-parameters.
    //
    // x and y lie on the node line L(x_m), through the last node parallel to L.
    // Critical parameters are positions of x0, so a point pt on L(x_m) maps to
    // t = (-1)^m (proj(pt) - proj(A[m])), proj being projection onto dir.
    // `params` is the sorted output of compute_critical_params, memoised by the
    // caller (Lemma 8, p173).
    Vector last_lever_offset(const std::vector<Point> &pivots, const Point &origin) {
        // A_0 = 0; A_i = 2*(p_i - origin) - A_{i-1}, folded down to A_m.
        Vector acc(0, 0);
        for (const auto &pivot : pivots)
            acc = Vector(2 * (pivot.x() - origin.x()) - acc.x(),
                         2 * (pivot.y() - origin.y()) - acc.y());
        return acc;
    }

    bool equivalent_given_params(const Line &L, const std::size_t m, const Vector &last_offset,
                                 const std::vector<FT> &params, const Point &x, const Point &y) {
        const Point o = L.point(0);
        const Vector dir = L.to_vector();
        const FT dsq = dir * dir;
        if (dsq == FT(0)) return x == y;  // degenerate line (defensive)

        const FT a_proj = (last_offset * dir) / dsq;
        const FT sign = (m % 2 == 0) ? FT(1) : FT(-1);

        const FT tx = sign * (((x - o) * dir) / dsq - a_proj);
        const FT ty = sign * (((y - o) * dir) / dsq - a_proj);
        FT lo = tx, hi = ty;
        if (lo > hi) std::swap(lo, hi);

        // Not equivalent iff some param lies in (lo, hi): the smallest param > lo
        // separates x and y iff it is also < hi.
        const auto it = std::upper_bound(params.begin(), params.end(), lo);
        return !(it != params.end() && *it < hi);
    }

    std::optional<ChainEntry> try_length_0_chain(const std::vector<std::size_t> &R,
                                                 const std::vector<Chord> &Xi, const std::size_t a,
                                                 const std::size_t b, const Polygon &poly,
                                                 const std::vector<bool> &reflex_flags) {
        const Chord &c = Xi[a];
        const Chord &c_prime = Xi[b];

        const Line L0(poly[c.i], poly[c.j]);
        const Line L1(poly[c_prime.i], poly[c_prime.j]);

        // The chain has length 0 precisely where R(C, C') ∩ C+ ∩ C'+ is empty
        if (!R_intersect_halfplanes_empty(R, c, c_prime, poly, reflex_flags)) return std::nullopt;

        ChainEntry entry;

        // Length 0: a single node x0 = xm = L0 ∩ L1. Parallel or coincident lines
        // give no point, hence no length-0 chain.
        const auto inter = CGAL::intersection(L0, L1);
        if (!inter) return std::nullopt;

        const Point *pt = std::get_if<Point>(&*inter);
        if (!pt) return std::nullopt;

        const Point &pt_ref = *pt;

        // p173 admissibility — x0, xm lie on the bounding chords.
        // [CY86] p174: "A(C,C') is the area of the unique balanced (C,C')-chain,
        // if it is ADMISSIBLE. Otherwise A(C,C') = -infinity." An inadmissible
        // chain is therefore absent, not stored.
        if (!on_chord_segment(pt_ref, c)) return std::nullopt;
        if (!on_chord_segment(pt_ref, c_prime)) return std::nullopt;

        const Point c_ref = reference_point(c, poly);
        const Point c_prime_ref = reference_point(c_prime, poly);

        // The area of a (C, C')-chain is AREA((c, x_0, x_1, ..., x_m, c')), where c
        // and c' are the reference points of C and C' and the x_i are its nodes.
        entry.area = triangle_area_signed(c_ref, pt_ref, c_prime_ref);

        entry.x0 = pt_ref;
        entry.xm = pt_ref;

        entry.chord_c_idx = a;  // store C, C' indices to recover L0, L1
        entry.chord_c_prime_idx = b;

        return entry;
    }

    std::optional<std::pair<Point, Point>>
    find_balanced_chord(const Point &p,  // reflex corner — pivot
                        const Line &L0,  // line through bounding chord C
                        const Line &L1)  // line through bounding chord C'
    {
        // x1 = 2p - x0 must lie on L1
        // Parametrise L0: x0 = A + t*(B - A)
        // Then 2p - x0 must satisfy L1's equation → solve for t

        // As x0 moves along L0, x1 = 2p - x0 traces L0 reflected through p,
        // the line through 2p - A and 2p - B for two points A, B on L0.

        const Point A = L0.point(0);
        const Point B = L0.point(1);

        const auto A_ref = Point(2 * p.x() - A.x(), 2 * p.y() - A.y());
        const auto B_ref = Point(2 * p.x() - B.x(), 2 * p.y() - B.y());

        const Line L0_reflected(A_ref, B_ref);

        // x1 = intersection of L0_reflected and L1
        const auto result = CGAL::intersection(L0_reflected, L1);
        if (!result) return std::nullopt;  // parallel lines

        const Point *x1 = std::get_if<Point>(&*result);
        if (!x1) return std::nullopt;  // coincident lines

        // x0 = 2p - x1
        Point x0(2 * p.x() - x1->x(), 2 * p.y() - x1->y());

        return std::make_pair(x0, *x1);
    }

    FT length_1_chain_area(const Point &c_ref,  // reference point of C
                           const Point &x0, const Point &x1,
                           const Point &c_prime_ref)  // reference point of C'
    {
        // Area of quadrilateral (c, x0, x1, c')
        // Split into two triangles
        return triangle_area_signed(c_ref, x0, x1) + triangle_area_signed(c_ref, x1, c_prime_ref);
    }

    std::optional<ChainEntry> try_length_1_chain(const std::vector<std::size_t> &R,
                                                 const std::vector<Chord> &Xi, const std::size_t a,
                                                 const std::size_t b, const Polygon &poly,
                                                 const std::vector<bool> &reflex_flags) {
        const Chord &c = Xi[a];
        const Chord &c_prime = Xi[b];

        const Line L0(poly[c.i], poly[c.j]);
        const Line L1(poly[c_prime.i], poly[c_prime.j]);
        const Point c_ref = reference_point(c, poly);
        const Point c_prime_ref = reference_point(c_prime, poly);

        std::optional<ChainEntry> best = std::nullopt;

        for (const std::size_t idx : R) {
            if (!reflex_flags[idx]) {
                continue;
            }
            const Point &p = poly[idx];

            // The pivot must lie strictly inside C+ and C'+ (§4, open ρ); a pivot
            // on a bounding line is rejected here.
            if (!in_C_plus(p, c, poly, reflex_flags)) {
                continue;
            }
            if (!in_C_plus(p, c_prime, poly, reflex_flags)) {
                continue;
            }

            auto chord = find_balanced_chord(p, L0, L1);
            if (!chord) {
                continue;
            }
            auto [x0, xm] = *chord;

            // Degenerate constructions are not chords (a chord has positive
            // length and crosses both lines): x0 == xm, or a chord collinear
            // with a bounding line.
            if (x0 == xm) {
                continue;
            }
            if (CGAL::orientation(poly[c.i], poly[c.j], xm) == CGAL::COLLINEAR
                || CGAL::orientation(poly[c_prime.i], poly[c_prime.j], x0) == CGAL::COLLINEAR) {
                continue;
            }

            // p173 admissibility, per [CY86] p174's -infinity rule.
            if (!on_chord_segment(x0, c) || !on_chord_segment(xm, c_prime)) {
                continue;
            }

            // The core (finite-region) side of the chord. The finite region of
            // the wedge L0+ ∩ L1+ is on the side of the apex L0 ∩ L1. The
            // reference points are only a fallback, when the apex is collinear
            // with the chord or L0 ∥ L1: c_ref is a valid witness only when x0
            // does not separate it from the apex.
            CGAL::Orientation core_side = CGAL::COLLINEAR;
            if (const auto apex = CGAL::intersection(L0, L1))
                if (const Point *o = std::get_if<Point>(&*apex))
                    core_side = CGAL::orientation(x0, xm, *o);
            if (core_side == CGAL::COLLINEAR) core_side = CGAL::orientation(x0, xm, c_ref);
            if (core_side == CGAL::COLLINEAR) core_side = CGAL::orientation(x0, xm, c_prime_ref);
            if (core_side == CGAL::COLLINEAR) {
                continue;
            }

            bool chain_is_empty = true;
            for (const std::size_t r_idx : R) {
                const Point &r = poly[r_idx];
                if (r == p) continue;

                // §4 open ρ: a corner on a bounding line never obstructs, which
                // the strict in_C_plus handles.
                if (!in_C_plus(r, c, poly, reflex_flags)
                    || !in_C_plus(r, c_prime, poly, reflex_flags))
                    continue;
                const auto r_side = CGAL::orientation(x0, xm, r);
                if (r_side != CGAL::COLLINEAR && r_side != core_side) {
                    chain_is_empty = false;
                    break;
                }
            }
            // Emptiness is what Theorem 5's uniqueness rests on.
            if (!chain_is_empty) {
                continue;
            }

            if (!is_contained(poly, Segment(x0, xm))) {
                continue;
            }

            // Signed, as in the paper.
            FT area = length_1_chain_area(c_ref, x0, xm, c_prime_ref);

            // Theorem 5 promises at most one survivor; keep the largest.
            if (!best.has_value() || area > best->area) {
                ChainEntry entry;
                entry.area = area;
                entry.x0 = x0;
                entry.xm = xm;
                entry.pivot_indices = {idx};
                entry.chord_c_idx = a;
                entry.chord_c_prime_idx = b;
                best = entry;
            }
        }

        return best;
    }

#pragma endregion

#pragma region Step 6
    // ordering check "x'_m', p, q, x''_0 occur in
    // that order in L2" — Section 5 case A
    bool ordered_along_line(const Point &a, const Point &b, const Point &c, const Point &d) {
        // All four points lie on L2 (the line of C'') by construction.
        return CGAL::are_ordered_along_line(a, b, c) && CGAL::are_ordered_along_line(b, c, d);
    }

    // Case B version — five points
    // "x'_m', p, q, x''_0, x'_{m'+2} occur in
    // that order on L2" — Section 5 case B
    bool ordered_along_line_5(const Point &a, const Point &b, const Point &c, const Point &d,
                              const Point &e) {
        return CGAL::are_ordered_along_line(a, b, c) && CGAL::are_ordered_along_line(b, c, d)
               && CGAL::are_ordered_along_line(c, d, e);
    }

    // "verify that the two chains do not intersect"
    // — Section 5 cases (A) and (B), Fig. 16.
    //
    // `left` is the (C, C_mid)-chain and `right` the (C_mid, C')-chain. They share
    // no node, so any intersection between their polylines is a crossing.
    // O(m'·m'').
    bool chains_do_not_intersect(const ChainEntry &left,  // (C, C_mid)-chain
                                 const ChainEntry &right, const Polygon &poly) {
        // (C_mid, C')-chain — L2/poly implicit via nodes
        const std::vector<Point> ln = reconstruct_chain_nodes(left, poly);
        const std::vector<Point> rn = reconstruct_chain_nodes(right, poly);

        for (std::size_t i = 0; i + 1 < ln.size(); ++i) {
            const Segment s1(ln[i], ln[i + 1]);
            for (std::size_t j = 0; j + 1 < rn.size(); ++j)
                if (CGAL::do_intersect(s1, Segment(rn[j], rn[j + 1]))) return false;
        }
        return true;
    }

    // Section 5 case (A) junction construction
    std::optional<Point>
    compute_junction_node_A(const Point &x_prime_next,    // x'_{m'+1} = x'_{m'} | p_{m'+1}
                            const Point &x_double_first,  // x''_0
                            const Line &L0, const Line &L1) {
        const Vector L0_dir = L0.to_vector();
        const Vector L1_dir = L1.to_vector();

        // L_0(x'_{m'+1}) — Section 3 notation
        const Line L0_parallel(x_prime_next, x_prime_next + L0_dir);

        // L_1(x''_0) — Section 3 notation
        const Line L1_parallel(x_double_first, x_double_first + L1_dir);

        // x_{m'+1} = L_1(x''_0) ∩ L_0(x'_{m'+1})
        // — Section 5
        const auto inter = CGAL::intersection(L0_parallel, L1_parallel);
        if (!inter) return std::nullopt;
        const Point *pt = std::get_if<Point>(&*inter);
        if (!pt) return std::nullopt;
        return *pt;
    }

    // "Set x'_{m'+1} = x'_m|p and x'_{m'+2} =
    // x'_{m'+1}|q" — Section 5 case B
    // "x|p" is the reflection of x through p (Section 3).
    std::pair<Point, Point> compute_junction_nodes_B(const Point &x_prime_last,  // x'_{m'}
                                                     const Point &pivot_p,       // p
                                                     const Point &pivot_q)       // q
    {
        // x'_{m'+1} = reflection of x'_{m'} through p
        Point x_next(2 * pivot_p.x() - x_prime_last.x(), 2 * pivot_p.y() - x_prime_last.y());

        // x'_{m'+2} = reflection of x'_{m'+1} through q
        Point x_next2(2 * pivot_q.x() - x_next.x(), 2 * pivot_q.y() - x_next.y());

        return {x_next, x_next2};
    }

    // Memoised path data for the two composition-compatibility checks (Lemma 8,
    // p173), keyed by chord indices that fully determine (L, pivots, R):
    //   L0 side: (left_idx, mid_idx, variant)  variant 0=+p, 1=+q (case A), 2=+p,q (case B)
    //   L1 side: (mid_idx, right_idx)           reversed right pivots (shared by A and B)
    const CritEntry &cached_l0_entry(CritParamCache &cache, const std::size_t left_idx,
                                     const std::size_t mid_idx, const int variant, const Line &,
                                     const std::function<CritEntry()> &build) {
        if (auto *hit = cache.get_l0(left_idx, mid_idx, variant)) {
            return *hit;
        }
        return cache.set_l0(left_idx, mid_idx, variant, build());
    }

    const CritEntry &cached_l1_entry(CritParamCache &cache, const std::size_t mid_idx,
                                     const std::size_t right_idx, const Line &,
                                     const std::function<CritEntry()> &build) {
        if (auto *hit = cache.get_l1(mid_idx, right_idx)) {
            return *hit;
        }
        return cache.set_l1(mid_idx, right_idx, build());
    }

    VerificationResult verify_case_A_cached(const ChainEntry &left, const ChainEntry &right,
                                            const Chord &C_mid, const Chord &C_left,
                                            const Chord &C_right, const Polygon &poly,
                                            CritParamCache &pcache, const std::size_t left_idx,
                                            const std::size_t mid_idx,
                                            const std::size_t right_idx) {
        const Point &p = poly[C_mid.i];
        const Point &q = poly[C_mid.j];
        const Line L0(poly[C_left.i], poly[C_left.j]);
        const Line L1(poly[C_right.i], poly[C_right.j]);
        const Line L2(p, q);

        // CHECK 1 — ordering x'_{m'}, p, q, x''_0 on L_2
        // Section 5 case (A)
        if (!ordered_along_line(left.xm, p, q, right.x0)) {
            return {false, FT(0)};
        }

        // CHECK 2 — chains do not intersect
        // Section 5
        if (!chains_do_not_intersect(left, right, poly)) {
            return {false, FT(0)};
        }

        // p_{m'+1} = p if m'' even, q if m'' odd
        // — Section 5 case (A) pivot assignment
        const bool m_double_prime_even = (right.pivot_indices.size() % 2 == 0);
        const Point &junction_pivot = m_double_prime_even ? p : q;

        // x'_{m'+1} = x'_{m'} | p_{m'+1}
        // — Section 5
        const Point x_prime_next(2 * junction_pivot.x() - left.xm.x(),
                                 2 * junction_pivot.y() - left.xm.y());

        // x_{m'+1} = L_1(x''_0) ∩ L_0(x'_{m'+1})
        // — Section 5
        const auto junction = compute_junction_node_A(x_prime_next, right.x0, L0, L1);
        if (!junction) {
            return {false, FT(0)};
        }

        // CHECK 3 — (L_0, p'_1, ..., p'_{m'+1})-equivalence
        // Paper: "Verify x_{m'+1} and x'_{m'+1} are (L_0, p'_1..p'_{m'+1})-
        // equivalent." Section 5 also notes the path must be extended to include
        // one of the pivots on L_1 (the junction pivot) — done here.
        {
            // variant 0 = extended by p (m'' even), 1 = extended by q (m'' odd)
            const int variant = m_double_prime_even ? 0 : 1;
            const std::size_t m_ext = left.pivot_indices.size() + 1;

            // Equivalent iff no critical point lies strictly between the two
            // nodes, using the memoised path data for (left, mid, variant).
            const CritEntry &ce = cached_l0_entry(pcache, left_idx, mid_idx, variant, L0, [&] {
                std::vector<Point> extended_pivots = pivots_of(left.pivot_indices, poly);
                extended_pivots.push_back(junction_pivot);
                const auto R_left = compute_R(C_left, C_mid, poly.size());
                return CritEntry{compute_critical_params(L0, extended_pivots, R_left, poly),
                                 last_lever_offset(extended_pivots, L0.point(0))};
            });
            if (!equivalent_given_params(L0, m_ext, ce.last_offset, ce.params, *junction,
                                         x_prime_next)) {
                return {false, FT(0)};
            }
        }

        // CHECK 4 — equivalence on L_1 with the right pivots reversed. Case (A)'s
        // text has them forward, an erratum: the path walks back from L_1, and
        // case (B) writes "(L_1, p''_{m''}..p''_1)".
        {
            // Memoised by (mid, right), shared with case B's CHECK 4.
            const CritEntry &ce = cached_l1_entry(pcache, mid_idx, right_idx, L1, [&] {
                const std::vector<Point> right_pts = pivots_of(right.pivot_indices, poly);
                const std::vector<Point> reversed_pivots(right_pts.rbegin(), right_pts.rend());
                const auto R_right = compute_R(C_mid, C_right, poly.size());
                return CritEntry{compute_critical_params(L1, reversed_pivots, R_right, poly),
                                 last_lever_offset(reversed_pivots, L1.point(0))};
            });
            if (!equivalent_given_params(L1, right.pivot_indices.size(), ce.last_offset, ce.params,
                                         *junction, right.x0)) {
                return {false, FT(0)};
            }
        }

        // composed area = left + right
        // — triangle term added in Step II — Section 6
        return {true, left.area + right.area, CompositionCase::CASE_A, junction_pivot};
    }

    // case B verification — Section 5
    VerificationResult verify_case_B_cached(const ChainEntry &left, const ChainEntry &right,
                                            const Chord &C_mid, const Chord &C_left,
                                            const Chord &C_right, const Polygon &poly,
                                            CritParamCache &pcache, const std::size_t left_idx,
                                            const std::size_t mid_idx,
                                            const std::size_t right_idx) {
        const Point &p = poly[C_mid.i];
        const Point &q = poly[C_mid.j];

        const Line L0(poly[C_left.i], poly[C_left.j]);
        const Line L1(poly[C_right.i], poly[C_right.j]);

        // Compute junction nodes for case B
        // "Set x'_{m'+1} = x'_m|p and
        // x'_{m'+2} = x'_{m'+1}|q" — Section 5
        auto [x_next, x_next2] = compute_junction_nodes_B(left.xm, p, q);

        // CHECK 1 — chains do not intersect
        // "verify that 𝒞' and 𝒞'' do not
        // intersect" — Section 5
        if (!chains_do_not_intersect(left, right, poly)) {
            return {false, FT(0)};
        }

        // CHECK 2 — ordering of five points on L2
        // "x'_{m'}, p, q, x''_0, x'_{m'+2}
        // occur in that order on L2" — Section 5
        if (!ordered_along_line_5(left.xm, p, q, right.x0, x_next2)) {
            return {false, FT(0)};
        }

        // x_{m'+2} = L_1(x''_0) ∩ L_0(x'_{m'+2}), as in case A's junction.
        const Vector L0_dir = L0.to_vector();
        const Vector L1_dir = L1.to_vector();

        const Line L0_through(x_next2, x_next2 + L0_dir);
        const Line L1_through(right.x0, right.x0 + L1_dir);

        const auto inter = CGAL::intersection(L0_through, L1_through);
        if (!inter) {
            return {false, FT(0)};
        }
        const Point *junction = std::get_if<Point>(&*inter);
        if (!junction) {
            return {false, FT(0)};
        }

        // CHECK 3 — equivalence on L0, on the path EXTENDED by the seam
        // pivots. The paper writes "(L0, p'_1..p'_{m'})-equivalent", but the
        // points x_{m'+2} and x'_{m'+2} lie on no node line of that path (each
        // seam reflection changes the distance from L0); the check is only
        // well-defined on the (L0, p'_1..p'_{m'}, p, q)-path, whose node-(m'+2)
        // line contains both — the case-B analogue of case A's junction-pivot
        // extension (and of the Section 5 note after Lemma 8).
        {
            // variant 2 = case B's double extension (+p, +q).
            const std::size_t m_ext = left.pivot_indices.size() + 2;
            const CritEntry &ce = cached_l0_entry(pcache, left_idx, mid_idx, 2, L0, [&] {
                std::vector<Point> extended_pivots = pivots_of(left.pivot_indices, poly);
                extended_pivots.push_back(p);
                extended_pivots.push_back(q);
                const auto R_left = compute_R(C_left, C_mid, poly.size());
                return CritEntry{compute_critical_params(L0, extended_pivots, R_left, poly),
                                 last_lever_offset(extended_pivots, L0.point(0))};
            });
            if (!equivalent_given_params(L0, m_ext, ce.last_offset, ce.params, *junction,
                                         x_next2)) {
                return {false, FT(0)};
            }
        }

        // CHECK 4 — (L_1, p''_{m''}, ..., p''_1)-equivalence
        // Paper: reversed pivot order for case (B), with L_1 as starting line.
        {
            // Same memoised (mid, right) data as case A's CHECK 4.
            const CritEntry &ce = cached_l1_entry(pcache, mid_idx, right_idx, L1, [&] {
                const std::vector<Point> right_pts_b = pivots_of(right.pivot_indices, poly);
                const std::vector<Point> reversed_pivots(right_pts_b.rbegin(), right_pts_b.rend());
                const auto R_right = compute_R(C_mid, C_right, poly.size());
                return CritEntry{compute_critical_params(L1, reversed_pivots, R_right, poly),
                                 last_lever_offset(reversed_pivots, L1.point(0))};
            });
            if (!equivalent_given_params(L1, right.pivot_indices.size(), ce.last_offset, ce.params,
                                         *junction, right.x0)) {
                return {false, FT(0)};
            }
        }

        return {true, left.area + right.area, CompositionCase::CASE_B, Point()};
    }

    // Lemma 8 — O(1) verification checking both
    // case A and case B — Section 5
    VerificationResult
    verify_decomposition_cached(const ChainEntry &left, const ChainEntry &right, const Chord &C_mid,
                                const Chord &C_left, const Chord &C_right, const Polygon &poly,
                                CritParamCache &pcache, const std::size_t left_idx,
                                const std::size_t mid_idx, const std::size_t right_idx) {
        if (auto result_A = verify_case_A_cached(left, right, C_mid, C_left, C_right, poly, pcache,
                                                 left_idx, mid_idx, right_idx);
            result_A.valid)
            return result_A;

        if (auto result_B = verify_case_B_cached(left, right, C_mid, C_left, C_right, poly, pcache,
                                                 left_idx, mid_idx, right_idx);
            result_B.valid)
            return result_B;

        return {false, FT(0), CompositionCase::CASE_A, Point()};
    }

    // A composed chain in two halves: compose_core gives the pivots and endpoints,
    // finish_composed the area. p173's endpoint gate needs only x0 and xm and
    // rejects most compositions, so the area is computed only for survivors.
    namespace {
        struct ComposedCore {
            std::vector<std::size_t> pivot_indices;
            Point x0, xm;
            bool solved{};  // false ⇒ the L0 ∥ L1 degenerate fallback was taken
        };

        ComposedCore compose_core(const ChainEntry &left, const ChainEntry &right, const Chord &C,
                                  const Chord &C_mid, const Chord &C_prime, const Polygon &poly,
                                  const VerificationResult &vr) {
            const Line L0(poly[C.i], poly[C.j]);
            const Line L1(poly[C_prime.i], poly[C_prime.j]);

            ComposedCore core;

            // x0 and xm are not inherited from the sub-chains: the junction node
            // is not a reflection of left.xm, so left.x0 does not parametrise the
            // composed chain. They are solved from (L0, L1, all pivots) below.
            core.pivot_indices = left.pivot_indices;

            if (vr.case_type == CompositionCase::CASE_A) {
                // case A inserts one junction pivot
                // p_{m'+1} = p if m'' even, q if m'' odd — Section 5
                const bool m_double_prime_even = (right.pivot_indices.size() % 2 == 0);
                core.pivot_indices.push_back(m_double_prime_even ? C_mid.i : C_mid.j);
            } else {
                // CASE_B
                // case B inserts two junction pivots
                // p_{m'+1} = p, p_{m'+2} = q always — Section 5
                core.pivot_indices.push_back(C_mid.i);
                core.pivot_indices.push_back(C_mid.j);
            }

            core.pivot_indices.insert(core.pivot_indices.end(), right.pivot_indices.begin(),
                                      right.pivot_indices.end());

            // The composed chain's x0 from its defining data (L0, L1, all pivots).
            // Theorem 5 makes this chain unique,
            // and Lemma 8 bounds reconstruction at O(n). Then derive xm by the lever
            // rule from that x0 so x0, the pivots and xm are mutually consistent.
            const std::vector<Point> composed_pivots = pivots_of(core.pivot_indices, poly);
            if (const auto x0_solved = solve_balanced_chain_x0(L0, L1, composed_pivots)) {
                core.solved = true;
                core.x0 = *x0_solved;
                Point prev = core.x0;
                for (const auto &pivot : composed_pivots)
                    prev = Point(2 * pivot.x() - prev.x(), 2 * pivot.y() - prev.y());
                core.xm = prev;
            } else {
                // L0 ∥ L1: the lever endpoint constraint is degenerate (no unique
                // x0). Fall back to the sub-chain endpoints; the seam-sum area is
                // applied by finish_composed, there being no consistent node set to
                // integrate over.
                core.solved = false;
                core.x0 = left.x0;
                core.xm = right.xm;
            }
            return core;
        }

        // The p174 area over the core's own nodes, plus the chord indices.
        ChainEntry finish_composed(const ComposedCore &core, const ChainEntry &left,
                                   const ChainEntry &right, const Chord &C, const Chord &C_mid,
                                   const Chord &C_prime, const Polygon &poly) {
            const Point c = reference_point(C, poly);
            const Point c_mid = reference_point(C_mid, poly);
            const Point c_prime = reference_point(C_prime, poly);

            ChainEntry composed;
            composed.pivot_indices = core.pivot_indices;
            composed.x0 = core.x0;
            composed.xm = core.xm;

            // The area is computed from the composed chain's own nodes, not by the
            // paper's O(1) shortcut left + right + AREA(△cc''c'). The shortcut is
            // exact only when the sub-chains are lever-equivalent to the composed
            // one, which the equivalence checks can accept without. A(C,C') is
            // defined as the area of the chain (p174), so this is correct by
            // construction, and the equivalence checks only gate existence.
            if (core.solved) {
                const std::vector<Point> composed_pivots = pivots_of(core.pivot_indices, poly);
                FT area = FT(0);
                Point prev = core.x0;
                for (const auto &pivot : composed_pivots) {
                    const Point next(2 * pivot.x() - prev.x(), 2 * pivot.y() - prev.y());
                    area += triangle_area_signed(c, prev, next);
                    prev = next;
                }
                area += triangle_area_signed(c, core.xm, c_prime);
                composed.area = area;
            } else {
                composed.area = left.area + right.area + triangle_area_signed(c, c_mid, c_prime);
            }

            // Section 5's extended-path critical intervals are computed where they
            // are used, in verify_case_A/B, not stored here.

            // chord indices stored for future
            // line reconstruction — outermost bounding chords
            composed.chord_c_idx = left.chord_c_idx;
            composed.chord_c_prime_idx = right.chord_c_prime_idx;
            return composed;
        }
    }  // namespace

    // O(n) stages × O(n^6) triples = O(n^7)
    // total — Section 6
    void compute_A_matrix_all_stages(AMatrix &A, const std::vector<Chord> &Xi, const Polygon &poly,
                                     const std::vector<bool> &reflex_flags,
                                     const std::size_t n_vertices)  // vertex count — stage bound
    {
        const std::size_t N = Xi.size();

        // Memoized critical parameters for the Lemma 8 equivalence checks,
        // shared across all stages. Keyed by chord/seam indices; valid across
        // stages because each A-entry is immutable once set (Theorem 5).
        CritParamCache pcache;

        // Sparse adjacency over A: row[a] lists c with A(a,c) present, col[b] the c
        // with A(c,b), both ascending. Walking the shorter and probing the other
        // gives the mids in ascending order, which matters: composition commits
        // to the first valid decomposition. Updated as entries are added.
        std::vector<std::vector<std::size_t>> row(N), col(N);
        for (const auto &kv : A.data) {
            row[kv.first.first].push_back(kv.first.second);
            col[kv.first.second].push_back(kv.first.first);
        }
        for (auto &v : row) std::sort(v.begin(), v.end());
        for (auto &v : col) std::sort(v.begin(), v.end());

        // O(n) stages — Section 6
        for (std::size_t stage = 1; stage <= n_vertices; stage++) {
            bool any_new = false;

            // iterate all triples (C, C'', C')
            // — Section 6
            for (std::size_t a = 0; a < N; a++) {
                for (std::size_t b = 0; b < N; b++) {
                    if (!is_admissible_pair(Xi[a], Xi[b], poly.size())) continue;

                    // Theorem 5: at most one valid chain per pair, so a filled
                    // slot is final. CritParamCache also relies on this.
                    if (A.is_valid(a, b)) continue;

                    // Walk the shorter side; probe the other. Same mids, same order.
                    const bool by_row = row[a].size() <= col[b].size();
                    const std::vector<std::size_t> &mids = by_row ? row[a] : col[b];
                    bool stored_here = false;
                    for (std::size_t mi = 0; mi < mids.size(); ++mi) {
                        const std::size_t c = mids[mi];
                        if (by_row) {
                            if (!A.is_valid(c, b)) continue;
                        } else {
                            if (!A.is_valid(a, c)) continue;
                        }

                        const ChainEntry &left = A.get_ref(a, c);
                        const ChainEntry &right = A.get_ref(c, b);

                        // Lemma 8 verification
                        // — Section 5. Indices: left=a, mid=c, right=b
                        // (C_mid=Xi[c], C_left=Xi[a], C_right=Xi[b]).
                        VerificationResult vr;
                        vr = verify_decomposition_cached(left, right, Xi[c], Xi[a], Xi[b], poly,
                                                         pcache, a, c, b);

                        if (!vr.valid) continue;

                        // "If so, construct the (C,C')-chain" (Section 6): build,
                        // gate and store one verification result. Endpoints first,
                        // area only if the endpoint gate passes.
                        const auto attempt = [&](const VerificationResult &v) -> bool {
                            ComposedCore core;
                            core = compose_core(left, right, Xi[a], Xi[c], Xi[b], poly, v);

                            // p173 admissibility for the COMPOSED (C,C')-chain: the
                            // first/last nodes x0, xm lie on the bounding chords.
                            if (!on_chord_segment(core.x0, Xi[a])) {
                                return false;
                            }
                            if (!on_chord_segment(core.xm, Xi[b])) {
                                return false;
                            }

                            const ChainEntry cc
                                = finish_composed(core, left, right, Xi[a], Xi[c], Xi[b], poly);

                            // §2/§6: every chord is a segment in P. A re-solved
                            // interior chord can cut across a concavity and leave P,
                            // which the corner-based emptiness test does not see.
                            {
                                const std::vector<Point> ck_nodes
                                    = reconstruct_chain_nodes(cc, poly);
                                for (std::size_t t = 0; t + 1 < ck_nodes.size(); ++t)
                                    if (!is_contained(poly,
                                                      Segment(ck_nodes[t], ck_nodes[t + 1]))) {
                                        return false;
                                    }
                            }

                            A.set(a, b, cc);
                            any_new = true;
                            stored_here = true;
                            return true;
                        };

                        bool stored_now = attempt(vr);

                        // Theorem 5 uniqueness —
                        // at most one valid composition
                        if (stored_now) break;
                    }

                    // Keep the adjacency in step with A, sorted. Outside the loop
                    // because `mids` aliases one of these vectors.
                    if (stored_here) {
                        row[a].insert(std::upper_bound(row[a].begin(), row[a].end(), b), b);
                        col[b].insert(std::upper_bound(col[b].begin(), col[b].end(), a), a);
                    }
                }
            }

            // monotonicity — no new sub-chains
            // means no future compositions possible
            if (!any_new) break;
        }

        // Report CritParamCache's size here, before it is freed. It is counted in
        // FTs rather than measured by RSS, which does not drop when it is freed.
        if (memprobe::enabled()) {
            std::size_t fts = 0;
            // +2 per entry for the cached A[m] Vector's two FTs.
            for (const auto &kv : pcache.l0) fts += kv.second.params.size() + 2;
            for (const auto &kv : pcache.l1) fts += kv.second.params.size() + 2;
            memprobe::note_crit(pcache.l0.size(), pcache.l1.size(), fts);
        }
    }
#pragma endregion

#pragma region Step 8
    // ── Chain-definition guards ──────────────────────────────────────────────
    // Section 3 defines a chain: chords whose consecutive intersections (the
    // nodes) form a convex polygon core, contained in the bounding wedge
    // L0+ ∩ L1+ (Section 4). Steps II/III say maximisation needs no
    // compatibility check (p175/p176), but with maximal-chord admissibility,
    // cells can assemble into multiply wound reference polygons whose signed
    // area exceeds the true optimum. These guards enforce the chain definition
    // compositionally: they are the conditions under which the signed-area sum
    // equals a clipped area.

    Vector chord_dir(const Chord &c, const Polygon &poly) {
        return poly[c.j] - poly[c.i];
    }

    // Left turn or collinear (CCW convention).
    bool left_turn_ok(const Vector &u, const Vector &v) {
        return u.x() * v.y() - u.y() * v.x() >= FT(0);
    }

    // Section 4: chains live in the bounding wedge. Under the CCW convention
    // the chain side of directed chord C_{i,j} is its LEFT closed half-plane.
    bool in_chain_wedge(const Point &x, const Chord &ca, const Chord &cb, const Polygon &poly) {
        return CGAL::orientation(poly[ca.i], poly[ca.j], x) != CGAL::RIGHT_TURN
               && CGAL::orientation(poly[cb.i], poly[cb.j], x) != CGAL::RIGHT_TURN;
    }

    // Section 3: the core is a convex CLOSED polygon, so the closure turns at
    // both endpoints must also be left turns. (Any contiguous subsequence of
    // a convex polygon's vertices is itself convex, so requiring this of every
    // sub-chain never rejects a piece of a genuine potato.)
    bool closed_core_ok(const Point &x_first, const Vector &d_first, const Point &x_last,
                        const Vector &d_last) {
        if (x_first == x_last) return true;
        const Vector close = x_first - x_last;
        return left_turn_ok(d_last, close) && left_turn_ok(close, d_first);
    }

    // Parameter of x along chord c's corner pair: t(v_i) = 0, t(v_j) = 1.
    FT chord_param(const Point &x, const Chord &c, const Polygon &poly) {
        const Vector d = poly[c.j] - poly[c.i];
        return ((x - poly[c.i]) * d) / (d * d);
    }

    // p174 compatibility at a seam C'', with Section-3 pivots (a single pivot,
    // p = q, is allowed — the truncated seam segment must contain at least one
    // of C'''s corners, not necessarily both): the traversal moves FORWARD
    // along the seam, and the turns onto and off the seam chord are left turns.
    bool seam_ok(const MEntry &left, const Chord &mid, const MEntry &right, const Polygon &poly) {
        const FT tm = chord_param(left.x_last, mid, poly);
        const FT t0 = chord_param(right.x_first, mid, poly);
        if (tm > t0) return false;                   // backward along the seam
        if (t0 < FT(0) || tm > FT(1)) return false;  // no seam corner inside
        const Vector ds = chord_dir(mid, poly);
        return left_turn_ok(left.d_last, ds) && left_turn_ok(ds, right.d_first);
    }

    // Build the M-leaf state for an A entry, or nullopt if the entry violates
    // the Section-3 chain definition: left turns at every node (in TRAVERSAL
    // directions — not reference-point wedges, which diverge from the true
    // turn for far nodes), a convex closed CCW core, no closed-loop
    // degeneracy (x0 == xm with m >= 1 is a fictional composed entry), and
    // wedge containment of every node.
    std::optional<MEntry> m_leaf_from_A(const ChainEntry &ae, const Chord &ca, const Chord &cb,
                                        const Polygon &poly) {
        const Vector da = chord_dir(ca, poly);
        const Vector db = chord_dir(cb, poly);

        if (ae.pivot_indices.empty()) {
            if (!left_turn_ok(da, db)) return std::nullopt;
            if (!in_chain_wedge(ae.x0, ca, cb, poly)) return std::nullopt;
            return MEntry{ae.area, std::nullopt, ae.x0, ae.xm, db, da};
        }

        if (ae.x0 == ae.xm) return std::nullopt;  // closed-loop phantom

        const auto nodes = reconstruct_chain_nodes(ae, poly);
        if (nodes.size() >= 3) {
            FT shoelace(0);
            for (std::size_t t = 1; t + 1 < nodes.size(); ++t)
                shoelace += triangle_area_signed(nodes[0], nodes[t], nodes[t + 1]);
            if (shoelace < FT(0)) return std::nullopt;  // core not CCW
        }

        const Vector d_first = nodes[1] - nodes[0];
        const Vector d_last = nodes[nodes.size() - 1] - nodes[nodes.size() - 2];
        if (!left_turn_ok(da, d_first)) return std::nullopt;
        if (!left_turn_ok(d_last, db)) return std::nullopt;
        for (std::size_t t = 1; t + 1 < nodes.size(); ++t)
            if (!left_turn_ok(nodes[t] - nodes[t - 1], nodes[t + 1] - nodes[t]))
                return std::nullopt;
        if (!closed_core_ok(nodes.front(), d_first, nodes.back(), d_last)) return std::nullopt;
        for (const auto &nd : nodes)
            if (!in_chain_wedge(nd, ca, cb, poly)) return std::nullopt;

        return MEntry{ae.area, std::nullopt, ae.x0, ae.xm, d_first, d_last};
    }

    // Step II, stage 1 (Section 6). Areas are signed (2·AREA = Σ x_i(y_{i+1} - y_i),
    // positive CCW); the recursion depends on the sign.
    void compute_M_stage_1(MMatrix &M, const std::vector<Chord> &Xi,
                           const std::vector<std::vector<std::size_t>> &Xi_by_vertex,
                           const Polygon &poly, const std::vector<bool> &reflex_flags) {
        const std::size_t n = poly.size();

        for (std::size_t i = 0; i < n; i++) {
            const std::size_t j = (i + 1) % n;

            for (const std::size_t a : Xi_by_vertex[i]) {
                if (Xi[a].j != j) continue;

                const Point c = reference_point(Xi[a], poly);
                const Point &vj = poly[j];
                for (const std::size_t b : Xi_by_vertex[j]) {
                    if (!is_admissible_pair(Xi[a], Xi[b], n)) continue;
                    // Stage 1: M(C,C') = AREA(△ c v_j c'). The convex core needs a
                    // left turn at v_j from C into C'; the sign of the area is
                    // not that test, and may be negative.
                    if (!left_turn_ok(chord_dir(Xi[a], poly), chord_dir(Xi[b], poly))) continue;
                    const Point c_prime = reference_point(Xi[b], poly);
                    const FT area = triangle_area_signed(c, vj, c_prime);
                    MEntry entry{area, std::nullopt,           vj,
                                 vj,   chord_dir(Xi[b], poly), chord_dir(Xi[a], poly)};
                    M.set(a, b, std::move(entry));
                }
            }
        }
    }

#pragma endregion

#pragma region Step 9
    // recursive formula for stage s > 1:
    // M(C,C') = max{A(C,C'),
    //               max_{C'' ∈ Ξ_k, i<k<j}
    //               {M(C,C'') + M(C'',C') + AREA(△cc''c')}}
    // — Section 6
    void compute_M_stage_s(MMatrix &M, const AMatrix &A, const std::vector<Chord> &Xi,
                           const std::vector<std::vector<std::size_t>> &Xi_by_vertex,
                           const Polygon &poly,
                           const std::size_t s)  // current stage
    {
        const std::size_t n = poly.size();

        // "stage s processes pairs (C∈Ξ_i, C'∈Ξ_j)
        // where j-i = s(mod n)" — Section 6
        for (std::size_t i = 0; i < n; i++) {
            const std::size_t j = (i + s) % n;

            for (std::size_t a : Xi_by_vertex[i]) {
                for (std::size_t b : Xi_by_vertex[j]) {
                    // skip non-admissible pairs
                    // M(C,C') = -inf if not admissible
                    // — Section 6
                    if (!is_admissible_pair(Xi[a], Xi[b], n)) continue;

                    Point c = reference_point(Xi[a], poly);
                    Point cp = reference_point(Xi[b], poly);

                    std::optional<MEntry> best = std::nullopt;

                    // first candidate is A(C,C') — balanced chain
                    // area — Section 6. Gated by the Section-3 chain
                    // definition (m_leaf_from_A): turns, convex closed CCW
                    // core, wedge containment, no closed-loop phantom.
                    if (A.is_valid(a, b)) {
                        if (auto leaf = m_leaf_from_A(A.get(a, b), Xi[a], Xi[b], poly)) {
                            best = std::move(*leaf);
                        }
                    }

                    // iterate C'' ∈ Ξ_k for
                    // i < k < j — Section 6
                    // "i < k < j" is the cyclic range between the
                    // two chords (consistent with admissibility) — the paper
                    // indexes cyclically throughout.
                    for (std::size_t k = (i + 1) % n; k != j; k = (k + 1) % n) {
                        for (std::size_t c_idx : Xi_by_vertex[k]) {
                            // Both sub-entries must exist
                            if (!M.is_valid(a, c_idx)) continue;
                            if (!M.is_valid(c_idx, b)) continue;

                            MEntry m_left = M.get(a, c_idx);
                            MEntry m_right = M.get(c_idx, b);

                            // The chain-definition guards at the seam
                            // (see above): p174 compatibility ordering,
                            // Section 3 core turns and closure, Section 4
                            // wedge containment.
                            if (!seam_ok(m_left, Xi[c_idx], m_right, poly)) continue;
                            if (!closed_core_ok(m_left.x_first, m_left.d_first, m_right.x_last,
                                                m_right.d_last))
                                continue;
                            if (!in_chain_wedge(m_left.x_first, Xi[a], Xi[b], poly)
                                || !in_chain_wedge(m_left.x_last, Xi[a], Xi[b], poly)
                                || !in_chain_wedge(m_right.x_first, Xi[a], Xi[b], poly)
                                || !in_chain_wedge(m_right.x_last, Xi[a], Xi[b], poly))
                                continue;

                            // reference point c''
                            // of C'' used in triangle area
                            // — Section 6
                            Point cpp = reference_point(Xi[c_idx], poly);

                            // AREA(△cc''c')
                            // — Section 6 recursive formula
                            FT tri = triangle_area_signed(c, cpp, cp);

                            // No sign filter on tri or candidate: the
                            // recursion maximizes SIGNED values; the shoelace
                            // identity keeps the sum exact regardless of sign.
                            FT candidate = m_left.area + m_right.area + tri;

                            // take the maximum over all C'' — Section 6
                            if (!best.has_value() || candidate > best->area) {
                                best = MEntry{candidate,      c_idx,          m_left.x_first,
                                              m_right.x_last, m_left.d_first, m_right.d_last};
                            }
                        }
                    }

                    // Store best found value
                    if (best.has_value()) {
                        M.set(a, b, std::move(*best));
                    }
                }  // end C' loop
            }  // end C loop
        }  // end i loop
    }
#pragma endregion

#pragma region Step 10
    // A chain's nodes from x0 and the pivots by the lever rule x_i = 2 p_i - x_{i-1}
    // (p_i is the midpoint of [x_{i-1}, x_i], Section 3): Lemma 8's O(n)
    // reconstruction.
    std::vector<Point> reconstruct_chain_nodes(const ChainEntry &chain, const Polygon &poly) {
        std::vector<Point> nodes;
        nodes.push_back(chain.x0);

        Point prev = chain.x0;

        for (const std::size_t pidx : chain.pivot_indices) {
            const Point &pivot = poly[pidx];
            Point xi(2 * pivot.x() - prev.x(), 2 * pivot.y() - prev.y());
            nodes.push_back(xi);
            prev = xi;
        }

        // The last node must equal the stored xm exactly.
        if (!chain.pivot_indices.empty()) {
            assert(nodes.back() == chain.xm
                   && "Chain node reconstruction mismatch — "
                      "pivot sequence or x0 storage error");
        }

        return nodes;
    }

    // Clip a polygon by the closed half-plane left of A→B (C+ for a directed
    // chord, Section 2). Sutherland-Hodgman.
    Polygon clip_by_halfplane(const Polygon &poly, const Point &A, const Point &B) {
        if (poly.size() == 0) return poly;

        Polygon result;
        const std::size_t n = poly.size();

        for (std::size_t i = 0; i < n; i++) {
            const Point &curr = poly[i];
            const Point &next = poly[(i + 1) % n];

            // C+ = NOT RIGHT_TURN — Section 2
            const bool curr_in = CGAL::orientation(A, B, curr) != CGAL::RIGHT_TURN;
            const bool next_in = CGAL::orientation(A, B, next) != CGAL::RIGHT_TURN;

            if (curr_in) result.push_back(curr);

            // The edge crosses the clip line: add the intersection.
            if (curr_in != next_in) {
                const Line clip_line(A, B);
                const Line edge_line(curr, next);
                const auto inter = CGAL::intersection(clip_line, edge_line);
                if (inter) {
                    if (const Point *pt = std::get_if<Point>(&*inter)) result.push_back(*pt);
                }
            }
        }

        // Remove consecutive duplicate vertices.
        if (result.size() > 1) {
            Polygon cleaned;
            for (std::size_t i = 0; i < result.size(); i++) {
                if (result[i] != result[(i + 1) % result.size()]) cleaned.push_back(result[i]);
            }
            result = cleaned;
        }

        return result;
    }

    // Orient a chord through a reflex corner so C+ is on its left. Section 2:
    // "C+ is the closed half-plane such that for a sufficiently small disc D
    // centered at u, D ∩ C+ ⊆ P"
    void orient_chord_for_clipping(Point &p1, Point &p2, const std::size_t reflex_idx,
                                   const Polygon &poly) {
        const std::size_t prev_idx = (reflex_idx + poly.size() - 1) % poly.size();
        const std::size_t next_idx = (reflex_idx + 1) % poly.size();

        // Interior direction at the reflex corner: the negated sum of the edge
        // vectors from it.
        const Vector v_in(-(poly[prev_idx].x() + poly[next_idx].x() - 2 * poly[reflex_idx].x()),
                          -(poly[prev_idx].y() + poly[next_idx].y() - 2 * poly[reflex_idx].y()));

        const Point test_pt(poly[reflex_idx].x() + v_in.x(), poly[reflex_idx].y() + v_in.y());

        if (CGAL::orientation(p1, p2, test_pt) == CGAL::RIGHT_TURN) {
            std::swap(p1, p2);
        }
    }

    // Orient an extremal chord so C+ is on its left (Section 2). False for a
    // convex–convex diagonal, which has no C+.
    bool orient_extremal_chord(Point &p1, Point &p2, const Chord &chord, const Polygon &poly,
                               const std::vector<bool> &reflex_flags) {
        if (reflex_flags[chord.i]) {
            orient_chord_for_clipping(p1, p2, chord.i, poly);
        } else if (reflex_flags[chord.j]) {
            orient_chord_for_clipping(p1, p2, chord.j, poly);
        } else {
            // Section 2 defines C+ for a chord through a reflex corner. The other
            // legal case is an edge of P: its interior side is left of the CCW
            // boundary edge. Edges must be kept: in a nonconvex P their extended
            // half-planes clip away distant parts of P.
            const std::size_t n = poly.size();
            if ((chord.i + 1) % n == chord.j) {
                return true;  // CCW edge v_i -> v_j : interior LEFT, keep as-is
            }
            if ((chord.j + 1) % n == chord.i) {
                std::swap(p1, p2);  // CCW edge v_j -> v_i : orient p1->p2 as CCW
                return true;
            }
            // A genuine non-consecutive convex–convex diagonal has no
            // well-defined C+ and must never bound the potato.
            return false;
        }

        return true;
    }

    void collect_chords(const MMatrix &M, const AMatrix &A, const std::vector<Chord> &Xi,
                        const Polygon &poly, const std::vector<bool> &reflex_flags, std::size_t a,
                        std::size_t b, std::vector<std::pair<Point, Point>> &chords) {
        if (!M.is_valid(a, b)) return;

        const MEntry &entry = M.get(a, b);

        if (!entry.best_c_mid_idx.has_value()) {
            if (!A.is_valid(a, b)) return;
            const auto ce = A.get(a, b);
            const auto nodes = reconstruct_chain_nodes(ce, poly);

            for (std::size_t i = 0; i + 1 < nodes.size(); i++) {
                Point p1 = nodes[i];
                Point p2 = nodes[i + 1];
                orient_chord_for_clipping(p1, p2, ce.pivot_indices[i], poly);
                chords.emplace_back(p1, p2);
            }
            return;
        }

        const std::size_t mid = entry.best_c_mid_idx.value();

        Point p1 = poly[Xi[mid].i];
        Point p2 = poly[Xi[mid].j];
        if (orient_extremal_chord(p1, p2, Xi[mid], poly, reflex_flags)) chords.push_back({p1, p2});

        collect_chords(M, A, Xi, poly, reflex_flags, a, mid, chords);
        collect_chords(M, A, Xi, poly, reflex_flags, mid, b, chords);
    }

    // Q = P ∩ (∩ C_i+) — Section 2
    Polygon reconstruct_potato(const MMatrix &M, const AMatrix &A, const std::vector<Chord> &Xi,
                               const Polygon &poly, const std::vector<bool> &reflex_flags,
                               std::size_t a, std::size_t b) {
        std::vector<std::pair<Point, Point>> chords;

        collect_chords(M, A, Xi, poly, reflex_flags, a, b, chords);
        collect_chords(M, A, Xi, poly, reflex_flags, b, a, chords);

        // The two bounding chords.
        {
            Point p1 = poly[Xi[a].i];
            Point p2 = poly[Xi[a].j];
            if (orient_extremal_chord(p1, p2, Xi[a], poly, reflex_flags))
                chords.push_back({p1, p2});
        }
        {
            Point p1 = poly[Xi[b].i];
            Point p2 = poly[Xi[b].j];
            if (orient_extremal_chord(p1, p2, Xi[b], poly, reflex_flags))
                chords.push_back({p1, p2});
        }

        // Q = P ∩ (⋂ C_i+) — Section 2. Two stages:
        //
        // (1) K = ⋂ C_i+ is convex, so build it by clipping a large bounding box
        //     by each half-plane; the intermediate stays convex, so
        //     Sutherland-Hodgman is exact. Clipping the nonconvex P directly
        //     could split it into several pieces that S-H joins into a
        //     self-intersecting ring.
        FT minx = poly[0].x(), maxx = poly[0].x();
        FT miny = poly[0].y(), maxy = poly[0].y();
        for (const auto &v : poly) {
            if (v.x() < minx) minx = v.x();
            if (v.x() > maxx) maxx = v.x();
            if (v.y() < miny) miny = v.y();
            if (v.y() > maxy) maxy = v.y();
        }
        const FT pad = (maxx - minx) + (maxy - miny) + FT(1);
        Polygon K;  // CCW bounding box
        K.push_back(Point(minx - pad, miny - pad));
        K.push_back(Point(maxx + pad, miny - pad));
        K.push_back(Point(maxx + pad, maxy + pad));
        K.push_back(Point(minx - pad, maxy + pad));
        for (const auto &[p1, p2] : chords) {
            K = clip_by_halfplane(K, p1, p2);
            if (K.size() < 3) return {};
        }

        // (2) Q = P ∩ K by a boolean intersection, which handles the nonconvex P
        //     and a disconnected result. For a valid winner K ⊆ P and Q = K;
        //     intersecting with P anyway means the result never leaves P.
        if (!K.is_simple()) return {};
        if (K.orientation() != CGAL::COUNTERCLOCKWISE) K.reverse_orientation();
        Polygon subject = poly;
        if (subject.orientation() != CGAL::COUNTERCLOCKWISE) subject.reverse_orientation();

        std::vector<Polygon_with_holes> pieces;
        CGAL::intersection(subject, K, std::back_inserter(pieces));

        // The largest piece; a valid skull is a single convex region.
        const Polygon *best = nullptr;
        FT best_area = FT(-1);
        for (const auto &pwh : pieces) {
            const Polygon &ob = pwh.outer_boundary();
            const FT ar = CGAL::abs(ob.area());
            if (ar > best_area) {
                best_area = ar;
                best = &ob;
            }
        }
        return best ? *best : Polygon{};
    }

    PotatoResult compute_potato_area(const MMatrix &M, const AMatrix &A,
                                     const std::vector<Chord> &Xi, const Polygon &poly) {
        const std::size_t N = Xi.size();

        PotatoResult best;
        best.area = FT(0);  // overwritten on first match (see `found`)
        best.c_idx = 0;
        best.c_prime_idx = 0;

        // Step III: max{ M(C,C') + M(C',C) : C∈Ξ_i, C'∈Ξ_j }. The junctions are
        // seams on the bounding chords themselves (on C from rev.x_last to
        // fwd.x_first, on C' from fwd.x_last to rev.x_first), so the Step II
        // seam and wedge guards apply here too.
        bool found = false;
        for (std::size_t a = 0; a < N; a++) {
            for (std::size_t b = 0; b < N; b++) {
                if (!M.is_valid(a, b)) continue;
                if (!M.is_valid(b, a)) continue;

                const MEntry fwd = M.get(a, b);
                const MEntry rev = M.get(b, a);

                if (!seam_ok(rev, Xi[a], fwd, poly)) continue;
                if (!seam_ok(fwd, Xi[b], rev, poly)) continue;
                if (!in_chain_wedge(fwd.x_first, Xi[a], Xi[b], poly)
                    || !in_chain_wedge(fwd.x_last, Xi[a], Xi[b], poly)
                    || !in_chain_wedge(rev.x_first, Xi[a], Xi[b], poly)
                    || !in_chain_wedge(rev.x_last, Xi[a], Xi[b], poly))
                    continue;

                const FT total = fwd.area + rev.area;

                if (!found || total > best.area) {
                    best.area = total;
                    best.c_idx = a;
                    best.c_prime_idx = b;
                    found = true;
                }
            }
        }

        return best;
    }
#pragma endregion

    Polygon pre_checks(const Polygon &input) {
        Polygon poly = input;

        if (!poly.is_simple()) {
            throw std::invalid_argument("Input polygon must be simple.");
        }

        if (poly.is_clockwise_oriented()) {
            poly.reverse_orientation();
        }

        remove_collinear_vertices(poly);

        if (!poly.is_simple()) {
            throw std::invalid_argument("Input polygon must remain simple after cleanup.");
        }

        return poly;
    }

    Solution convex_skull(const Polygon &polygon) {
        Polygon poly = pre_checks(polygon);

        if (poly.is_convex()) return {poly, std::nullopt, std::nullopt, std::nullopt};

        const auto &reflex_flags = reflex_vertices_flags(poly);
        memprobe::begin(poly.size());
        memprobe::mark("pre-checks");

        // Steps 1-2: Ξ_i = {C_{i,j} | C_{i,j} is defined}, the extremal chords
        // from corner v_i; Xi is their union.
        const std::vector<Chord> Xi = generate_extremal_chords(poly, reflex_flags);

        const std::size_t N = Xi.size();

        // Memory snapshots at phase boundaries, for solve_memory.py: stored
        // entries versus bucket counts of A and M.
        AMatrix *probe_A = nullptr;
        MMatrix *probe_M = nullptr;
        const auto mem_mark = [&](const char *phase) {
            if (!memprobe::enabled()) return;
            memprobe::Snapshot s;
            s.xi = N;
            if (probe_A) {
                s.a_stored = probe_A->data.size();
                s.a_buckets = probe_A->data.bucket_count();
            }
            if (probe_M) {
                s.m_stored = probe_M->data.size();
                s.m_buckets = probe_M->data.bucket_count();
            }
            memprobe::mark(phase, s);
        };
        mem_mark("chords-Xi");

        // Step 3: A(C,C') is the area of the unique balanced (C,C')-chain, if it is
        // admissible: (C_{i,j}, C_{k,l}) is admissible, and the first and last
        // nodes x_0 and x_m lie in C_{i,j} and C_{k,l}.
        AMatrix A;
        // Reserved from observed occupancy (A at most about 0.05*|Xi|^2 entries,
        // M about 0.16*|Xi|^2) with headroom, rather than the N^2 key space.
        A.reserve(N * N / 8 + 16);
        probe_A = &A;
        mem_mark("A-reserved");

        // "It is easy by a brute force method to compute all the (C_0, C_1)-chains
        // of lengths 0 or 1 in time O(n^6)."
        for (std::size_t a = 0; a < N; a++) {
            for (std::size_t b = 0; b < N; b++) {
                // Step 4: skip pairs that are not admissible.
                if (!is_admissible_pair(Xi[a], Xi[b], poly.size())) continue;

                // R(C_{i,j}, C_{k,l}): the corners of P from v_j to v_k inclusive.
                auto R = compute_R(Xi[a], Xi[b], poly.size());

                // Length-0 chain
                if (auto entry = try_length_0_chain(R, Xi, a, b, poly, reflex_flags);
                    entry.has_value()) {
                    A.set(a, b, entry.value());
                    continue;
                }

                // Length-1 chain. The two cases are mutually exclusive: length 0
                // needs R ∩ C+ ∩ C'+ empty, length 1 a pivot strictly inside it.
                if (auto chain_entry = try_length_1_chain(R, Xi, a, b, poly, reflex_flags);
                    chain_entry.has_value()) {
                    A.set(a, b, chain_entry.value());
                }
                // Otherwise the entry stays absent (-inf).
            }
        }

        mem_mark("A-atoms");

        // Step 6
        compute_A_matrix_all_stages(A, Xi, poly, reflex_flags, poly.size());
        mem_mark("A-composed");

        // M(C,C'): the maximum area of admissible (C,C')-chains (Section 6); an
        // absent key is -inf.
        MMatrix M;
        M.reserve(N * N / 4 + 16);
        probe_M = &M;
        mem_mark("M-reserved");

        // Ξ_i: the chords grouped by start vertex i.
        std::vector<std::vector<std::size_t>> Xi_by_vertex(poly.size());

        for (std::size_t idx = 0; idx < Xi.size(); idx++) Xi_by_vertex[Xi[idx].i].push_back(idx);

        // Stage s = 1, the base case (Section 6).
        compute_M_stage_1(M, Xi, Xi_by_vertex, poly, reflex_flags);
        mem_mark("M-stage-1");

        // Stages s = 2..n-1 (Section 6).
        for (std::size_t s = 2; s < poly.size(); s++) {
            compute_M_stage_s(M, A, Xi, Xi_by_vertex, poly, s);
        }
        mem_mark("M-all-stages");

        const auto &result = compute_potato_area(M, A, Xi, poly);
        mem_mark("step-III");

        auto potato_polygon
            = reconstruct_potato(M, A, Xi, poly, reflex_flags, result.c_idx, result.c_prime_idx);

        return {potato_polygon, N, A.data.size(), M.data.size()};
    }
}  // namespace potato_peeling::chang_and_yap

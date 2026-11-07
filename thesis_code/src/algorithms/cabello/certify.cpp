#include "potato_peeling_core/algorithms/cabello/certify.hpp"

#include <algorithm>
#include <variant>
#include <cstddef>
#include <iterator>

#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>

namespace potato_peeling::cabello {
    namespace {
        // Parameter of q along the ray a -> b, scaled by |b-a|^2 so it stays
        // rational: t(q) = (q - a) . (b - a), with t(a) = 0 and t(b) = |b-a|^2.
        FT param_along(const Point &a, const Vector &d, const Point &q) {
            return (q - a) * d;
        }
    }  // namespace

    bool point_in_polygon(const Polygon &p, const Point &a) {
        // Closed region: a point on the boundary counts as inside.
        return p.bounded_side(a) != CGAL::ON_UNBOUNDED_SIDE;
    }

    bool segment_in_polygon(const Polygon &p, const Point &a, const Point &b) {
        if (a == b) return point_in_polygon(p, a);

        const Segment s(a, b);
        const Vector d = b - a;
        const FT t_end = d * d;  // = |b-a|^2 > 0

        // Cut [a,b] at every point where it meets the boundary, then test the
        // midpoint of each resulting piece. A piece whose midpoint is inside the
        // closed region cannot leave it, because leaving would require another
        // boundary crossing, and every crossing is a cut point.
        std::vector<FT> cuts;
        cuts.reserve(2 * p.size() + 2);
        cuts.push_back(FT(0));
        cuts.push_back(t_end);

        for (auto e = p.edges_begin(); e != p.edges_end(); ++e) {
            const auto result = CGAL::intersection(s, *e);
            if (!result) continue;
            if (const Point *q = std::get_if<Point>(&*result)) {
                cuts.push_back(param_along(a, d, *q));
            } else if (const Segment *sub = std::get_if<Segment>(&*result)) {
                // The chord runs along this edge: both ends are cut points.
                cuts.push_back(param_along(a, d, sub->source()));
                cuts.push_back(param_along(a, d, sub->target()));
            }
        }

        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

        for (std::size_t i = 0; i + 1 < cuts.size(); ++i) {
            const FT tm = (cuts[i] + cuts[i + 1]) / FT(2);
            const Point mid = a + d * (tm / t_end);
            if (!point_in_polygon(p, mid)) return false;
        }
        return true;
    }

    Certificate certify(const Polygon &p, const std::vector<Point> &candidate) {
        Certificate cert;

        std::vector<Point> pts = candidate;
        std::sort(pts.begin(), pts.end());
        pts.erase(std::unique(pts.begin(), pts.end()), pts.end());

        if (pts.empty()) {
            cert.convex = true;
            cert.contained = true;
            cert.message = "empty candidate";
            return cert;
        }

        // A point or a segment is a legal zero-area answer.
        if (pts.size() == 1) {
            cert.convex = true;
            cert.contained = point_in_polygon(p, pts[0]);
            cert.polygon.push_back(pts[0]);
            if (!cert.contained) cert.message = "single point outside P";
            return cert;
        }

        std::vector<Point> hull;
        CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull));

        if (hull.size() < 3) {
            cert.convex = true;
            cert.contained = segment_in_polygon(p, hull.front(), hull.back());
            for (const auto &q : hull) cert.polygon.push_back(q);
            if (!cert.contained) cert.message = "degenerate segment leaves P";
            return cert;
        }

        for (const auto &q : hull) cert.polygon.push_back(q);

        // Redundant after convex_hull_2, but checked independently.
        if (!cert.polygon.is_simple()) {
            cert.message = "hull is not simple";
            return cert;
        }
        if (!cert.polygon.is_convex()) {
            cert.message = "hull is not convex";
            return cert;
        }
        cert.convex = true;

        // Convex and every edge inside P, so the whole polygon is (certify.hpp).
        for (auto e = cert.polygon.edges_begin(); e != cert.polygon.edges_end(); ++e) {
            if (!segment_in_polygon(p, e->source(), e->target())) {
                cert.message = "hull edge leaves P";
                return cert;
            }
        }
        cert.contained = true;

        if (cert.polygon.orientation() == CGAL::CLOCKWISE) cert.polygon.reverse_orientation();
        cert.area = CGAL::abs(cert.polygon.area());
        CGAL::exact(cert.area);
        return cert;
    }

}  // namespace potato_peeling::cabello

#include "potato_peeling_core/geometry.hpp"

#include <cstddef>

namespace potato_peeling {
    std::vector<bool> reflex_vertices_flags(const Polygon &polygon) {
        const std::size_t n = polygon.size();

        std::vector<bool> reflex_flags(n, false);

        const bool is_ccw = polygon.is_counterclockwise_oriented();
        const CGAL::Orientation reflex_turn = is_ccw ? CGAL::RIGHT_TURN : CGAL::LEFT_TURN;

        for (std::size_t i = 0; i < n; ++i) {
            const Point &prev = polygon[(i + n - 1) % n];
            const Point &curr = polygon[i];
            const Point &next = polygon[(i + 1) % n];

            if (CGAL::orientation(prev, curr, next) == reflex_turn) {
                reflex_flags[i] = true;
            }
        }

        return reflex_flags;
    }

    void remove_collinear_vertices(Polygon &polygon) {
        bool changed = true;

        while (changed) {
            changed = false;

            const std::size_t n = polygon.size();
            if (n <= 3) {
                return;
            }

            Polygon cleaned;

            for (std::size_t i = 0; i < n; ++i) {
                const Point &prev = polygon[(i + n - 1) % n];
                const Point &curr = polygon[i];
                const Point &next = polygon[(i + 1) % n];

                const bool removable = CGAL::collinear(prev, curr, next)
                                       && CGAL::collinear_are_ordered_along_line(prev, curr, next);

                if (!removable) {
                    cleaned.push_back(curr);
                } else {
                    changed = true;
                }
            }

            polygon = cleaned;
        }
    }

    Polygon counterclockwise(const Polygon &p) {
        Polygon q = p;
        if (q.is_clockwise_oriented()) q.reverse_orientation();
        return q;
    }

    FT path_area(const std::vector<Point> &path) {
        if (path.size() < 3) return FT(0);
        FT s{0};
        for (std::size_t i = 0; i < path.size(); ++i) {
            const Point &a = path[i];
            const Point &b = path[(i + 1) % path.size()];
            s += a.x() * b.y() - b.x() * a.y();
        }
        return CGAL::abs(s) / FT(2);
    }
}  // namespace potato_peeling

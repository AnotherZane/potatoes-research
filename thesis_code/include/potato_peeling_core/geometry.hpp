#pragma once

#include <vector>
#include "potato_peeling_core/types.hpp"

// Polygon helpers shared by the algorithms.

namespace potato_peeling {
    inline constexpr double PI = 3.14159265358979323846;

    // Orders points by their projection onto a direction.
    struct SortAlongLine {
        Direction direction;

        explicit SortAlongLine(const Direction &dir) : direction(dir) {}

        bool operator()(const Point &p1, const Point &p2) const {
            const auto t1 = direction.dx() * p1.x() + direction.dy() * p1.y();
            const auto t2 = direction.dx() * p2.x() + direction.dy() * p2.y();
            return t1 < t2;
        }
    };

    // One flag per vertex: true at reflex corners, for either orientation.
    std::vector<bool> reflex_vertices_flags(const Polygon &polygon);

    // Removes every vertex lying on the segment between its neighbours.
    void remove_collinear_vertices(Polygon &polygon);

    // Copy of p oriented counter-clockwise.
    Polygon counterclockwise(const Polygon &p);

    // |shoelace| / 2 of a closed path.
    FT path_area(const std::vector<Point> &path);
}  // namespace potato_peeling

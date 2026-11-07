#pragma once

// Exact certification of a candidate answer in EPECK. The candidate is a plain
// list of points and is checked from scratch.
//
// For a simple polygon P (no holes), a convex K whose boundary lies in P lies in
// P: a point of int(K) outside P could reach infinity without crossing P, hence
// without crossing boundary(K), which is impossible from inside a convex set. So
// "convex, and every hull edge lies in P" proves containment. The same argument
// is why Cabello's DP only needs to check consecutive boundary pairs.

#include <string>
#include <vector>

#include "potato_peeling_core/types.hpp"

namespace potato_peeling::cabello {

    struct Certificate {
        bool convex{false};
        // Convex and every edge inside P.
        bool contained{false};
        // Exact rational area of the certified polygon (0 for degenerate input).
        FT area{0};
        // CCW, collinear vertices removed. Empty when the input had < 3 distinct
        // points or was not certifiable.
        Polygon polygon;
        // Non-empty iff a check failed; names the first failure.
        std::string message;

        bool ok() const { return convex && contained; }
    };

    // `candidate` is treated as a point set and its convex hull is certified,
    // matching the paper's line 15, `return conv(best)`.
    Certificate certify(const Polygon &p, const std::vector<Point> &candidate);

    // True iff every point of [a,b] lies in the closed region bounded by the
    // simple polygon p. Touching or running along the boundary is allowed.
    bool segment_in_polygon(const Polygon &p, const Point &a, const Point &b);

    // Exact point-in-polygon (closed region: boundary counts as inside).
    bool point_in_polygon(const Polygon &p, const Point &a);

}  // namespace potato_peeling::cabello

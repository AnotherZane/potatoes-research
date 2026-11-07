#pragma once

// Line 1 of LargePotato [C17 p15] asks for A(P) with A(P) <= A*(P) <= C2 * A(P).
// The paper cites Hall-Holt et al. for this; we do not implement it. The lower
// half only needs a convex subset of P, which both functions below return. The
// upper half is replaced by U = area(P) >= A*(P) in large_potato.

#include "potato_peeling_core/types.hpp"

namespace potato_peeling::cabello {

    // A(P), in the same units as area(P): the polygon is never rescaled, so
    // wherever the paper assumes area(P) = 1 the ratio is carried explicitly.
    double a_bound(const Polygon &p);

    // Largest-area triangle with corners at vertices of P that is contained in P.
    // A triangle whose edges lie in P lies in P, so the test is three exact segment
    // queries over a precomputed visibility matrix: O(n^3) overall.
    double largest_vertex_triangle(const Polygon &p);

}  // namespace potato_peeling::cabello

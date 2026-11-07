#pragma once

#include <cstddef>
#include <optional>
#include <utility>

#include "potato_peeling_core/algorithms/aronov/mesh.hpp"
#include "potato_peeling_core/types.hpp"

namespace potato_peeling::aronov {

    // The predicates of Aronov §2 (p353-p356). All of them are either
    // construction-free (orientation tests) or a single exact rational evaluation.
    // Every one assumes the mesh has been sweep-normalised, so no edge is vertical
    // and no two vertices share an x-coordinate.

    // y of the supporting line of e at abscissa x. Total: after sweep normalisation
    // no edge of the mesh is vertical, so the divisor is never zero.
    FT edge_y_at(const Mesh &mesh, std::size_t e, const FT &x);

    // X(e, e') = X(e) intersect X(e'), as the OPEN interval of x-coordinates.
    // Returns nullopt when the spans are disjoint or meet in a single point: a
    // degenerate trapezoid carries no area and is never a usable pair.
    std::optional<std::pair<FT, FT>> overlap(const Mesh &mesh, std::size_t e, std::size_t e_prime);

    // area(T(e, e')): the trapezoid bounded above by e, below by e', and left and
    // right by the vertical lines through the ends of X(e, e') (p353).
    //
    // Closed form, so no polygon is ever constructed. Precondition: e is above e'.
    FT trapezoid_area(const Mesh &mesh, std::size_t e, std::size_t e_prime);

    // Admissibility conditions 2 and 3 (p354): a collection admissible for (e, e')
    // must contain the triangle BELOW e and the triangle ABOVE e'. If either side is
    // the exterior, no admissible collection exists and Q(e, e') is -infinity.
    //
    // The mirrored halves of those conditions ("but not the triangle above e", "not
    // the one below e'") need no test: the DP builds the region out of vertical
    // strips bounded above by e and below by e', so those triangles are never added.
    bool pair_precondition(const Mesh &mesh, std::size_t e, std::size_t e_prime);

    // Number of leading entries f of E_L(v) for which f followed by e makes no right
    // turn at v, where e is the outgoing bottom-chain edge and v = left(e);
    // equivalently 1 + index of opp(e) (p356), or 0 when opp(e) is undefined. A
    // linear scan: deg(v) is small.
    std::size_t bottom_valid_prefix(const Mesh &mesh, std::size_t e);

    // Index of the first entry f of E_L(v) for which e followed by f makes no right
    // turn at v, where e is the incoming top-chain edge and v = left(e);
    // E_L(v).size() when there is none.
    std::size_t top_valid_suffix(const Mesh &mesh, std::size_t e);
}  // namespace potato_peeling::aronov

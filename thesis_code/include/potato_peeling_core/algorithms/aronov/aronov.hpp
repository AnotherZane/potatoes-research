#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "boost/unordered_map.hpp"
#include "potato_peeling_core/algorithms/aronov/mesh.hpp"
#include "potato_peeling_core/types.hpp"

namespace potato_peeling::aronov {

    // An ordered pair of edges (e, e') with e ABOVE e' — the index of the paper's
    // Q-table (p354).
    using EdgePair = std::pair<std::size_t, std::size_t>;

    struct QEntry {
        FT area;
        // The sub-pair whose Q achieved the maximum; nullopt means the base case,
        // where e and e' share their left endpoint (p356).
        std::optional<EdgePair> from;
    };

    // Sparse: only promising pairs carry a value, and absent means the paper's
    // -infinity (Lazy_exact_nt has no -infinity to store).
    using QTable = boost::unordered_map<EdgePair, QEntry>;

    struct MeshedSkull {
        FT area{0};
        // Triangles of M forming the region.
        std::vector<std::size_t> faces;
        // CCW boundary, in the CALLER's coordinate frame (unsheared).
        Polygon boundary;
        // The Step-III pair that won, i.e. the pair sharing a right endpoint whose
        // Q was maximal. Absent only for an empty result.
        std::optional<EdgePair> winner;
    };

    // Step III (p354): the answer is the maximum of Q over the promising pairs that
    // share their RIGHT endpoint — only then does the clipping line x = max X(e,e')
    // pass through that endpoint, making U(C) the whole union rather than a clipped
    // one.
    std::optional<EdgePair> best_pair(const Mesh &mesh, const QTable &q);

    // Walk the backpointers of `winner` and materialise the region: its two convex
    // chains, the boundary polygon, and the triangles enclosed by them.
    MeshedSkull reconstruct(const Mesh &mesh, const QTable &q, EdgePair winner);

    // Checks that the DP value equals the reconstructed area exactly, the boundary is
    // simple and convex, and the region is made of mesh triangles.
    struct SkullCheck {
        bool ok{true};
        std::string message;

        explicit operator bool() const { return ok; }
    };

    SkullCheck check_skull(const Mesh &mesh, const MeshedSkull &skull);

    // Largest-area convex subpolygon of M that is a union of triangles (Theorem 1).
    MeshedSkull largest_convex_subpolygon(const Mesh &mesh);

    // Validate and clean a polygon: simple, collinear vertices removed, still simple.
    Polygon pre_checks(const Polygon &input);

    // The largest convex union of triangles of P's default mesh. A convex P is
    // returned as is.
    Polygon convex_skull(const Polygon &polygon);
}  // namespace potato_peeling::aronov

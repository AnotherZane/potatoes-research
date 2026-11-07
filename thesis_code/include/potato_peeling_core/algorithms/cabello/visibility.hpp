#pragma once

// The visibility oracle and the visibility-graph builder.
//
// The paper builds G(P,R) with Ben-Moshe et al. (Theorem 13 [C17 p15]) and
// answers segment queries with Chazelle et al. [C17 p16]. We use a CGAL AABB tree
// over the boundary edges instead: O(n log n) to build, O(log n + k) per query.
//
// `visible(a,b)` requires both endpoints strictly inside P and rejects any
// segment that meets the boundary. It constructs nothing, so it is exact. It is
// conservative: a segment that grazes a reflex vertex while staying inside is
// rejected. That has probability zero for sampled points and can only shrink the
// answer, never let it leave P. certify.hpp has the exact test that allows
// boundary contact.

#include <memory>
#include <vector>

#include "potato_peeling_core/algorithms/cabello/graph.hpp"
#include "potato_peeling_core/algorithms/cabello/kernel.hpp"

namespace potato_peeling::cabello {

    // Immutable after construction; shared across large_potato's worker threads.
    class VisibilityOracle {
    public:
        explicit VisibilityOracle(const Poly &p);

        ~VisibilityOracle();

        VisibilityOracle(VisibilityOracle &&) noexcept;

        VisibilityOracle &operator=(VisibilityOracle &&) noexcept;

        VisibilityOracle(const VisibilityOracle &) = delete;

        VisibilityOracle &operator=(const VisibilityOracle &) = delete;

        // Strictly inside: the boundary does not count.
        bool inside(const Pt &q) const;

        // See the header comment. False if either endpoint is not strictly inside.
        bool visible(const Pt &a, const Pt &b) const;

    private:
        struct Impl;
        Poly poly_;
        std::unique_ptr<Impl> impl_;
    };

    // G(P, points) [C17 p2]. O(|points|^2 log n).
    PointGraph visibility_graph(const VisibilityOracle &oracle, const std::vector<Pt> &points);

}  // namespace potato_peeling::cabello

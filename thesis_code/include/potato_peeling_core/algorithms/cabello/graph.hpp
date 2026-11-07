#pragma once

// The visibility graph of a point set, and its convex cliques.
//
// [C17 p2]: "The visibility graph of R, denoted G(P,R), has R as vertex set and,
// for any two points x and y in R, the edge xy is in G(P,R) whenever the segment
// xy is contained in P."
//
// [C17 p9]: "A set of vertices U from H is a convex clique if: (i) there is an
// edge between any two vertices of U, and (ii) the points of U are in convex
// position."

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "potato_peeling_core/algorithms/cabello/kernel.hpp"

namespace potato_peeling::cabello {

    class PointGraph {
    public:
        PointGraph() = default;

        explicit PointGraph(std::vector<Pt> pts)
            : points_(std::move(pts)), adj_(points_.size() * points_.size(), 0) {}

        std::size_t size() const { return points_.size(); }

        const std::vector<Pt> &points() const { return points_; }

        const Pt &point(const std::size_t i) const { return points_[i]; }

        bool edge(const std::size_t i, const std::size_t j) const {
            return adj_[i * points_.size() + j] != 0;
        }

        void add_edge(const std::size_t i, const std::size_t j) {
            adj_[i * points_.size() + j] = 1;
            adj_[j * points_.size() + i] = 1;
        }

    private:
        std::vector<Pt> points_;
        std::vector<std::uint8_t> adj_;  // row-major, symmetric, zero diagonal
    };

    // The result of phi(H, s): a convex clique, its area, and its vertices in
    // CCW boundary order with `s` last.
    struct Clique {
        double area{0};
        std::vector<std::size_t> vertices;  // indices into the graph

        bool empty() const { return vertices.empty(); }
    };

    // The order that decides which point is "highest" [C17 p10]. Ties in y, which
    // the paper ignores as probability zero, are broken by x.
    inline bool higher(const Pt &a, const Pt &b) {
        if (a.y() != b.y()) return a.y() > b.y();
        return a.x() > b.x();
    }

}  // namespace potato_peeling::cabello

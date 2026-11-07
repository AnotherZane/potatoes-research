#pragma once

// Uniform sampling in a simple polygon, and in the parallelogram Gamma.
//
// [C17 p14]: "By the linear-time algorithm of Chazelle, we compute a
// triangulation of P, giving triangles T_1..T_{n-2}. We then compute the prefix
// sums S_i = area(T_1) + ... + area(T_i). To sample a point, we select a random
// number x in [0,1], perform a binary search to find the smallest index j such
// that x <= S_j, and sample a random point inside T_j. A random point inside T_j
// can be generated using a random point inside a parallelogram that contains two
// congruent copies of T_j."
//
// We use CGAL's constrained Delaunay triangulation instead of Chazelle's, as the
// paper's footnote 2 [C17 p14] suggests for practical use.

#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <random>
#include <vector>

#include "potato_peeling_core/algorithms/cabello/kernel.hpp"

namespace potato_peeling::cabello {

    struct Triangulation {
        std::vector<std::array<Pt, 3>> triangles;
        // prefix[k] = sum of areas of triangles[0..k]. prefix.back() == total_area.
        std::vector<double> prefix;
        double total_area{0};

        std::size_t size() const { return triangles.size(); }

        bool empty() const { return triangles.empty(); }
    };

    // Triangulate a simple polygon using only its own vertices.
    Triangulation triangulate(const Poly &p);

    // Uniform in a triangle, by the paper's parallelogram-reflection trick.
    Pt sample_in_triangle(const std::array<Pt, 3> &t, std::mt19937_64 &gen);

    // Uniform in P: pick a triangle with probability proportional to its area,
    // then a point inside it.
    Pt sample_uniform(const Triangulation &tri, std::mt19937_64 &gen);

    // `count` uniform points accepted by `inside`. A point that rounding puts
    // outside P is re-drawn from the same stream; `retries` counts the re-draws.
    std::vector<Pt> sample_uniform(const Triangulation &tri, std::size_t count,
                                   std::mt19937_64 &gen,
                                   const std::function<bool(const Pt &)> &inside,
                                   std::size_t *retries = nullptr);

    // The parallelogram Gamma(a, b, A) of [C17 p8], as origin + u*e1 + v*e2 with
    // u, v in [0,1].
    //
    //   a' = 2a - b,  b' = 2b - a,  d = 2A / |y(a) - y(b)|
    //
    // so |a'b'| = 3|ab| and area = 12A. (The 2014 conference version states
    // |a'b'| = 2|ab|, which is wrong.)
    struct Parallelogram {
        Pt origin;
        Vec e1, e2;
        double area{0};
    };

    // nullopt iff y(a) == y(b), where Gamma is undefined. This has probability 0
    // for sampled points; the caller skips the edge.
    std::optional<Parallelogram> gamma(const Pt &a, const Pt &b, double A);

    Pt sample_uniform(const Parallelogram &g, std::mt19937_64 &gen);

}  // namespace potato_peeling::cabello

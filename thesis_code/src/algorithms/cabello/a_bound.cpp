#include "potato_peeling_core/algorithms/cabello/a_bound.hpp"

#include <cstddef>
#include <vector>

#include "potato_peeling_core/algorithms/aronov/aronov.hpp"
#include "potato_peeling_core/algorithms/cabello/certify.hpp"

namespace potato_peeling::cabello {

    double largest_vertex_triangle(const Polygon &p) {
        const std::size_t n = p.size();
        if (n < 3) return 0.0;

        // Pairwise visibility among the vertices, exact. O(n^2) queries.
        std::vector<char> vis(n * n, 0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j)
                if (segment_in_polygon(p, p[i], p[j])) {
                    vis[i * n + j] = 1;
                    vis[j * n + i] = 1;
                }

        FT best(0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j) {
                if (!vis[i * n + j]) continue;
                for (std::size_t k = j + 1; k < n; ++k) {
                    if (!vis[i * n + k] || !vis[j * n + k]) continue;
                    // Every edge lies in P, so the triangle does (certify.hpp).
                    if (FT a = CGAL::abs(CGAL::area(p[i], p[j], p[k])); a > best) {
                        CGAL::exact(a);
                        best = a;
                    }
                }
            }
        return CGAL::to_double(best);
    }

    double a_bound(const Polygon &p) {
        // Aronov's convex skull is a convex subset of P, so a lower bound on A*(P).
        const Polygon skull = aronov::convex_skull(p);
        if (skull.size() < 3) return largest_vertex_triangle(p);
        FT a = CGAL::abs(skull.area());
        CGAL::exact(a);
        return CGAL::to_double(a);
    }

}  // namespace potato_peeling::cabello

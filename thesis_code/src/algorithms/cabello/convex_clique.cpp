#include "potato_peeling_core/algorithms/cabello/convex_clique.hpp"
#include "potato_peeling_core/algorithms/cabello/kernel.hpp"

#include <algorithm>
#include <vector>

namespace potato_peeling::cabello {
    namespace {
        constexpr std::size_t NONE = static_cast<std::size_t>(-1);

        // Angular half of a direction: [0, pi) -> 0, [pi, 2pi) -> 1. Pure
        // predicate on the sign bits, so it is exact.
        int half(const Vec &v) {
            if (v.y() < 0) return 1;
            if (v.y() > 0) return 0;
            return v.x() < 0 ? 1 : 0;
        }

        // Strict angular order on directions, measured from the +x axis. Exact:
        // `half` is a sign test, and within a half the order is an exact
        // orientation. The three-point orientation is used because CGAL filters it
        // faster than the two-vector form, and this is phi's innermost operation.
        bool ang_less(const Vec &a, const Vec &b) {
            const int ha = half(a), hb = half(b);
            if (ha != hb) return ha < hb;
            return CGAL::orientation(Pt(CGAL::ORIGIN), CGAL::ORIGIN + a, CGAL::ORIGIN + b)
                   == CGAL::LEFT_TURN;
        }

        // The candidates for phi(H, s), sorted CCW radially from s.
        std::vector<std::size_t> radial_candidates(const PointGraph &h, const std::size_t s) {
            std::vector<std::size_t> cand;
            for (std::size_t k = 0; k < h.size(); ++k) {
                if (k == s || !h.edge(k, s)) continue;
                if (higher(h.point(k), h.point(s))) continue;
                cand.push_back(k);
            }
            const Pt &sp = h.point(s);
            // Must be exact: the DP assumes position in this array increases with
            // angle around s. An inexact cross product can call nearly collinear
            // candidates collinear, order them by distance instead, and let phi
            // return a clique that winds twice while every turn is still left.
            std::sort(cand.begin(), cand.end(), [&](const std::size_t a, const std::size_t b) {
                const CGAL::Orientation o = CGAL::orientation(sp, h.point(a), h.point(b));
                if (o != CGAL::COLLINEAR) return o == CGAL::LEFT_TURN;
                return CGAL::compare_distance_to_point(sp, h.point(a), h.point(b)) == CGAL::SMALLER;
            });
            return cand;
        }

        // Range argmax over a fixed array, O(m log m) build, O(1) query.
        class SparseTable {
        public:
            explicit SparseTable(const std::vector<double> &v) : v_(v) {
                const std::size_t m = v.size();
                if (m == 0) return;
                std::size_t levels = 1;
                while ((1u << levels) <= m) ++levels;
                table_.assign(levels, std::vector<std::size_t>(m));
                for (std::size_t i = 0; i < m; ++i) table_[0][i] = i;
                for (std::size_t k = 1; k < levels; ++k) {
                    const std::size_t len = 1u << k;
                    for (std::size_t i = 0; i + len <= m; ++i) {
                        const std::size_t a = table_[k - 1][i];
                        const std::size_t b = table_[k - 1][i + len / 2];
                        table_[k][i] = v_[a] >= v_[b] ? a : b;
                    }
                }
            }

            // argmax over [lo, hi); NONE if empty.
            std::size_t query(const std::size_t lo, const std::size_t hi) const {
                if (lo >= hi) return NONE;
                std::size_t k = 0;
                while ((1u << (k + 1)) <= hi - lo) ++k;
                const std::size_t a = table_[k][lo];
                const std::size_t b = table_[k][hi - (1u << k)];
                return v_[a] >= v_[b] ? a : b;
            }

        private:
            const std::vector<double> &v_;
            std::vector<std::vector<std::size_t>> table_;
        };
    }  // namespace

    Clique phi(const PointGraph &h, const std::size_t s) {
        Clique best;
        double best_val = 0.0;
        const std::vector<std::size_t> x = radial_candidates(h, s);
        const std::size_t n = x.size();
        if (n < 2) return best;

        const Pt &sp = h.point(s);

        std::vector<double> opt(n * n, 0.0);
        std::vector<char> has(n * n, 0);
        std::vector<std::size_t> parent(n * n, NONE);

        // Per-i scratch: the arrivals h < i with x_h x_i in E(H), sorted by the
        // angle of x_i - x_h, with their Opt[h,i] values.
        std::vector<std::size_t> arrivals;
        std::vector<Vec> dirs;
        std::vector<double> vals;

        for (std::size_t i = 0; i < n; ++i) {
            arrivals.clear();
            for (std::size_t g = 0; g < i; ++g)
                if (has[g * n + i] && h.edge(x[g], x[i])) arrivals.push_back(g);

            std::sort(
                arrivals.begin(), arrivals.end(), [&](const std::size_t a, const std::size_t b) {
                    return ang_less(h.point(x[i]) - h.point(x[a]), h.point(x[i]) - h.point(x[b]));
                });

            const std::size_t m = arrivals.size();
            dirs.clear();
            vals.clear();
            dirs.reserve(m);
            vals.reserve(m);
            for (const std::size_t g : arrivals) {
                dirs.push_back(h.point(x[i]) - h.point(x[g]));
                vals.push_back(opt[g * n + i]);
            }
            const SparseTable range_max(vals);

            for (std::size_t j = i + 1; j < n; ++j) {
                if (!h.edge(x[i], x[j])) continue;

                double best_sub = 0.0;
                std::size_t best_h = NONE;

                if (m > 0) {
                    // Valid arrivals are those with cross(dir, w) > 0, where
                    // w = x_j - x_i: the open arc of directions strictly between
                    // angle(-w) and angle(w). Contiguous in angular order, so two
                    // binary searches bound it (it may wrap past 2pi).
                    const Vec w = h.point(x[j]) - h.point(x[i]);
                    const Vec minus_w = -w;

                    const auto lo = static_cast<std::size_t>(
                        std::upper_bound(dirs.begin(), dirs.end(), minus_w, ang_less)
                        - dirs.begin());
                    const auto hi = static_cast<std::size_t>(
                        std::lower_bound(dirs.begin(), dirs.end(), w, ang_less) - dirs.begin());

                    // The arc wraps past angle 0 exactly when w is in the upper half,
                    // since angle(-w) and angle(w) differ by pi. Comparing lo and hi
                    // instead would miss a wrapped arc that starts at index 0.
                    const bool wraps = half(w) == 0;

                    std::size_t arg = NONE;
                    if (!wraps) {
                        arg = range_max.query(lo, hi);
                    } else {
                        const std::size_t a = range_max.query(lo, m);
                        const std::size_t b = range_max.query(0, hi);
                        if (a == NONE)
                            arg = b;
                        else if (b == NONE)
                            arg = a;
                        else
                            arg = vals[a] >= vals[b] ? a : b;
                    }
                    if (arg != NONE && vals[arg] > 0.0) {
                        best_sub = vals[arg];
                        best_h = arrivals[arg];
                    }
                }

                // Opt[i,j] is the fan area, so a missing predecessor contributes 0
                // and the triangle (s, x_i, x_j) is the base.
                const double area
                    = CGAL::to_double(CGAL::abs(CGAL::area(sp, h.point(x[i]), h.point(x[j]))));
                opt[i * n + j] = area + best_sub;
                has[i * n + j] = 1;
                parent[i * n + j] = best_h;

                const double scored = opt[i * n + j];
                if (scored > best_val) {
                    best_val = scored;
                    best.vertices = {i, j};
                }
            }
        }

        best.area = best_val;
        if (best.vertices.empty()) return best;

        std::vector<std::size_t> chain;
        std::size_t i = best.vertices[0], j = best.vertices[1];
        chain.push_back(j);
        chain.push_back(i);
        while (parent[i * n + j] != NONE) {
            const std::size_t g = parent[i * n + j];
            chain.push_back(g);
            j = i;
            i = g;
        }
        std::reverse(chain.begin(), chain.end());

        best.vertices.clear();
        best.vertices.reserve(chain.size() + 1);
        for (const std::size_t k : chain) best.vertices.push_back(x[k]);
        best.vertices.push_back(s);
        return best;
    }

}  // namespace potato_peeling::cabello

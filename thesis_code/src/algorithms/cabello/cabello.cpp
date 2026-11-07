#include "potato_peeling_core/algorithms/cabello/cabello.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <atomic>
#include <thread>

#include <utility>

#include "potato_peeling_core/algorithms/cabello/rng.hpp"
#include "potato_peeling_core/algorithms/cabello/graph.hpp"
#include "potato_peeling_core/algorithms/cabello/a_bound.hpp"
#include "potato_peeling_core/algorithms/cabello/certify.hpp"
#include "potato_peeling_core/algorithms/cabello/convex_clique.hpp"
#include "potato_peeling_core/algorithms/cabello/sampling.hpp"
#include "potato_peeling_core/algorithms/cabello/visibility.hpp"
#include "potato_peeling_core/geometry.hpp"

namespace potato_peeling::cabello {
    namespace {
        std::size_t sized(const double literal, const double scale, const std::size_t floor_value) {
            const double v = std::ceil(literal * scale);
            if (!std::isfinite(v) || v < static_cast<double>(floor_value)) return floor_value;
            // Cap sizes that could not be allocated.
            constexpr double cap = 5e7;
            return static_cast<std::size_t>(std::min(v, cap));
        }

        struct EdgeResult {
            double area{0};
            std::vector<Pt> points;
            std::size_t phi_calls{0};
            std::size_t r_in{0}, s_in{0};
            bool skipped_horizontal{false};
        };

        // Lines 9-14 for one edge ab of G(P,R). Draws only from substreams
        // addressed by (iteration, edge), so the result does not depend on which
        // thread runs it or when.
        EdgeResult search_edge(const Options &opt, const VisibilityOracle &oracle, const Pt &a,
                               const Pt &b, const double gamma_area, const std::size_t r_ab,
                               const std::size_t s_ab, const std::uint64_t iteration,
                               const std::uint64_t edge_index) {
            EdgeResult out;

            const auto box = gamma(a, b, gamma_area);
            if (!box) {
                // y(a) == y(b) leaves Gamma undefined. Probability zero.
                out.skipped_horizontal = true;
                return out;
            }

            // Lines 9-11: sample in Gamma, keep what lands in P.
            auto gen_r = substream(opt.seed, Stream::SampleRab, iteration, edge_index);
            auto gen_s = substream(opt.seed, Stream::SampleSab, iteration, edge_index);

            std::vector<Pt> pts;
            pts.reserve(r_ab + s_ab);
            for (std::size_t k = 0; k < r_ab; ++k) {
                if (Pt z = sample_uniform(*box, gen_r); oracle.inside(z)) pts.push_back(z);
            }
            out.r_in = pts.size();

            const std::size_t n_r = pts.size();  // R_ab occupies pts[0, n_r)

            for (std::size_t k = 0; k < s_ab; ++k) {
                if (Pt z = sample_uniform(*box, gen_s); oracle.inside(z)) pts.push_back(z);
            }
            out.s_in = pts.size() - n_r;
            if (pts.size() < 3) return out;

            const PointGraph g = visibility_graph(oracle, pts);

            // Lines 12-14: the highest vertex ranges over S_ab only [C17 p16].
            for (std::size_t s = n_r; s < pts.size(); ++s) {
                const Clique u = phi(g, s);
                ++out.phi_calls;
                if (u.area > out.area) {
                    out.area = u.area;
                    out.points.clear();
                    out.points.reserve(u.vertices.size());
                    for (const std::size_t v : u.vertices) out.points.push_back(g.point(v));
                }
            }
            return out;
        }
    }  // namespace

    Polygon pre_checks(const Polygon &input) {
        Polygon poly = input;

        if (!poly.is_simple()) {
            throw std::invalid_argument("Input polygon must be simple.");
        }

        remove_collinear_vertices(poly);

        if (!poly.is_simple()) {
            throw std::invalid_argument("Input polygon must remain simple after cleanup.");
        }

        return poly;
    }

    Result large_potato(const Polygon &p, const Options &opt) {
        const auto t0 = std::chrono::steady_clock::now();

        Result res;
        res.seed = opt.seed;

        if (!(opt.epsilon > 0.0 && opt.epsilon < 1.0))
            throw std::invalid_argument("cabello: epsilon must lie in (0,1)");
        if (!(opt.delta > 0.0 && opt.delta < 1.0))
            throw std::invalid_argument("cabello: delta must lie in (0,1)");
        if (opt.sample_scale <= 0.0 || opt.s_scale <= 0.0)
            throw std::invalid_argument("cabello: sample scales must be positive");

        const Polygon cleaned = pre_checks(p);
        if (cleaned.size() < 3) return res;

        // A convex P is its own answer.
        if (cleaned.is_convex()) {
            res.polygon = cleaned;
            res.area = CGAL::abs(cleaned.area());
            CGAL::exact(res.area);
            res.convex = true;
            res.contained = true;
            res.a_of_p = CGAL::to_double(res.area);
            res.upper_bound = res.a_of_p;
            res.upper_bound_certified = true;
            res.c2_effective = res.a_of_p > 0.0 ? res.upper_bound / res.a_of_p : 0.0;
            res.seconds
                = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            return res;
        }

        const Poly q = to_inexact(cleaned);
        const double area_p = std::abs(CGAL::to_double(q.area()));
        if (area_p <= 0.0) return res;

        const VisibilityOracle oracle(q);
        const Triangulation tri = triangulate(q);
        if (tri.empty()) return res;

        // Line 1.
        const double A = a_bound(cleaned);
        if (A <= 0.0) return res;
        // U = area(P) replaces the paper's C2*A(P). Every use of C2*A(P) only needs
        // an upper bound on A*(P), and area(P) is one, so C2 becomes U/A(P).
        const double U = area_p;
        const double C2 = U / A;
        const double C3 = c3(C2);
        res.a_of_p = A;
        res.upper_bound = U;
        res.upper_bound_certified = true;
        res.c2_effective = C2;

        // The paper's sizes. P is not rescaled to unit area, so area(P) appears
        // explicitly where the paper has 1.
        res.literal_r = 60.0 * area_p / A;                                     // line 2
        res.literal_r_ab = 96.0 * C1 * C2 / std::pow(opt.epsilon / 2.0, 1.5);  // line 9
        res.literal_s_ab = 288.0 * C2 / opt.epsilon;                           // line 10

        const std::size_t r = sized(res.literal_r, 1.0, 3);
        const std::size_t r_ab = sized(res.literal_r_ab, opt.sample_scale, 3);
        const std::size_t s_ab = sized(res.literal_s_ab, opt.s_scale, 1);
        res.r = r;
        res.r_ab = r_ab;
        res.s_ab = s_ab;

        // Line 4: repeat 3*log2(1/delta) times.
        const auto iterations
            = static_cast<std::size_t>(std::ceil(3.0 * std::log2(1.0 / opt.delta)));

        const double edge_budget = C3 * static_cast<double>(cleaned.size());

        double best_area = 0;
        std::vector<Pt> best_points;

        for (std::size_t it = 0; it < iterations; ++it) {
            ++res.iterations;

            // Line 5.
            auto gen_r = substream(opt.seed, Stream::SampleR, it);
            std::size_t retries = 0;
            const std::vector<Pt> R = sample_uniform(
                tri, r, gen_r, [&](const Pt &z) { return oracle.inside(z); }, &retries);

            // Lines 6-7. With C3 around 10^7 the edge budget is never hit in
            // practice, but it bounds the worst case.
            const PointGraph g = visibility_graph(oracle, R);

            std::vector<std::pair<std::size_t, std::size_t>> edges;
            for (std::size_t i = 0; i < R.size(); ++i)
                for (std::size_t j = i + 1; j < R.size(); ++j)
                    if (g.edge(i, j)) edges.emplace_back(i, j);
            res.graph_edges += edges.size();

            if (static_cast<double>(edges.size()) > edge_budget) {
                res.edge_budget_exceeded = true;
                continue;  // line 6 fails: this iteration does no work
            }

            // Search a uniformly random subset of at most max_edges edges. The edge
            // list is built in a fixed order, so the subset depends only on the seed.
            if (opt.max_edges && edges.size() > *opt.max_edges) {
                auto gen_e = substream(opt.seed, Stream::EdgeSubsample, it);
                for (std::size_t k = 0; k < *opt.max_edges; ++k) {
                    const std::size_t pick = k + gen_e() % (edges.size() - k);
                    std::swap(edges[k], edges[pick]);
                }
                edges.resize(*opt.max_edges);
            }
            res.edges_searched += edges.size();

            // Lines 8-14.
            std::vector<EdgeResult> per_edge(edges.size());
            const unsigned threads = std::max(1u, opt.threads);
            if (threads == 1) {
                for (std::size_t e = 0; e < edges.size(); ++e)
                    per_edge[e] = search_edge(opt, oracle, R[edges[e].first], R[edges[e].second], U,
                                              r_ab, s_ab, it, e);
            } else {
                std::vector<std::thread> pool;
                std::atomic<std::size_t> next{0};
                for (unsigned t = 0; t < threads; ++t) {
                    pool.emplace_back([&] {
                        for (;;) {
                            const std::size_t e = next.fetch_add(1);
                            if (e >= edges.size()) return;
                            per_edge[e] = search_edge(opt, oracle, R[edges[e].first],
                                                      R[edges[e].second], U, r_ab, s_ab, it, e);
                        }
                    });
                }
                for (auto &th : pool) th.join();
            }

            // Reduce in index order; strict > keeps the earliest edge on a tie.
            for (std::size_t e = 0; e < per_edge.size(); ++e) {
                const EdgeResult &er = per_edge[e];
                res.phi_calls += er.phi_calls;
                if (er.skipped_horizontal) ++res.skipped_horizontal;
                if (er.area > best_area) {
                    best_area = er.area;
                    best_points = er.points;
                }
            }

            // Early exit: U >= A*(P), so reaching (1-eps)*U is already good enough.
            // Checked between iterations only, so threads do not affect it.
            if (best_area >= (1.0 - opt.epsilon) * U) break;
        }

        std::vector<Point> exact_points;
        exact_points.reserve(best_points.size());
        for (const auto &z : best_points) exact_points.push_back(to_exact(z));

        res.clique = best_points;

        const Certificate cert = certify(cleaned, exact_points);
        res.polygon = cert.polygon;
        res.area = cert.area;
        res.convex = cert.convex;
        res.contained = cert.contained;
        res.message = cert.message;

        res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return res;
    }

}  // namespace potato_peeling::cabello

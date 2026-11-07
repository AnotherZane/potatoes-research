#include <cstdint>
#include <cstddef>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <sstream>

#include "../casters.hpp"
#include "potato_peeling_core/algorithms/cabello/cabello.hpp"
#include "potato_peeling_core/types.hpp"

namespace nb = nanobind;

namespace {
    using namespace potato_peeling;
    using namespace potato_peeling::cabello;

    Options make_options(const double epsilon, const double delta, const std::uint64_t seed,
                         const double sample_scale, const double s_scale,
                         const nb::object &max_edges, const unsigned threads) {
        Options opt;
        opt.epsilon = epsilon;
        opt.delta = delta;
        opt.seed = seed;
        opt.sample_scale = sample_scale;
        opt.s_scale = s_scale;
        opt.threads = threads;

        if (!max_edges.is_none())
            opt.max_edges = nb::cast<std::size_t>(max_edges);
        else
            opt.max_edges.reset();
        return opt;
    }

    nb::dict result_to_dict(const Result &res) {
        nb::list verts;
        for (const auto &v : res.polygon)
            verts.append(nb::make_tuple(CGAL::to_double(v.x()), CGAL::to_double(v.y())));

        FT area = res.area;
        CGAL::exact(area);
        std::ostringstream os;
        os << area.exact();

        // phi's output before certify takes its hull; see Result::clique.
        nb::list clique;
        for (const auto &v : res.clique)
            clique.append(nb::make_tuple(CGAL::to_double(v.x()), CGAL::to_double(v.y())));

        // Exact vertices as rational strings, so a caller can re-certify the answer.
        auto rat = [](const potato_peeling::FT &v) {
            potato_peeling::FT c = v;
            CGAL::exact(c);
            std::ostringstream o;
            o << c.exact();
            return o.str();
        };
        nb::list verts_exact;
        for (const auto &v : res.polygon)
            verts_exact.append(nb::make_tuple(rat(v.x()), rat(v.y())));

        nb::dict out;
        out["vertices"] = verts;
        out["vertices_exact"] = verts_exact;
        out["clique"] = clique;
        out["area"] = CGAL::to_double(area);
        out["area_exact"] = os.str();
        out["convex"] = res.convex;
        out["contained"] = res.contained;
        out["message"] = res.message;
        out["seed"] = res.seed;
        out["iterations"] = res.iterations;
        out["r"] = res.r;
        out["r_ab"] = res.r_ab;
        out["s_ab"] = res.s_ab;
        out["graph_edges"] = res.graph_edges;
        out["edges_searched"] = res.edges_searched;
        out["phi_calls"] = res.phi_calls;
        out["skipped_horizontal"] = res.skipped_horizontal;
        out["a_of_p"] = res.a_of_p;
        out["upper_bound"] = res.upper_bound;
        out["upper_bound_certified"] = res.upper_bound_certified;
        out["c2_effective"] = res.c2_effective;
        out["seconds"] = res.seconds;
        out["edge_budget_exceeded"] = res.edge_budget_exceeded;
        out["literal_r"] = res.literal_r;
        out["literal_r_ab"] = res.literal_r_ab;
        out["literal_s_ab"] = res.literal_s_ab;
        return out;
    }
}  // namespace

void bind_cabello(nb::module_ &m) {
    m.def(
        "large_potato",
        [](const Polygon &polygon, const double epsilon, const double delta,
           const std::uint64_t seed, const double sample_scale, const double s_scale,
           const nb::object &max_edges, const unsigned threads) {
            const Options opt
                = make_options(epsilon, delta, seed, sample_scale, s_scale, max_edges, threads);
            Result res;
            {
                nb::gil_scoped_release release;
                res = large_potato(polygon, opt);
            }
            return result_to_dict(res);
        },
        "large_potato(polygon, ...) -> dict: Algorithm LargePotato [C17 p15]. Returns the "
        "certified polygon, its exact area, the certificates, the seed and full telemetry.",
        nb::arg("polygon"), nb::arg("epsilon") = 0.1, nb::arg("delta") = 1.0 / 3.0,
        nb::arg("seed") = 0, nb::arg("sample_scale") = 3e-4, nb::arg("s_scale") = 3e-2,
        nb::arg("max_edges").none() = nb::int_(64), nb::arg("threads") = 1u);
}

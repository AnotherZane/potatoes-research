#include <sstream>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "casters.hpp"
#include "potato_peeling_core/refinement/refinement.hpp"
#include "potato_peeling_core/types.hpp"

namespace nb = nanobind;

void bind_refinement(nb::module_ &m) {
    m.def(
        "refine",
        [](const potato_peeling::Polygon &polygon, const potato_peeling::Polygon &start,
           const unsigned max_sweeps, const unsigned rot_steps, const double rot_span,
           const unsigned split_steps, const unsigned adjoin_rounds) {
            potato_peeling::refinement::Options opts;
            opts.max_sweeps = max_sweeps;
            opts.rot_steps = rot_steps;
            opts.rot_span = rot_span;
            opts.split_steps = split_steps;
            opts.adjoin_rounds = adjoin_rounds;
            const auto r = potato_peeling::refinement::refine(polygon, start, opts);
            // The vertices are constructed points, so the doubles can round
            // outside P; `vertices_exact` has the rationals.
            nb::list verts, verts_exact;
            for (auto v = r.polygon.vertices_begin(); v != r.polygon.vertices_end(); ++v) {
                verts.append(nb::make_tuple(CGAL::to_double(v->x()), CGAL::to_double(v->y())));
                std::ostringstream xs, ys;
                xs << v->x().exact();
                ys << v->y().exact();
                verts_exact.append(nb::make_tuple(xs.str(), ys.str()));
            }
            std::ostringstream os;
            os << r.area.exact();
            nb::dict out;
            out["ok"] = r.ok;
            out["area"] = CGAL::to_double(r.area);
            out["area_exact"] = os.str();
            out["start_area"] = CGAL::to_double(r.start_area);
            out["vertices"] = verts;
            out["vertices_exact"] = verts_exact;
            out["sweeps"] = r.sweeps;
            out["evaluations"] = r.evaluations;
            out["halfplanes_start"] = r.halfplanes_start;
            out["halfplanes_end"] = r.halfplanes_end;
            out["splits_accepted"] = r.splits_accepted;
            out["adjoins_accepted"] = r.adjoins_accepted;
            out["message"] = r.message;
            return out;
        },
        "refine(polygon, start, max_sweeps=12, rot_steps=9, rot_span=0.30, "
        "split_steps=0, adjoin_rounds=0) -> dict:\n"
        "local ascent on a convex subset of P. `start` must be convex and inside P.\n"
        "Every state refinement accepted was certified convex-and-contained in exact EPECK,\n"
        "so `area_exact` is a genuine LOWER bound on A*(P).\n"
        "\n"
        "CAUTION: `vertices` are rounded to double. They are CONSTRUCTED points, so the\n"
        "rounding can move a vertex OUTSIDE P and re-certifying them reports\n"
        "contained = false. That is a property of the rounding, not of the answer --\n"
        "use `vertices_exact` (num/den strings) to feed the result back into an exact\n"
        "check, e.g. as a seed for another pass.",
        nb::arg("polygon"), nb::arg("start"), nb::arg("max_sweeps") = 12, nb::arg("rot_steps") = 9,
        nb::arg("rot_span") = 0.30, nb::arg("split_steps") = 0, nb::arg("adjoin_rounds") = 0);
}

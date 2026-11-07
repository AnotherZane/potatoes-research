#include <sstream>
#include <boost/functional/hash.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include "../casters.hpp"
#include <CGAL/Polygon_with_holes_2.h>
#include <vector>
#include "potato_peeling_core/algorithms/aronov/aronov.hpp"
#include "potato_peeling_core/algorithms/aronov/mesh.hpp"
#include "potato_peeling_core/types.hpp"

namespace nb = nanobind;

void bind_aronov(nb::module_ &m) {
    m.def(
        "skull_summary",
        [](const potato_peeling::Polygon &polygon) {
            const potato_peeling::aronov::Mesh mesh = potato_peeling::aronov::build_mesh(polygon);
            const auto mesh_check = validate(mesh);
            const potato_peeling::aronov::MeshedSkull skull
                = potato_peeling::aronov::largest_convex_subpolygon(mesh);
            const auto skull_check = check_skull(mesh, skull);

            nb::dict out;
            out["vertices"] = mesh.num_vertices();
            out["faces"] = mesh.num_faces();
            out["edges"] = mesh.num_edges();
            out["mesh_area"] = CGAL::to_double(mesh.total_area());
            out["area"] = CGAL::to_double(skull.area);
            // Exact predicates, decided here. CGAL::to_double of a Lazy_exact_nt
            // reflects the lazy interval, not just the value, so two FTs that are
            // exactly equal can convert to doubles that differ by an ULP — never
            // re-derive these comparisons on the Python side.
            out["area_positive"] = (skull.area > potato_peeling::FT(0));
            out["area_within_mesh"] = (skull.area <= mesh.total_area());
            out["is_whole_mesh"] = (skull.area == mesh.total_area());
            out["used_faces"] = skull.faces.size();
            out["boundary"] = skull.boundary;
            // Exact area and boundary: "area" and "boundary" above are rounded, and
            // rounded vertices cannot be re-verified against P exactly.
            {
                auto rat = [](const potato_peeling::FT &v) {
                    potato_peeling::FT c = v;
                    CGAL::exact(c);
                    std::ostringstream o;
                    o << c.exact();
                    return o.str();
                };
                potato_peeling::FT a = skull.area;
                CGAL::exact(a);
                std::ostringstream oa;
                oa << a.exact();
                out["area_exact"] = oa.str();
                nb::list b_exact;
                for (const auto &pt : skull.boundary)
                    b_exact.append(nb::make_tuple(rat(pt.x()), rat(pt.y())));
                out["boundary_exact"] = b_exact;
            }
            out["mesh_valid"] = mesh_check.ok;
            out["skull_valid"] = skull_check.ok;
            out["message"] = mesh_check.ok ? skull_check.message : mesh_check.message;
            return out;
        },
        "skull_summary(polygon) -> dict: solve, then run the O1-O3 oracles.", nb::arg("polygon"));
}

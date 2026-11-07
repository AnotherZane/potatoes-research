#include <sstream>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <string>
#include "../casters.hpp"
#include "potato_peeling_core/algorithms/chang_and_yap/chang_and_yap.hpp"
#include "potato_peeling_core/algorithms/chang_and_yap/memprobe.hpp"
// certify.hpp's exact segment_in_polygon, for the containment check below.
#include "potato_peeling_core/algorithms/cabello/certify.hpp"
#include "potato_peeling_core/types.hpp"

namespace nb = nanobind;

void bind_chang_and_yap(nb::module_ &m) {
    m.def("pre_checks", potato_peeling::chang_and_yap::pre_checks,
          "pre_checks(polygon: list[tuple[float, float]]) -> list[tuple[float, float]]",
          nb::arg("polygon"));

    // The skull with its exact area. An area recomputed from the rounded double
    // vertices suffers cancellation; here it is computed from the exact polygon.
    m.def(
        "convex_skull_exact",
        [](const potato_peeling::Polygon &polygon) {
            const auto solution = potato_peeling::chang_and_yap::convex_skull(polygon);
            const potato_peeling::Polygon &skull = solution.skull;
            nb::list verts;
            for (const auto &pt : skull)
                verts.append(nb::make_tuple(CGAL::to_double(pt.x()), CGAL::to_double(pt.y())));

            // The vertices as exact rationals, "num/den" per coordinate. Every
            // construction is rational, so nothing is lost, and a caller can
            // test contact with dP exactly.
            auto rat = [](const potato_peeling::FT &v) {
                potato_peeling::FT c = v;
                CGAL::exact(c);
                std::ostringstream o;
                o << c.exact();
                return o.str();
            };
            nb::list verts_exact;
            for (const auto &pt : skull)
                verts_exact.append(nb::make_tuple(rat(pt.x()), rat(pt.y())));

            // The cleaned input ring the skull was computed against.
            const potato_peeling::Polygon cleaned_ring
                = potato_peeling::chang_and_yap::pre_checks(polygon);
            nb::list poly_exact;
            for (const auto &pt : cleaned_ring)
                poly_exact.append(nb::make_tuple(rat(pt.x()), rat(pt.y())));

            nb::dict out;
            out["vertices"] = verts;
            out["vertices_exact"] = verts_exact;
            out["polygon_exact"] = poly_exact;
            out["n_xi"] = solution.n_xi;
            out["a_entries"] = solution.a_entries;
            out["m_entries"] = solution.m_entries;
            if (skull.size() < 3) {
                out["area"] = 0.0;
                out["area_exact"] = std::string("0");
                out["convex"] = true;  // a point or segment is trivially convex
                out["contained"] = true;
                return out;
            }
            potato_peeling::FT a = CGAL::abs(skull.area());
            CGAL::exact(a);
            std::ostringstream os;
            os << a.exact();
            out["area"] = CGAL::to_double(a);
            out["area_exact"] = os.str();

            // Both predicates are decided on the exact polygon, before conversion
            // to double, where a short edge's turn can change sign.
            out["convex"] = skull.is_convex();
            // A convex K whose edges all lie in a simple P lies in P
            // (certify.hpp). Boundary contact is allowed: the optimum's chords
            // graze reflex corners.
            const potato_peeling::Polygon &cleaned = cleaned_ring;
            bool contained = true;
            for (auto e = skull.edges_begin(); e != skull.edges_end(); ++e)
                if (!potato_peeling::cabello::segment_in_polygon(cleaned, e->source(),
                                                                 e->target())) {
                    contained = false;
                    break;
                }
            out["contained"] = contained;
            return out;
        },
        "convex_skull_exact(polygon) -> dict: vertices, vertices_exact "
        "(each coordinate as an exact 'num/den' rational), polygon_exact (the CLEANED "
        "input ring, same encoding), area (the exact rational rounded once), area_exact, "
        "and the EXACT predicates convex / contained, decided on the EPECK polygon "
        "before export, plus n_xi, a_entries and m_entries (None for a convex polygon).",
        nb::arg("polygon"));

    // Memory instrumentation (memprobe.hpp), for solve_memory.py. Off by default;
    // one solve per process.
    m.def("set_mem_probe", potato_peeling::chang_and_yap::memprobe::set_enabled,
          "set_mem_probe(on: bool) -> None: collect per-phase memory snapshots "
          "during convex_skull. Benchmarking only; one solve per process.",
          nb::arg("on"));

    m.def(
        "mem_report",
        []() {
            const auto &r = potato_peeling::chang_and_yap::memprobe::last();
            nb::list phases;
            for (const auto &s : r.snapshots) {
                nb::dict d;
                d["phase"] = s.phase;
                d["rss"] = s.rss_bytes;
                d["peak_rss"] = s.peak_rss_bytes;
                d["seconds"] = s.seconds;
                d["xi"] = s.xi;
                d["a_stored"] = s.a_stored;
                d["a_buckets"] = s.a_buckets;
                d["m_stored"] = s.m_stored;
                d["m_buckets"] = s.m_buckets;
                d["crit_l0"] = s.crit_l0;
                d["crit_l1"] = s.crit_l1;
                d["crit_ft"] = s.crit_ft;
                d["heap_live"] = s.heap_live;
                d["heap_freed"] = s.heap_freed;
                d["heap_arena"] = s.heap_arena;
                phases.append(d);
            }
            nb::dict out;
            out["n"] = r.n;
            out["baseline_rss"] = r.baseline_rss_bytes;
            out["peak_rss"] = r.peak_rss_bytes;
            out["phases"] = phases;
            return out;
        },
        "mem_report() -> dict: the last solve's memory report -- n, baseline_rss, "
        "peak_rss (VmHWM, what the OOM killer acts on), and a per-phase list "
        "carrying rss/peak_rss plus |Xi| and the stored-vs-bucket counts of A and M.");
}

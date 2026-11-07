#pragma once

// nanobind type casters for the CGAL types.
//
// These MUST be visible in every translation unit that binds a function taking or
// returning one of them: nanobind picks the caster when `m.def(...)` is
// instantiated, which happens in the TU holding the bind_* function, not in
// module.cpp. A specialization defined only in module.cpp silently leaves the
// other TUs using the generic (opaque class) caster, and the resulting functions
// reject plain Python tuples/lists at call time with "incompatible function
// arguments". Hence: one shared header, included by every binding TU.
//
// Two nanobind rules shape everything below and are worth stating once:
//
//  * `from_python` / `from_cpp` are `noexcept`. Nothing here may throw, so the
//    conversions delegate to nanobind's own sub-casters (which report failure by
//    returning false / an invalid handle) rather than to the throwing
//    `nb::cast<T>()` front end.
//  * A Point and a Segment cross the boundary as a TUPLE, not a list. The Python
//    harnesses compare them against tuples and use them as dict keys, so the
//    std::pair caster -- which produces a tuple -- is the right delegate; the
//    std::array one would produce a list.

#include <cstdint>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>

#include "potato_peeling_core/types.hpp"

namespace nanobind::detail {
    template <> struct type_caster<potato_peeling::FT> {
        NB_TYPE_CASTER(potato_peeling::FT, const_name("float"))

        bool from_python(handle src, std::uint32_t flags, cleanup_list *cleanup) noexcept {
            make_caster<double> caster;
            if (!caster.from_python(src, flags_for_local_caster<double>(flags), cleanup))
                return false;
            value = potato_peeling::FT(caster.operator cast_t<double>());
            return true;
        }

        static handle from_cpp(const potato_peeling::FT &src, const rv_policy policy,
                               cleanup_list *cleanup) noexcept {
            return make_caster<double>::from_cpp(CGAL::to_double(src), policy, cleanup);
        }
    };

    template <> struct type_caster<potato_peeling::Point> {
        NB_TYPE_CASTER(potato_peeling::Point, const_name("tuple[float, float]"))

        using Pair = std::pair<double, double>;
        using PairCaster = make_caster<Pair>;

        // Python (x, y) tuple → CGAL Point
        bool from_python(handle src, std::uint32_t flags, cleanup_list *cleanup) noexcept {
            PairCaster caster;
            if (!caster.from_python(src, flags_for_local_caster<Pair>(flags), cleanup))
                return false;
            const Pair xy = caster.operator cast_t<Pair>();
            value = potato_peeling::Point(xy.first, xy.second);
            return true;
        }

        // CGAL Point → Python (x, y) tuple
        static handle from_cpp(const potato_peeling::Point &src, const rv_policy policy,
                               cleanup_list *cleanup) noexcept {
            return PairCaster::from_cpp(Pair(CGAL::to_double(src.x()), CGAL::to_double(src.y())),
                                        policy, cleanup);
        }
    };

    template <> struct type_caster<potato_peeling::Segment> {
        NB_TYPE_CASTER(potato_peeling::Segment,
                       const_name("tuple[tuple[float, float], tuple[float, float]]"))

        using Pair = std::pair<potato_peeling::Point, potato_peeling::Point>;
        using PairCaster = make_caster<Pair>;

        // Python ((x1, y1), (x2, y2)) tuple → Segment
        bool from_python(handle src, std::uint32_t flags, cleanup_list *cleanup) noexcept {
            PairCaster caster;
            if (!caster.from_python(src, flags_for_local_caster<Pair>(flags), cleanup))
                return false;
            const Pair pq = caster.operator cast_t<Pair>();
            value = potato_peeling::Segment(pq.first, pq.second);
            return true;
        }

        // Segment → Python ((x1, y1), (x2, y2)) tuple
        static handle from_cpp(const potato_peeling::Segment &src, const rv_policy policy,
                               cleanup_list *cleanup) noexcept {
            return PairCaster::from_cpp(Pair(src.start(), src.end()), policy, cleanup);
        }
    };

    template <> struct type_caster<potato_peeling::Polygon> {
        // io_name: a polygon ARRIVES as any sequence of (x, y) pairs and LEAVES
        // as a list of tuples. Encoding the asymmetry here is what makes the
        // generated .pyi say `Sequence[...]` for a parameter and `list[...]` for
        // a return value, instead of overstating what callers may pass.
        NB_TYPE_CASTER(potato_peeling::Polygon,
                       io_name("Sequence[tuple[float, float]]", "list[tuple[float, float]]"))

        using Points = std::vector<potato_peeling::Point>;
        using PointsCaster = make_caster<Points>;

        // Python list of (x, y) tuples → CGAL Polygon
        bool from_python(handle src, std::uint32_t flags, cleanup_list *cleanup) noexcept {
            PointsCaster caster;
            if (!caster.from_python(src, flags_for_local_caster<Points>(flags), cleanup))
                return false;
            const Points pts = caster.operator cast_t<Points>();
            value = potato_peeling::Polygon(pts.begin(), pts.end());
            return true;
        }

        // CGAL Polygon → Python list of (x, y) tuples
        static handle from_cpp(const potato_peeling::Polygon &src, const rv_policy policy,
                               cleanup_list *cleanup) noexcept {
            return PointsCaster::from_cpp(Points(src.vertices_begin(), src.vertices_end()), policy,
                                          cleanup);
        }
    };
}  // namespace nanobind::detail

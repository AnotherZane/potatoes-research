#pragma once

// The working kernel, and conversions to and from EPECK.
//
// The algorithm runs in EPICK without losing exactness where it matters. Every
// point it reasons about is an input vertex or a sampled point, and both are
// pairs of doubles, on which EPICK's predicates are exact. Its only inexact
// constructions are areas, which rank candidates but never decide legality, and
// the corners of the sampling parallelogram Gamma. The answer's vertices are
// sample points, so they are exact doubles. Certification of the answer runs in
// EPECK (certify.hpp).

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include "potato_peeling_core/types.hpp"

namespace potato_peeling::cabello {

    using CK = CGAL::Exact_predicates_inexact_constructions_kernel;

    using Pt = CK::Point_2;
    using Vec = CK::Vector_2;
    using Seg = CK::Segment_2;
    using Tri = CK::Triangle_2;
    using Poly = CGAL::Polygon_2<CK>;

    // EPECK -> EPICK. Lossy in general, exact for coordinates that are doubles,
    // which all input coordinates are (pre_checks only deletes vertices).
    using ToInexact = CGAL::Cartesian_converter<K, CK>;

    // EPICK -> EPECK. Always exact.
    using ToExact = CGAL::Cartesian_converter<CK, K>;

    inline Pt to_inexact(const Point &p) {
        return ToInexact()(p);
    }

    inline Point to_exact(const Pt &p) {
        return ToExact()(p);
    }

    Poly to_inexact(const Polygon &poly);

}  // namespace potato_peeling::cabello

#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Aff_transformation_2.h>

namespace potato_peeling {
    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    using Vector = K::Vector_2;
    using Point = K::Point_2;
    using Line = K::Line_2;
    using FT = K::FT;
    using Polygon = CGAL::Polygon_2<K>;
    using Polygon_with_holes = CGAL::Polygon_with_holes_2<K>;
    using Aff_transformation = CGAL::Aff_transformation_2<K>;
    using Segment = K::Segment_2;
    using Direction = K::Direction_2;
}  // namespace potato_peeling

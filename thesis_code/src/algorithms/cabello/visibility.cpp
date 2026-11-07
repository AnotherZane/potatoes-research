#include "potato_peeling_core/algorithms/cabello/visibility.hpp"

#include <cstddef>
#include <list>

#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_tree.h>

namespace potato_peeling::cabello {
    namespace {
        using Iterator = std::list<Seg>::const_iterator;
        using Primitive = CGAL::AABB_segment_primitive_2<CK, Iterator>;
        using Traits = CGAL::AABB_traits_2<CK, Primitive>;
        using Tree = CGAL::AABB_tree<Traits>;
    }  // namespace

    struct VisibilityOracle::Impl {
        // A list, not a vector: AABB_segment_primitive_2 stores iterators, and a
        // vector's would dangle on any reallocation.
        std::list<Seg> edges;
        Tree tree;
    };

    VisibilityOracle::VisibilityOracle(const Poly &p) : poly_(p), impl_(std::make_unique<Impl>()) {
        for (auto e = p.edges_begin(); e != p.edges_end(); ++e) impl_->edges.push_back(*e);
        impl_->tree.insert(impl_->edges.begin(), impl_->edges.end());
        impl_->tree.build();
        impl_->tree.accelerate_distance_queries();
    }

    VisibilityOracle::VisibilityOracle(VisibilityOracle &&) noexcept = default;

    VisibilityOracle &VisibilityOracle::operator=(VisibilityOracle &&) noexcept = default;

    VisibilityOracle::~VisibilityOracle() = default;

    bool VisibilityOracle::inside(const Pt &q) const {
        // Strictly inside; on the boundary is not inside. Exact under EPICK for
        // double coordinates -- bounded_side is a predicate, not a construction.
        return poly_.bounded_side(q) == CGAL::ON_BOUNDED_SIDE;
    }

    bool VisibilityOracle::visible(const Pt &a, const Pt &b) const {
        if (!inside(a) || !inside(b)) return false;
        if (a == b) return true;

        // do_intersect only: no intersection point is ever constructed, so this
        // stays exact on double coordinates.
        return !impl_->tree.do_intersect(Seg(a, b));
    }

    PointGraph visibility_graph(const VisibilityOracle &oracle, const std::vector<Pt> &points) {
        PointGraph g(points);
        for (std::size_t i = 0; i < points.size(); ++i)
            for (std::size_t j = i + 1; j < points.size(); ++j)
                if (oracle.visible(points[i], points[j])) g.add_edge(i, j);
        return g;
    }

}  // namespace potato_peeling::cabello

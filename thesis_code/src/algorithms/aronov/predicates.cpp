#include "potato_peeling_core/algorithms/aronov/predicates.hpp"

namespace potato_peeling::aronov {
    FT edge_y_at(const Mesh &mesh, const std::size_t e, const FT &x) {
        const Point &l = mesh.left_point(e);
        const Point &r = mesh.right_point(e);
        return l.y() + (r.y() - l.y()) * (x - l.x()) / (r.x() - l.x());
    }

    namespace {
        // Arrive at v along e_in, leave along e_out.
        //
        // The bottom chain runs left-to-right and the top chain right-to-left, but
        // both are traversed in CCW boundary order, so the convexity test is the
        // same predicate in both roles — only the E_L / E_R membership of the two
        // arguments differs. The paper's <(e_in, e_out) >= 0 admits the straight
        // case, hence "not a right turn" rather than "is a left turn".
        bool no_right_turn(const Mesh &mesh, const std::size_t e_in, const std::size_t v,
                           const std::size_t e_out) {
            return CGAL::orientation(mesh.point(mesh.other(e_in, v)), mesh.point(v),
                                     mesh.point(mesh.other(e_out, v)))
                   != CGAL::RIGHT_TURN;
        }
    }  // namespace

    std::optional<std::pair<FT, FT>> overlap(const Mesh &mesh, const std::size_t e,
                                             const std::size_t e_prime) {
        // By value: Epeck's Point_2::x() returns an FT by value, so binding these to
        // const references would rely on lifetime extension and is easy to get wrong.
        const FT le = mesh.left_point(e).x();
        const FT lp = mesh.left_point(e_prime).x();
        const FT re = mesh.right_point(e).x();
        const FT rp = mesh.right_point(e_prime).x();

        const FT a = (le > lp) ? le : lp;
        const FT b = (re < rp) ? re : rp;

        // A single-point touch carries no area and is never a usable pair.
        if (a < b) return std::make_pair(a, b);
        return std::nullopt;
    }

    FT trapezoid_area(const Mesh &mesh, const std::size_t e, const std::size_t e_prime) {
        const auto span = overlap(mesh, e, e_prime);
        if (!span) return FT(0);

        const auto &[a, b] = *span;
        const FT height_a = edge_y_at(mesh, e, a) - edge_y_at(mesh, e_prime, a);
        const FT height_b = edge_y_at(mesh, e, b) - edge_y_at(mesh, e_prime, b);
        return (b - a) * (height_a + height_b) / FT(2);
    }

    bool pair_precondition(const Mesh &mesh, const std::size_t e, const std::size_t e_prime) {
        return mesh.face_below(e) != NO_FACE && mesh.face_above(e_prime) != NO_FACE;
    }

    std::size_t bottom_valid_prefix(const Mesh &mesh, const std::size_t e) {
        const std::size_t v = mesh.left(e);
        const auto &el = mesh.vertex_left_edges[v];

        // E_L is ordered topmost first, i.e. by increasing direction angle in
        // (pi/2, 3pi/2). Writing d_in = v - other(e'') (which points rightwards, at
        // angle psi - pi) and d_out for e, the condition cross(d_in, d_out) >= 0
        // reads psi <= angle(d_out) + pi — a leading run of that order.
        std::size_t k = 0;
        while (k < el.size() && no_right_turn(mesh, el[k], v, e)) ++k;
        return k;
    }

    std::size_t top_valid_suffix(const Mesh &mesh, const std::size_t e) {
        const std::size_t v = mesh.left(e);
        const auto &el = mesh.vertex_left_edges[v];

        // Mirror of the above: here both directions point leftwards, so the
        // condition reads psi >= angle(d_in) — a trailing run of the same order.
        std::size_t k = 0;
        while (k < el.size() && !no_right_turn(mesh, e, v, el[k])) ++k;
        return k;
    }
}  // namespace potato_peeling::aronov

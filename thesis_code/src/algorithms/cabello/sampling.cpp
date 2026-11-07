#include "potato_peeling_core/algorithms/cabello/sampling.hpp"
#include "potato_peeling_core/algorithms/cabello/rng.hpp"

#include <algorithm>
#include <cmath>
#include <list>
#include <stdexcept>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

namespace potato_peeling::cabello {
    namespace {
        struct FaceInfo {
            int nesting_level{-1};

            bool in_domain() const { return nesting_level % 2 == 1; }
        };

        using Vb = CGAL::Triangulation_vertex_base_2<CK>;
        using Fbb = CGAL::Triangulation_face_base_with_info_2<FaceInfo, CK>;
        using Fb = CGAL::Constrained_triangulation_face_base_2<CK, Fbb>;
        using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
        using CDT = CGAL::Constrained_Delaunay_triangulation_2<CK, TDS, CGAL::Exact_predicates_tag>;

        // Flood fill from the infinite face, incrementing the nesting level at
        // every constrained crossing; odd levels are inside P. The polygon is
        // simple and every constraint lies on its boundary, so parity is exactly
        // inside/outside.
        void mark(CDT &ct, const CDT::Face_handle start, const int index,
                  std::list<CDT::Edge> &border) {
            if (start->info().nesting_level != -1) return;

            std::list<CDT::Face_handle> queue{start};
            while (!queue.empty()) {
                const CDT::Face_handle fh = queue.front();
                queue.pop_front();
                if (fh->info().nesting_level != -1) continue;

                fh->info().nesting_level = index;
                for (int i = 0; i < 3; ++i) {
                    if (const CDT::Edge e(fh, i); fh->neighbor(i)->info().nesting_level == -1) {
                        if (ct.is_constrained(e))
                            border.push_back(e);
                        else
                            queue.push_back(fh->neighbor(i));
                    }
                }
            }
        }

        void mark_domains(CDT &cdt) {
            for (auto f = cdt.all_faces_begin(); f != cdt.all_faces_end(); ++f)
                f->info().nesting_level = -1;

            std::list<CDT::Edge> border;
            mark(cdt, cdt.infinite_face(), 0, border);
            while (!border.empty()) {
                const CDT::Edge e = border.front();
                border.pop_front();
                if (e.first->neighbor(e.second)->info().nesting_level == -1)
                    mark(cdt, e.first->neighbor(e.second), e.first->info().nesting_level + 1,
                         border);
            }
        }
    }  // namespace

    Triangulation triangulate(const Poly &p) {
        Triangulation out;
        if (p.size() < 3) return out;

        CDT cdt;
        cdt.insert_constraint(p.vertices_begin(), p.vertices_end(), true);
        mark_domains(cdt);

        double running = 0;
        for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
            if (!f->info().in_domain()) continue;
            const std::array<Pt, 3> t{f->vertex(0)->point(), f->vertex(1)->point(),
                                      f->vertex(2)->point()};
            const double a = std::abs(CGAL::to_double(CGAL::area(t[0], t[1], t[2])));
            if (a <= 0) continue;  // a zero-area face carries no probability mass
            running += a;
            out.triangles.push_back(t);
            out.prefix.push_back(running);
        }
        out.total_area = running;
        return out;
    }

    Pt sample_in_triangle(const std::array<Pt, 3> &t, std::mt19937_64 &gen) {
        // [C17 p14]: sample the parallelogram spanned by two edges, which holds two
        // congruent copies of the triangle, and reflect the far copy back.
        double u = canonical(gen);
        double v = canonical(gen);
        if (u + v > 1.0) {
            u = 1.0 - u;
            v = 1.0 - v;
        }
        const Vec e1 = t[1] - t[0];
        const Vec e2 = t[2] - t[0];
        return t[0] + u * e1 + v * e2;
    }

    Pt sample_uniform(const Triangulation &tri, std::mt19937_64 &gen) {
        if (tri.empty()) throw std::invalid_argument("sample_uniform: empty triangulation");

        const double x = canonical(gen) * tri.total_area;
        // Smallest j with x <= prefix[j], the paper's binary search.
        const auto it = std::lower_bound(tri.prefix.begin(), tri.prefix.end(), x);
        const std::size_t j = std::min<std::size_t>(
            static_cast<std::size_t>(it - tri.prefix.begin()), tri.size() - 1);
        return sample_in_triangle(tri.triangles[j], gen);
    }

    std::vector<Pt> sample_uniform(const Triangulation &tri, const std::size_t count,
                                   std::mt19937_64 &gen,
                                   const std::function<bool(const Pt &)> &inside,
                                   std::size_t *retries) {
        std::vector<Pt> out;
        out.reserve(count);

        // t0 + u*e1 + v*e2 lies in the triangle exactly, but its rounded value can
        // land just outside P. Such a point would fail every visibility test, so
        // it is re-drawn from the same substream.
        std::size_t wasted = 0;
        while (out.size() < count) {
            Pt p = sample_uniform(tri, gen);
            if (inside(p)) {
                out.push_back(p);
            } else if (++wasted > 64 + 4 * count) {
                throw std::runtime_error(
                    "sample_uniform: too many points fell outside P; the triangulation and "
                    "the containment predicate disagree");
            }
        }
        if (retries) *retries = wasted;
        return out;
    }

    std::optional<Parallelogram> gamma(const Pt &a, const Pt &b, const double A) {
        const double dy = CGAL::to_double(a.y() - b.y());
        if (dy == 0.0) return std::nullopt;  // Gamma is undefined

        // [C17 p8]. a' = 2a - b, b' = 2b - a, translated horizontally by
        // +-d with d = 2A/|y(a) - y(b)|.
        const double d = 2.0 * A / std::abs(dy);
        const Pt ap = a + (a - b);  // 2a - b
        const Pt bp = b + (b - a);  // 2b - a

        Parallelogram g;
        g.origin = Pt(ap.x() - d, ap.y());
        g.e1 = Vec(2.0 * d, 0.0);  // a'-(d,0) -> a'+(d,0)
        g.e2 = bp - ap;            // = 3(b - a)
        // |e1 x e2| = 2d * 3|y(b)-y(a)| = 12A, by construction.
        g.area = std::abs(CGAL::to_double(g.e1.x() * g.e2.y() - g.e1.y() * g.e2.x()));
        return g;
    }

    Pt sample_uniform(const Parallelogram &g, std::mt19937_64 &gen) {
        const double u = canonical(gen);
        const double v = canonical(gen);
        return g.origin + u * g.e1 + v * g.e2;
    }

}  // namespace potato_peeling::cabello

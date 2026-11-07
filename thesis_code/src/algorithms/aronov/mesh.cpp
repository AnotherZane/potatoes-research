#include "potato_peeling_core/algorithms/aronov/mesh.hpp"

#include <algorithm>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <utility>

#include "potato_peeling_core/algorithms/aronov/predicates.hpp"
#include "potato_peeling_core/geometry.hpp"

namespace potato_peeling::aronov {
    namespace {
        // ── Angular helpers ───────────────────────────────────────────────────
        //
        // All exact: only sign tests on cross products, no constructed angles.

        FT cross(const Vector &a, const Vector &b) {
            return a.x() * b.y() - a.y() * b.x();
        }

        // 0 for a direction with angle in [0, pi), 1 for [pi, 2pi).
        int angular_half(const Vector &d) {
            if (d.y() > FT(0)) return 0;
            if (d.y() < FT(0)) return 1;
            return d.x() > FT(0) ? 0 : 1;
        }

        // Strict CCW order by angle from the +x axis.
        bool angle_less(const Vector &a, const Vector &b) {
            const int ha = angular_half(a);
            if (const int hb = angular_half(b); ha != hb) return ha < hb;
            return cross(a, b) > FT(0);
        }

        // ── CDT with per-face nesting level, for domain marking ───────────────

        struct FaceInfo {
            int nesting_level{-1};

            [[nodiscard]] bool in_domain() const { return nesting_level % 2 == 1; }
        };

        using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K>;
        using Fbb = CGAL::Triangulation_face_base_with_info_2<FaceInfo, K>;
        using Fbc = CGAL::Constrained_triangulation_face_base_2<K, Fbb>;
        using TDS = CGAL::Triangulation_data_structure_2<Vb, Fbc>;
        using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS>;

        void mark_domains(CDT &ct, const CDT::Face_handle start, const int index,
                          std::list<CDT::Edge> &border) {
            if (start->info().nesting_level != -1) return;

            std::list<CDT::Face_handle> queue{start};
            while (!queue.empty()) {
                const CDT::Face_handle fh = queue.front();
                queue.pop_front();
                if (fh->info().nesting_level != -1) continue;

                fh->info().nesting_level = index;
                for (int i = 0; i < 3; ++i) {
                    const CDT::Edge e(fh, i);
                    if (const CDT::Face_handle n = fh->neighbor(i); n->info().nesting_level == -1) {
                        if (ct.is_constrained(e))
                            border.push_back(e);
                        else
                            queue.push_back(n);
                    }
                }
            }
        }

        // Flood fill from the infinite face, incrementing the nesting level on every
        // constrained crossing; odd levels are inside P. Handles holes unchanged.
        // Valid because every constraint lies on dP.
        void mark_domains(CDT &cdt) {
            for (auto f = cdt.all_faces_begin(); f != cdt.all_faces_end(); ++f)
                f->info().nesting_level = -1;

            std::list<CDT::Edge> border;
            mark_domains(cdt, cdt.infinite_face(), 0, border);
            while (!border.empty()) {
                const CDT::Edge e = border.front();
                border.pop_front();
                if (const CDT::Face_handle n = e.first->neighbor(e.second);
                    n->info().nesting_level == -1)
                    mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
            }
        }

        [[noreturn]] void fail(const std::string &what) {
            throw std::invalid_argument("aronov::build_mesh: " + what);
        }

        // ── Shared construction ──────────────────────────────────────────────

        Mesh finalize(const std::vector<Point> &original_points,
                      std::vector<std::array<std::size_t, 3>> faces) {
            Mesh mesh;

            const std::size_t n = original_points.size();
            if (n < 3) fail("a mesh needs at least 3 vertices");
            if (faces.empty()) fail("a mesh needs at least one triangle");

            // 1. Sweep normalisation. det == 1, so areas are untouched.
            const FT s = separating_shear(original_points);
            mesh.shear = Aff_transformation(1, s, 0, 0, 1, 0);
            mesh.unshear = Aff_transformation(1, -s, 0, 0, 1, 0);

            std::vector<Point> sheared;
            sheared.reserve(n);
            for (const Point &p : original_points) sheared.push_back(mesh.shear(p));

            // 2. Order vertices by x. separating_shear guarantees strictness.
            std::vector<std::size_t> order(n);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
                return sheared[a].x() < sheared[b].x();
            });

            std::vector<std::size_t> rank(n);
            mesh.vertices.reserve(n);
            for (std::size_t k = 0; k < n; ++k) {
                rank[order[k]] = k;
                mesh.vertices.push_back(sheared[order[k]]);
            }
            for (std::size_t k = 1; k < n; ++k)
                if (!(mesh.vertices[k - 1].x() < mesh.vertices[k].x()))
                    fail("shear failed to separate x-coordinates (internal error)");

            // 3. Remap, orient CCW, and measure the faces.
            std::set<std::array<std::size_t, 3>> seen;
            mesh.faces.reserve(faces.size());
            mesh.face_area.reserve(faces.size());
            for (auto &f : faces) {
                for (auto &v : f) {
                    if (v >= n) fail("face references a vertex id out of range");
                    v = rank[v];
                }
                if (f[0] == f[1] || f[1] == f[2] || f[0] == f[2]) fail("face repeats a vertex");

                const Point &a = mesh.vertices[f[0]];
                const Point &b = mesh.vertices[f[1]];
                const Point &c = mesh.vertices[f[2]];

                const CGAL::Orientation o = CGAL::orientation(a, b, c);
                if (o == CGAL::COLLINEAR) fail("face has zero area");
                if (o == CGAL::CLOCKWISE) std::swap(f[1], f[2]);

                auto key = f;
                std::sort(key.begin(), key.end());
                if (!seen.insert(key).second) fail("duplicate face");

                mesh.faces.push_back(f);
                const Point &p0 = mesh.vertices[f[0]];
                mesh.face_area.push_back(
                    CGAL::abs(cross(mesh.vertices[f[1]] - p0, mesh.vertices[f[2]] - p0)) / FT(2));
            }

            // 4. Edges, with {above, below} taken relative to left -> right.
            std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge_id;
            for (std::size_t fi = 0; fi < mesh.faces.size(); ++fi) {
                const auto &f = mesh.faces[fi];
                for (int k = 0; k < 3; ++k) {
                    const std::size_t u = f[k];
                    const std::size_t w = f[(k + 1) % 3];
                    const std::size_t third = f[(k + 2) % 3];

                    const auto key = std::minmax(u, w);
                    auto [it, inserted] = edge_id.emplace(key, mesh.edges.size());
                    if (inserted) {
                        // Vertices are sorted by x, so the smaller id is the left one.
                        mesh.edges.push_back({key.first, key.second});
                        mesh.edge_faces.push_back({NO_FACE, NO_FACE});
                    }
                    const std::size_t e = it->second;

                    // Left of the directed segment left -> right is above it.
                    const CGAL::Orientation o
                        = CGAL::orientation(mesh.vertices[mesh.edges[e][0]],
                                            mesh.vertices[mesh.edges[e][1]], mesh.vertices[third]);
                    const int side = (o == CGAL::LEFT_TURN) ? 0 : 1;

                    if (mesh.edge_faces[e][side] != NO_FACE)
                        fail("two faces lie on the same side of an edge (overlapping triangles)");
                    mesh.edge_faces[e][side] = fi;
                }
            }

            // 5. Incident edges per vertex: CCW overall, then E_L / E_R topmost first.
            mesh.vertex_edges.assign(n, {});
            for (std::size_t e = 0; e < mesh.edges.size(); ++e) {
                mesh.vertex_edges[mesh.edges[e][0]].push_back(e);
                mesh.vertex_edges[mesh.edges[e][1]].push_back(e);
            }

            mesh.vertex_left_edges.assign(n, {});
            mesh.vertex_right_edges.assign(n, {});
            for (std::size_t v = 0; v < n; ++v) {
                auto dir = [&](const std::size_t e) {
                    return mesh.vertices[mesh.other(e, v)] - mesh.vertices[v];
                };

                auto &around = mesh.vertex_edges[v];
                std::sort(around.begin(), around.end(),
                          [&](const std::size_t a, const std::size_t b) {
                              return angle_less(dir(a), dir(b));
                          });

                // E_L directions all lie in (pi/2, 3pi/2), so the CCW order already
                // runs from the topmost edge downwards — exactly the order the
                // paper's list L(e) wants (p357).
                for (const std::size_t e : around) {
                    if (mesh.other(e, v) < v)
                        mesh.vertex_left_edges[v].push_back(e);
                    else
                        mesh.vertex_right_edges[v].push_back(e);
                }

                // E_R directions straddle angle 0, so sort them explicitly. Both lie
                // in the right half-plane, so a single cross product orders any pair.
                auto &er = mesh.vertex_right_edges[v];
                std::sort(er.begin(), er.end(), [&](const std::size_t a, const std::size_t b) {
                    return cross(dir(b), dir(a)) > FT(0);  // a is above b
                });
            }

            // 6. opp() (p356), now that E_L / E_R are in place.
            mesh.bottom_prefix_len.assign(mesh.edges.size(), 0);
            mesh.top_suffix_begin.assign(mesh.edges.size(), 0);
            for (std::size_t e = 0; e < mesh.edges.size(); ++e) {
                mesh.bottom_prefix_len[e] = bottom_valid_prefix(mesh, e);
                mesh.top_suffix_begin[e] = top_valid_suffix(mesh, e);
            }

            return mesh;
        }
    }  // namespace

    // ── Mesh accessors ───────────────────────────────────────────────────────

    std::size_t Mesh::other(const std::size_t e, const std::size_t v) const {
        if (edges[e][0] == v) return edges[e][1];
        if (edges[e][1] == v) return edges[e][0];
        throw std::invalid_argument("aronov::Mesh::other: v is not an endpoint of e");
    }

    std::size_t Mesh::find_edge(const std::size_t u, const std::size_t v) const {
        for (const std::size_t e : vertex_edges[u])
            if (other(e, u) == v) return e;
        return NO_FACE;
    }

    FT Mesh::total_area() const {
        FT total(0);
        for (const FT &a : face_area) total += a;
        return total;
    }

    // ── Shear ────────────────────────────────────────────────────────────────

    FT separating_shear(const std::vector<Point> &points) {
        const std::size_t n = points.size();

        // Fast path: the input already satisfies the paper's general-position
        // assumption, so the identity works and the O(n^2) scan is skipped.
        {
            std::vector<std::size_t> order(n);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
                return points[a].x() < points[b].x();
            });
            bool distinct = true;
            for (std::size_t k = 1; k < n && distinct; ++k)
                if (points[order[k - 1]].x() == points[order[k]].x()) distinct = false;
            if (distinct) return FT(0);
        }

        // x_i + s*y_i == x_j + s*y_j  <=>  s == (x_i - x_j) / (y_j - y_i).
        // Collect every s that collides, then pick one that does not.
        std::vector<FT> bad;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const FT dy = points[j].y() - points[i].y();
                if (dy == FT(0)) {
                    if (points[i].x() == points[j].x())
                        throw std::invalid_argument("aronov::separating_shear: duplicate vertex");
                    continue;  // no shear ever collides these two
                }
                if (const FT s = (points[i].x() - points[j].x()) / dy; s >= FT(0)) bad.push_back(s);
            }
        }

        std::sort(bad.begin(), bad.end());
        bad.erase(std::unique(bad.begin(), bad.end()), bad.end());

        if (bad.empty() || bad.front() != FT(0)) return FT(0);
        // bad is sorted, distinct and starts at 0; the midpoint of the first gap is
        // strictly between two consecutive bad values, hence not itself bad.
        if (bad.size() == 1) return FT(1);
        return bad[1] / FT(2);
    }

    // ── Builders ─────────────────────────────────────────────────────────────

    Mesh build_mesh(const Polygon &outer) {
        Polygon poly = outer;

        if (!poly.is_simple()) fail("input polygon must be simple");
        remove_collinear_vertices(poly);
        if (!poly.is_simple()) fail("input polygon must remain simple after cleanup");
        if (poly.size() < 3) fail("input polygon must have at least 3 vertices");
        if (poly.is_clockwise_oriented()) poly.reverse_orientation();

        // Triangulate in the caller's frame, so the mesh does not depend on our
        // internal choice of shear; the shear is applied afterwards by finalize and
        // changes no combinatorics.
        const std::vector<Point> boundary(poly.vertices_begin(), poly.vertices_end());
        const std::size_t n = boundary.size();

        CDT cdt;
        std::vector<CDT::Vertex_handle> handles;
        handles.reserve(n);
        for (const Point &p : boundary) handles.push_back(cdt.insert(p));
        for (std::size_t i = 0; i < n; ++i) cdt.insert_constraint(handles[i], handles[(i + 1) % n]);

        if (cdt.number_of_vertices() != n)
            fail("triangulation changed the vertex set (degenerate input polygon)");

        mark_domains(cdt);

        // Index every finite vertex, keeping only those used by an in-domain face.
        std::vector<std::array<std::size_t, 3>> raw_faces;
        std::size_t next = 0;
        for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v)
            v->info() = next++;
        std::vector<char> used(next, 0);
        for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
            if (!f->info().in_domain()) continue;
            const std::array<std::size_t, 3> t{f->vertex(0)->info(), f->vertex(1)->info(),
                                               f->vertex(2)->info()};
            for (const std::size_t i : t) used[i] = 1;
            raw_faces.push_back(t);
        }

        std::vector<Point> points;
        std::vector<std::size_t> remap(next, 0);
        for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
            const std::size_t id = v->info();
            if (!used[id]) continue;
            remap[id] = points.size();
            points.push_back(v->point());
        }
        std::vector<std::array<std::size_t, 3>> faces;
        faces.reserve(raw_faces.size());
        for (auto &t : raw_faces) faces.push_back({remap[t[0]], remap[t[1]], remap[t[2]]});

        return finalize(points, std::move(faces));
    }

    // ── Validation ───────────────────────────────────────────────────────────

    MeshCheck validate(const Mesh &mesh, const std::optional<FT> &expected_area) {
        std::ostringstream problems;
        int count = 0;
        auto report = [&](const std::string &what) {
            if (count++) problems << "; ";
            problems << what;
        };

        // The shear must be area-preserving, or every area below is meaningless.
        if (mesh.shear.m(0, 0) * mesh.shear.m(1, 1) - mesh.shear.m(0, 1) * mesh.shear.m(1, 0)
            != FT(1))
            report("shear is not area-preserving");

        for (std::size_t v = 1; v < mesh.num_vertices(); ++v)
            if (!(mesh.vertices[v - 1].x() < mesh.vertices[v].x())) {
                report("vertices are not strictly ordered by x at id " + std::to_string(v));
                break;
            }

        for (std::size_t f = 0; f < mesh.num_faces(); ++f) {
            const auto &t = mesh.faces[f];
            if (CGAL::orientation(mesh.vertices[t[0]], mesh.vertices[t[1]], mesh.vertices[t[2]])
                != CGAL::LEFT_TURN) {
                report("face " + std::to_string(f) + " is not CCW");
                break;
            }
        }

        // Every edge borders at least one triangle, and each face is registered
        // exactly once against each of its three edges.
        std::vector<int> incidence(mesh.num_faces(), 0);
        for (std::size_t e = 0; e < mesh.num_edges(); ++e) {
            const std::size_t above = mesh.face_above(e), below = mesh.face_below(e);
            if (above == NO_FACE && below == NO_FACE) {
                report("edge " + std::to_string(e) + " borders no triangle");
                break;
            }
            if (above != NO_FACE) ++incidence[above];
            if (below != NO_FACE) ++incidence[below];
        }
        for (std::size_t f = 0; f < mesh.num_faces(); ++f)
            if (incidence[f] != 3) {
                report("face " + std::to_string(f) + " has " + std::to_string(incidence[f])
                       + " edge incidences, expected 3");
                break;
            }

        // Boundary edges (exactly one triangle) must form closed cycles, so every
        // vertex carries an even number of them. True with holes as well.
        std::vector<int> boundary_degree(mesh.num_vertices(), 0);
        for (std::size_t e = 0; e < mesh.num_edges(); ++e)
            if ((mesh.face_above(e) == NO_FACE) != (mesh.face_below(e) == NO_FACE)) {
                ++boundary_degree[mesh.edges[e][0]];
                ++boundary_degree[mesh.edges[e][1]];
            }
        for (std::size_t v = 0; v < mesh.num_vertices(); ++v)
            if (boundary_degree[v] % 2 != 0) {
                report("vertex " + std::to_string(v) + " has an odd number of boundary edges");
                break;
            }

        // E_L and E_R must partition the incident edges of every vertex.
        for (std::size_t v = 0; v < mesh.num_vertices(); ++v)
            if (mesh.vertex_left_edges[v].size() + mesh.vertex_right_edges[v].size()
                != mesh.vertex_edges[v].size()) {
                report("E_L / E_R do not partition the incident edges of vertex "
                       + std::to_string(v));
                break;
            }

        // M is assumed connected (p350); the DP silently relies on it.
        if (mesh.num_faces() > 0) {
            std::vector<bool> seen(mesh.num_faces(), false);
            std::vector<std::size_t> stack{0};
            seen[0] = true;
            std::size_t reached = 1;
            while (!stack.empty()) {
                const std::size_t f = stack.back();
                stack.pop_back();
                const auto &t = mesh.faces[f];
                for (int k = 0; k < 3; ++k) {
                    const std::size_t e = mesh.find_edge(t[k], t[(k + 1) % 3]);
                    // validate() is the safety net, so it must survive a mesh whose
                    // edge table does not cover its own faces rather than index NO_FACE.
                    if (e == NO_FACE) continue;
                    for (const std::size_t nb : {mesh.face_above(e), mesh.face_below(e)})
                        if (nb != NO_FACE && !seen[nb]) {
                            seen[nb] = true;
                            ++reached;
                            stack.push_back(nb);
                        }
                }
            }
            if (reached != mesh.num_faces())
                report("mesh is not connected (" + std::to_string(reached) + " of "
                       + std::to_string(mesh.num_faces()) + " faces reachable)");
        }

        // T-vertex check: a vertex inside a neighbouring triangle's edge passes every
        // check above, exact area included, but makes the DP treat that edge as
        // exterior and return a suboptimal answer.
        //
        // Vertex ids are the x-order, so the vertices strictly inside an edge's
        // x-span are the ids strictly between its endpoints, and one collinear with
        // the edge lies on the open segment. Sign tests only.
        for (std::size_t e = 0; e < mesh.num_edges() && count == 0; ++e) {
            const std::size_t lo = mesh.left(e), hi = mesh.right(e);
            for (std::size_t v = lo + 1; v < hi; ++v)
                if (CGAL::orientation(mesh.left_point(e), mesh.right_point(e), mesh.vertices[v])
                    == CGAL::COLLINEAR) {
                    report("edge " + std::to_string(e) + " has vertex " + std::to_string(v)
                           + " in its interior (T-vertex: the mesh is not conforming)");
                    break;
                }
        }

        if (expected_area && mesh.total_area() != *expected_area)
            report("mesh area does not equal the expected area exactly");

        if (count == 0) return {};
        return {false, problems.str()};
    }
}  // namespace potato_peeling::aronov

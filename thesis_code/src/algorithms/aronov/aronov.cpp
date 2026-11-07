#include "potato_peeling_core/algorithms/aronov/aronov.hpp"

#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>

#include "potato_peeling_core/algorithms/aronov/sweep.hpp"
#include "potato_peeling_core/geometry.hpp"

namespace potato_peeling::aronov {

    Polygon pre_checks(const Polygon &input) {
        Polygon poly = input;

        if (!poly.is_simple()) {
            throw std::invalid_argument("Input polygon must be simple.");
        }

        remove_collinear_vertices(poly);

        if (!poly.is_simple()) {
            throw std::invalid_argument("Input polygon must remain simple after cleanup.");
        }

        return poly;
    }

    std::optional<EdgePair> best_pair(const Mesh &mesh, const QTable &q) {
        std::optional<EdgePair> winner;
        std::optional<FT> best;

        for (const auto &[pair, entry] : q) {
            // Step III (p354): e and e' must share their right endpoint, so that the
            // clipping line x = max X(e, e') passes through it and U(C) is the whole
            // union rather than a clipped one.
            if (mesh.right(pair.first) != mesh.right(pair.second)) continue;
            if (!best || entry.area > *best) {
                best = entry.area;
                winner = pair;
            }
        }
        return winner;
    }

    MeshedSkull reconstruct(const Mesh &mesh, const QTable &q, const EdgePair winner) {
        // Walk the backpointers from the winning pair down to the base case,
        // peeling one edge off whichever chain the DP extended at each step.
        std::vector<std::size_t> top, bottom;  // both collected right to left
        EdgePair cur = winner;
        for (;;) {
            const auto it = q.find(cur);
            if (it == q.end())
                throw std::runtime_error("aronov::reconstruct: dangling backpointer");

            if (!it->second.from) {
                // Base case: the chains meet, so both edges belong to the region.
                top.push_back(cur.first);
                bottom.push_back(cur.second);
                break;
            }

            const EdgePair sub = *it->second.from;
            // The sub-pair keeps one edge and replaces the other; the replaced one is
            // the edge this step consumed.
            if (sub.first == cur.first)
                bottom.push_back(cur.second);
            else
                top.push_back(cur.first);
            cur = sub;
        }

        std::reverse(top.begin(), top.end());  // now left to right
        std::reverse(bottom.begin(), bottom.end());

        MeshedSkull skull;
        skull.winner = winner;
        skull.area = q.find(winner)->second.area;

        // CCW boundary: leftmost vertex, along the bottom chain to the rightmost
        // vertex, then back along the top chain. The top chain's leftmost vertex is
        // the same point we started from, so it is not repeated.
        Polygon sheared;
        sheared.push_back(mesh.point(mesh.left(bottom.front())));
        for (const std::size_t e : bottom) sheared.push_back(mesh.point(mesh.right(e)));
        for (std::size_t i = top.size(); i-- > 1;) sheared.push_back(mesh.point(mesh.left(top[i])));

        for (auto v = sheared.vertices_begin(); v != sheared.vertices_end(); ++v)
            skull.boundary.push_back(mesh.to_original(*v));

        // The triangles: flood fill inwards from just below the top chain, blocked by
        // the boundary edges. Cheaper and more direct than a point-in-polygon sweep,
        // and it fails loudly (below) if the chains do not actually enclose a region.
        std::set<std::size_t> boundary_edges(top.begin(), top.end());
        boundary_edges.insert(bottom.begin(), bottom.end());

        const std::size_t seed = mesh.face_below(top.back());
        if (seed != NO_FACE) {
            std::vector<bool> seen(mesh.num_faces(), false);
            std::vector<std::size_t> stack{seed};
            seen[seed] = true;
            while (!stack.empty()) {
                const std::size_t f = stack.back();
                stack.pop_back();
                skull.faces.push_back(f);

                const auto &t = mesh.faces[f];
                for (int k = 0; k < 3; ++k) {
                    const std::size_t e = mesh.find_edge(t[k], t[(k + 1) % 3]);
                    if (e == NO_FACE || boundary_edges.count(e)) continue;
                    const std::size_t nb
                        = (mesh.face_above(e) == f) ? mesh.face_below(e) : mesh.face_above(e);
                    // Escaping through a non-boundary edge into the exterior means the
                    // chains did not enclose the region; check_skull then reports an
                    // area mismatch (O1).
                    if (nb == NO_FACE || seen[nb]) continue;
                    seen[nb] = true;
                    stack.push_back(nb);
                }
            }
            std::sort(skull.faces.begin(), skull.faces.end());
        }

        return skull;
    }

    SkullCheck check_skull(const Mesh &mesh, const MeshedSkull &skull) {
        std::ostringstream problems;
        int count = 0;
        auto report = [&](const std::string &what) {
            if (count++) problems << "; ";
            problems << what;
        };

        if (!skull.winner) {
            report("no winning pair");
            return {false, problems.str()};
        }

        // O1 — the DP value and the reconstruction must agree EXACTLY. A divergence
        // is always a bug, never a tolerance issue.
        FT reconstructed(0);
        for (const std::size_t f : skull.faces) {
            if (f >= mesh.num_faces()) {
                report("face id out of range");
                break;
            }
            reconstructed += mesh.face_area[f];
        }
        if (reconstructed != skull.area) report("O1: DP area != reconstructed area");

        // O3 — a region made of mesh triangles is inside P by construction, so the
        // content of this check is that the region is non-empty and well formed.
        if (skull.faces.empty()) report("O3: region contains no triangles");

        // O2 — simple, correctly oriented, and convex. Collinear vertices are legal:
        // the paper's angle condition admits a straight continuation, so convexity
        // must be tested as "no right turn", not "every turn is a left turn".
        const Polygon &b = skull.boundary;
        if (b.size() < 3) {
            report("O2: boundary has fewer than 3 vertices");
        } else {
            if (!b.is_simple()) report("O2: boundary is not simple");
            if (!b.is_counterclockwise_oriented()) report("O2: boundary is not CCW");
            for (std::size_t i = 0; i < b.size(); ++i) {
                const Point &p = b[i];
                const Point &c = b[(i + 1) % b.size()];
                const Point &n = b[(i + 2) % b.size()];
                if (CGAL::orientation(p, c, n) == CGAL::RIGHT_TURN) {
                    report("O2: boundary is not convex at vertex "
                           + std::to_string((i + 1) % b.size()));
                    break;
                }
            }
            // The boundary encloses exactly the reported triangles.
            if (CGAL::abs(b.area()) != skull.area) report("O1: boundary area != DP area");
        }

        if (count == 0) return {};
        return {false, problems.str()};
    }

    MeshedSkull largest_convex_subpolygon(const Mesh &mesh) {
        return sweep::solve(mesh);
    }

    Polygon convex_skull(const Polygon &polygon) {
        // A convex polygon is its own convex skull, exactly: the largest convex
        // subset of a convex set is the set itself. O(n), and it skips building
        // the mesh and running the DP entirely.
        //
        // This is not merely an optimisation here. Aronov's answer is restricted to
        // unions of mesh triangles, so on a convex polygon whose triangulation has
        // an interior structure the DP can return LESS than the whole polygon; the
        // shortcut makes the answer exact on convex input rather than mesh-limited.
        if (polygon.is_simple() && polygon.is_convex()) {
            Polygon out = polygon;
            if (out.orientation() == CGAL::CLOCKWISE) out.reverse_orientation();
            return out;
        }

        const Mesh mesh = build_mesh(polygon);
        // The shear is area-preserving (validate checks det == 1), so the mesh's total
        // area must equal the cleaned polygon's exactly.
        if (const MeshCheck check = validate(mesh, CGAL::abs(pre_checks(polygon).area()));
            !check.ok)
            throw std::runtime_error("aronov: mesh validation failed: " + check.message);
        const MeshedSkull skull = largest_convex_subpolygon(mesh);
        if (const SkullCheck check = check_skull(mesh, skull); !check.ok)
            throw std::runtime_error("aronov: skull validation failed: " + check.message);
        return skull.boundary;
    }

}  // namespace potato_peeling::aronov

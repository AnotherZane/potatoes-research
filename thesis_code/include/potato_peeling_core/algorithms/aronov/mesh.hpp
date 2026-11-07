#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "potato_peeling_core/types.hpp"

namespace potato_peeling::aronov {

    // "No incident triangle on this side", i.e. the exterior of P — the unbounded
    // face or a hole. The mesh M covers exactly int(P) (Aronov p352), so every side
    // of an edge that is not a stored triangle IS the exterior. Triangles are closed
    // and the exterior is open (p353), so touching an exterior side is legal.
    inline constexpr std::size_t NO_FACE = static_cast<std::size_t>(-1);

    // Indexed, sweep-normalised triangular mesh.
    //
    // Everything is stored in SHEARED coordinates: `shear` is an area-preserving
    // (det == 1) affine map chosen so that no two vertices share an x-coordinate,
    // which is the paper's standing general-position assumption (p353). Because the
    // shear is affine, convexity, orientation and the mesh combinatorics are
    // unchanged, and because det == 1 all areas are exact in either frame — only
    // the reported boundary needs `unshear`.
    struct Mesh {
        // Vertices, sorted by x. All x-coordinates are distinct, hence no edge of
        // the mesh is vertical and the sweep order is a plain x-comparison.
        std::vector<Point> vertices;

        // Triangles of M, CCW-oriented, as vertex ids.
        std::vector<std::array<std::size_t, 3>> faces;
        std::vector<FT> face_area;

        // Edges of M as {left, right} vertex ids, with x(left) < x(right).
        std::vector<std::array<std::size_t, 2>> edges;

        // {above, below} face ids relative to the directed segment left -> right.
        // NO_FACE means the exterior of P lies on that side.
        std::vector<std::array<std::size_t, 2>> edge_faces;

        // Incident edges per vertex, CCW by the angle of the direction v -> other,
        // measured from the +x axis.
        std::vector<std::vector<std::size_t>> vertex_edges;

        // E_L(v) and E_R(v) (p356): the incident edges whose other endpoint lies to
        // the left / to the right of v. Both are ordered TOPMOST FIRST, the order
        // the paper's list L(e) is built in ("the first entry corresponds to the
        // topmost edge (smallest slope)", p357).
        std::vector<std::vector<std::size_t>> vertex_left_edges;
        std::vector<std::vector<std::size_t>> vertex_right_edges;

        // opp() (p356), precomputed per edge e against v = left(e), so both are
        // indexed by edge id and address E_L(v) = vertex_left_edges[left(e)].
        //
        // The convex continuations of e are a contiguous run of E_L(v) in
        // topmost-first order (`bottom_valid_prefix` / `top_valid_suffix` in
        // predicates.hpp), which is what makes the paper's prefix-maximum lookup
        // (Observation 1, p357) legitimate.
        //
        // bottom: e is the OUTGOING edge of the bottom chain at v; the valid
        //         continuations are E_L(v)[0 .. bottom_prefix_len[e]).
        std::vector<std::size_t> bottom_prefix_len;
        // top:    e is the INCOMING edge of the top chain at v; the valid
        //         continuations are E_L(v)[top_suffix_begin[e] .. end).
        std::vector<std::size_t> top_suffix_begin;

        Aff_transformation shear{1, 0, 0, 0, 1, 0};
        Aff_transformation unshear{1, 0, 0, 0, 1, 0};

        [[nodiscard]] std::size_t num_vertices() const { return vertices.size(); }
        [[nodiscard]] std::size_t num_faces() const { return faces.size(); }
        [[nodiscard]] std::size_t num_edges() const { return edges.size(); }

        [[nodiscard]] const Point &point(const std::size_t v) const { return vertices[v]; }

        [[nodiscard]] std::size_t left(const std::size_t e) const { return edges[e][0]; }
        [[nodiscard]] std::size_t right(const std::size_t e) const { return edges[e][1]; }

        [[nodiscard]] const Point &left_point(const std::size_t e) const {
            return vertices[edges[e][0]];
        }
        [[nodiscard]] const Point &right_point(const std::size_t e) const {
            return vertices[edges[e][1]];
        }

        [[nodiscard]] std::size_t face_above(const std::size_t e) const { return edge_faces[e][0]; }
        [[nodiscard]] std::size_t face_below(const std::size_t e) const { return edge_faces[e][1]; }

        // The endpoint of e that is not v. Precondition: v is an endpoint of e.
        [[nodiscard]] std::size_t other(std::size_t e, std::size_t v) const;

        // Edge id joining u and v, or NO_FACE if they are not adjacent. O(deg(u)).
        [[nodiscard]] std::size_t find_edge(std::size_t u, std::size_t v) const;

        // Total area of M, i.e. |area(P)|.
        [[nodiscard]] FT total_area() const;

        // Map a sheared point back to the caller's frame.
        [[nodiscard]] Point to_original(const Point &p) const { return unshear(p); }
    };

    // The shear parameter s of (x, y) -> (x + s*y, y) that separates all
    // x-coordinates. Returns 0 when the input already has distinct x-coordinates.
    // Throws std::invalid_argument on a duplicated point.
    FT separating_shear(const std::vector<Point> &points);

    // Build M as the constrained Delaunay triangulation of a simple polygon. No
    // Steiner points are added: the mesh vertices are exactly the polygon vertices
    // (after collinear-vertex removal).
    Mesh build_mesh(const Polygon &outer);

    // Mesh validation. `expected_area`, when given, is checked against total_area()
    // exactly.
    struct MeshCheck {
        bool ok{true};
        std::string message;

        explicit operator bool() const { return ok; }
    };

    MeshCheck validate(const Mesh &mesh, const std::optional<FT> &expected_area = std::nullopt);
}  // namespace potato_peeling::aronov

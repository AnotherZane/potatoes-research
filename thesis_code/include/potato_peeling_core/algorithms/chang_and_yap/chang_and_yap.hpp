#pragma once

#include <cassert>
#include <cstddef>
#include <optional>
#include <vector>
#include <functional>
#include <tuple>
#include <utility>

#include "boost/unordered_map.hpp"
#include "potato_peeling_core/types.hpp"

namespace potato_peeling::chang_and_yap {

    struct Chord {
        // Endpoints of the maximal chord in P (may extend past reflex vertices).
        Point start;
        Point end;

        // Indices of the defining corners v_i, v_j.
        std::size_t i{}, j{};

        Chord() = default;
        Chord(Point s, Point e, const std::size_t i_, const std::size_t j_)
            : start(std::move(s)), end(std::move(e)), i(i_), j(j_) {}

        bool operator==(const Chord &other) const { return (i == other.i && j == other.j); }
    };

    struct ChainEntry {
        FT area;

        Point x0;
        Point xm;

        // Pivots as polygon corner indices; every pivot is a corner of P. Indices
        // rather than Points keep entries small, and give the neighbours needed
        // to determine C+ (Section 2).
        std::vector<std::size_t> pivot_indices;

        // Indices of the bounding chords C and C', so L0 and L1 are available to
        // the composition checks. The paper does not say how to store them.
        std::size_t chord_c_idx;
        std::size_t chord_c_prime_idx;
    };

    // An entry's pivot points, from its indices.
    inline std::vector<Point> pivots_of(const std::vector<std::size_t> &idx, const Polygon &poly) {
        std::vector<Point> out;
        out.reserve(idx.size());
        for (const std::size_t i : idx) out.push_back(poly[i]);
        return out;
    }

    // What a Lemma 8 composition-compatibility check (p173) needs about one
    // (L, pivots, R) path. Both members depend only on the cache key: the key fixes
    // the chords, hence the pivots and L's origin, and each A-entry is unique per
    // pair (Theorem 5) and immutable once set.
    //
    // `last_offset` is A[m], the only lever offset equivalent_given_params reads,
    // so a cache hit costs Theta(log n) and allocates nothing.
    struct CritEntry {
        std::vector<FT> params;    // sorted critical parameters along L
        Vector last_offset{0, 0};  // A[m] of the (L, p_1..p_m)-path
    };

    // Memoised CritEntry per path, shared across the whole A-matrix build. Hashed:
    // nothing needs key order. Node-based, so references from set_* survive rehash.
    struct CritParamCache {
        boost::unordered_map<std::tuple<std::size_t, std::size_t, int>, CritEntry> l0;
        boost::unordered_map<std::pair<std::size_t, std::size_t>, CritEntry> l1;

        // A pointer into the map, or nullptr if absent.
        CritEntry *get_l0(const std::size_t left_idx, const std::size_t mid_idx,
                          const int variant) {
            const auto key = std::make_tuple(left_idx, mid_idx, variant);
            if (const auto it = l0.find(key); it != l0.end()) return &it->second;
            return nullptr;
        }

        CritEntry &set_l0(const std::size_t left_idx, const std::size_t mid_idx, const int variant,
                          CritEntry val) {
            return l0.emplace(std::make_tuple(left_idx, mid_idx, variant), std::move(val))
                .first->second;
        }

        CritEntry *get_l1(const std::size_t mid_idx, const std::size_t right_idx) {
            const auto key = std::make_pair(mid_idx, right_idx);
            if (const auto it = l1.find(key); it != l1.end()) return &it->second;
            return nullptr;
        }

        CritEntry &set_l1(const std::size_t mid_idx, const std::size_t right_idx, CritEntry val) {
            return l1.emplace(std::make_pair(mid_idx, right_idx), std::move(val)).first->second;
        }
    };

    // Collapse the lazy expression DAG of a stored value. A Lazy_exact_nt keeps its
    // whole unevaluated DAG, which dominates a solve's memory; CGAL::exact() forces
    // the value and releases the DAG. Values are unchanged; this trades time for
    // memory.
    namespace detail {
        inline void prune(const FT &v) {
            (void)CGAL::exact(v);
        }

        inline void prune(const Point &p) {
            (void)CGAL::exact(p.x());
            (void)CGAL::exact(p.y());
        }

        inline void prune(const Vector &v) {
            (void)CGAL::exact(v.x());
            (void)CGAL::exact(v.y());
        }
    }  // namespace detail

    using ChordPair = std::pair<std::size_t, std::size_t>;

    struct AMatrix {
        using Map = boost::unordered_map<ChordPair, ChainEntry>;

        Map data;

        void reserve(std::size_t n) { data.reserve(n); }

        [[nodiscard]] bool is_valid(std::size_t a, std::size_t b) const {
            return data.count({a, b}) > 0;
        }

        [[nodiscard]] ChainEntry get(std::size_t a, std::size_t b) const {
            const auto it = data.find({a, b});
            // Callers must guard with is_valid(): an absent entry is the paper's
            // -infinity, which FT cannot represent.
            assert(it != data.end() && "AMatrix::get on absent key — guard with is_valid()");
            return it->second;
        }

        // Like get() without the copy. References survive rehash (node-based map),
        // and an A-entry is immutable once set (Theorem 5).
        [[nodiscard]] const ChainEntry &get_ref(std::size_t a, std::size_t b) const {
            const auto it = data.find({a, b});
            assert(it != data.end() && "AMatrix::get_ref on absent key — guard with is_valid()");
            return it->second;
        }

        void set(std::size_t a, std::size_t b, ChainEntry entry) {
            // The entry is final, so collapse its DAGs.
            detail::prune(entry.area);
            detail::prune(entry.x0);
            detail::prune(entry.xm);
            data[{a, b}] = std::move(entry);
        }
    };

    struct MEntry {
        FT area;
        // The middle chord that achieved the maximum in the recurrence; nullopt
        // if the entry came directly from A(C,C').
        std::optional<std::size_t> best_c_mid_idx;
        // Chain-definition state (Section 3): the entry's chain endpoints and
        // the traversal directions leaving x_first / arriving at x_last. Used
        // by the seam-compatibility and convex-core guards in Step II/III.
        Point x_first{0, 0};
        Point x_last{0, 0};
        Vector d_first{0, 0};
        Vector d_last{0, 0};
    };

    struct MMatrix {
        using Map = boost::unordered_map<ChordPair, MEntry>;

        Map data;

        void reserve(std::size_t n) { data.reserve(n); }

        [[nodiscard]] bool is_valid(std::size_t a, std::size_t b) const {
            return data.count({a, b}) > 0;
        }

        [[nodiscard]] MEntry get(std::size_t a, std::size_t b) const {
            const auto it = data.find({a, b});
            // Callers must guard with is_valid(); see AMatrix::get.
            assert(it != data.end() && "MMatrix::get on absent key — guard with is_valid()");
            return it->second;
        }

        // Store the value if the entry is absent or this is a strict improvement.
        // Returns whether it was stored.
        bool set(std::size_t a, std::size_t b, const FT &val,
                 const std::optional<std::size_t> c_mid_idx) {
            auto it = data.find({a, b});
            if (it == data.end() || val > it->second.area) {
                // Prune only what is kept.
                detail::prune(val);
                data[{a, b}] = {val, c_mid_idx};
                return true;
            }
            return false;
        }

        // Full-entry overload carrying the chain-definition state.
        bool set(std::size_t a, std::size_t b, MEntry entry) {
            auto it = data.find({a, b});
            if (it == data.end() || entry.area > it->second.area) {
                detail::prune(entry.area);
                detail::prune(entry.x_first);
                detail::prune(entry.x_last);
                detail::prune(entry.d_first);
                detail::prune(entry.d_last);
                data[{a, b}] = std::move(entry);
                return true;
            }
            return false;
        }
    };

    enum class CompositionCase { CASE_A, CASE_B };

    struct VerificationResult {
        bool valid;
        FT composed_area;
        // Which case verified; it decides the composed pivot sequence (Section 5).
        CompositionCase case_type;
        // Case A's junction pivot depends on the parity of m'' (Section 5).
        Point junction_pivot_a;  // for case A: p or q
    };

    // Step III (Section 6): max{M(C,C') + M(C',C) : C ∈ Ξ_i, C' ∈ Ξ_j}.
    struct PotatoResult {
        FT area;                  // maximum area found
        std::size_t c_idx;        // index of first bounding chord
        std::size_t c_prime_idx;  // index of second bounding chord
        // The chord indices allow reconstructing the polygon, for which the paper
        // gives no procedure ("easy to modify to compute potato itself").
    };

    bool is_contained(const Polygon &polygon, const Segment &segment);

    bool extends_past_corner(const Polygon &polygon, std::size_t k, const Vector &d);

    bool on_chord_segment(const Point &x, const Chord &c);

    Point reference_point(const Chord &C, const Polygon &poly);

    FT triangle_area_signed(const Point &a, const Point &b, const Point &c);

    std::pair<Point, Point> directed_C_plus(const Chord &c, const Polygon &P,
                                            const std::vector<bool> &reflex_flags);

    std::vector<Vector> compute_lever_offsets(const std::vector<Point> &pivots,
                                              const Point &origin);

    std::optional<Point> solve_balanced_chain_x0(const Line &L0, const Line &L1,
                                                 const std::vector<Point> &pivots);

    std::optional<FT> compute_critical_param_type_a(const Point &r,
                                                    std::size_t i,  // strip/node index (1-based)
                                                    const std::vector<Vector> &A,
                                                    const Point &origin, const Vector &dir,
                                                    int sign_i);

    std::optional<FT> compute_critical_param_type_b(std::size_t k, const std::vector<Vector> &A,
                                                    const Point &origin, const Vector &dir,
                                                    int sign_k);

    std::vector<FT> compute_critical_params(const Line &L0, const std::vector<Point> &pivots,
                                            const std::vector<std::size_t> &R, const Polygon &poly);

    // A[m], the only lever offset the equivalence test reads. O(m), paid once per
    // cache key.
    Vector last_lever_offset(const std::vector<Point> &pivots, const Point &origin);

    // The equivalence test given the memoised path data. Takes the pivot count
    // (only its parity is read) and the cached A[m]: Theta(log n), no allocation.
    bool equivalent_given_params(const Line &L, std::size_t m, const Vector &last_offset,
                                 const std::vector<FT> &params, const Point &x, const Point &y);

    // Get or compute the memoised path data for the two composition checks;
    // `build` runs only on a miss.
    const CritEntry &cached_l0_entry(CritParamCache &cache, std::size_t left_idx,
                                     std::size_t mid_idx, int variant, const Line &L,
                                     const std::function<CritEntry()> &build);

    const CritEntry &cached_l1_entry(CritParamCache &cache, std::size_t mid_idx,
                                     std::size_t right_idx, const Line &L,
                                     const std::function<CritEntry()> &build);

    std::optional<std::pair<Point, Point>>
    find_balanced_chord(const Point &p,  // reflex corner — pivot
                        const Line &L0,  // line through bounding chord C
                        const Line &L1);

    FT length_1_chain_area(const Point &c_ref,  // reference point of C
                           const Point &x0, const Point &x1, const Point &c_prime_ref);

    bool ordered_along_line(const Point &a, const Point &b, const Point &c, const Point &d);

    bool ordered_along_line_5(const Point &a, const Point &b, const Point &c, const Point &d,
                              const Point &e);

    bool
    chains_do_not_intersect(const ChainEntry &left, const ChainEntry &right,
                            const Polygon &poly);  // (C_mid, C')-chain — L2/poly implicit via nodes

    std::optional<Point>
    compute_junction_node_A(const Point &x_prime_next,    // x'_{m'+1} = x'_{m'} | p_{m'+1}
                            const Point &x_double_first,  // x''_0
                            const Line &L0, const Line &L1);

    std::pair<Point, Point> compute_junction_nodes_B(const Point &x_prime_last,  // x'_{m'}
                                                     const Point &pivot_p,       // p
                                                     const Point &pivot_q);

    // Lemma 8 verification of a composition, sharing one CritParamCache across
    // the A-matrix build. left_idx/mid_idx/right_idx are the Xi indices of
    // C_left/C_mid/C_right.
    VerificationResult verify_case_A_cached(const ChainEntry &left, const ChainEntry &right,
                                            const Chord &C_mid, const Chord &C_left,
                                            const Chord &C_right, const Polygon &poly,
                                            CritParamCache &pcache, std::size_t left_idx,
                                            std::size_t mid_idx, std::size_t right_idx);

    VerificationResult verify_case_B_cached(const ChainEntry &left, const ChainEntry &right,
                                            const Chord &C_mid, const Chord &C_left,
                                            const Chord &C_right, const Polygon &poly,
                                            CritParamCache &pcache, std::size_t left_idx,
                                            std::size_t mid_idx, std::size_t right_idx);

    VerificationResult verify_decomposition_cached(const ChainEntry &left, const ChainEntry &right,
                                                   const Chord &C_mid, const Chord &C_left,
                                                   const Chord &C_right, const Polygon &poly,
                                                   CritParamCache &pcache, std::size_t left_idx,
                                                   std::size_t mid_idx, std::size_t right_idx);

    std::vector<Point> reconstruct_chain_nodes(const ChainEntry &chain, const Polygon &poly);

    Polygon clip_by_halfplane(const Polygon &poly, const Point &A, const Point &B);

    void orient_chord_for_clipping(Point &p1, Point &p2, std::size_t reflex_idx,
                                   const Polygon &poly);

    bool orient_extremal_chord(Point &p1, Point &p2, const Chord &chord, const Polygon &poly,
                               const std::vector<bool> &reflex_flags);

    void collect_chords(const MMatrix &M, const AMatrix &A, const std::vector<Chord> &Xi,
                        const Polygon &poly, const std::vector<bool> &reflex_flags, std::size_t a,
                        std::size_t b, std::vector<std::pair<Point, Point>> &chords);

    bool in_C_plus(const Point &r, const Chord &c, const Polygon &P,
                   const std::vector<bool> &reflex_flags);

    bool R_intersect_halfplanes_empty(const std::vector<std::size_t> &R, const Chord &c,
                                      const Chord &c_prime, const Polygon &P,
                                      const std::vector<bool> &reflex_flags);

    Polygon pre_checks(const Polygon &input);

    std::vector<Chord> generate_extremal_chords(const Polygon &polygon,
                                                const std::vector<bool> &reflex_flags);

    bool is_admissible_pair(const Chord &c, const Chord &c_prime, const std::size_t &n);

    std::vector<std::size_t> compute_R(const Chord &c, const Chord &c_prime, const std::size_t &n);

    std::optional<ChainEntry> try_length_0_chain(const std::vector<std::size_t> &R,
                                                 const std::vector<Chord> &Xi, std::size_t a,
                                                 std::size_t b, const Polygon &poly,
                                                 const std::vector<bool> &reflex_flags);

    std::optional<ChainEntry> try_length_1_chain(const std::vector<std::size_t> &R,
                                                 const std::vector<Chord> &Xi, std::size_t a,
                                                 std::size_t b, const Polygon &poly,
                                                 const std::vector<bool> &reflex_flags);

    void compute_A_matrix_all_stages(AMatrix &A, const std::vector<Chord> &Xi, const Polygon &poly,
                                     const std::vector<bool> &reflex_flags, std::size_t n_vertices);

    Vector chord_dir(const Chord &c, const Polygon &poly);

    bool left_turn_ok(const Vector &u, const Vector &v);

    bool in_chain_wedge(const Point &x, const Chord &ca, const Chord &cb, const Polygon &poly);

    bool closed_core_ok(const Point &x_first, const Vector &d_first, const Point &x_last,
                        const Vector &d_last);

    FT chord_param(const Point &x, const Chord &c, const Polygon &poly);

    bool seam_ok(const MEntry &left, const Chord &mid, const MEntry &right, const Polygon &poly);

    std::optional<MEntry> m_leaf_from_A(const ChainEntry &ae, const Chord &ca, const Chord &cb,
                                        const Polygon &poly);

    void compute_M_stage_1(MMatrix &M, const std::vector<Chord> &Xi,
                           const std::vector<std::vector<std::size_t>> &Xi_by_vertex,
                           const Polygon &poly, const std::vector<bool> &reflex_flags);

    void compute_M_stage_s(MMatrix &M, const AMatrix &A, const std::vector<Chord> &Xi,
                           const std::vector<std::vector<std::size_t>> &Xi_by_vertex,
                           const Polygon &poly, std::size_t s);

    PotatoResult compute_potato_area(const MMatrix &M, const AMatrix &A,
                                     const std::vector<Chord> &Xi, const Polygon &poly);

    Polygon reconstruct_potato(const MMatrix &M, const AMatrix &A, const std::vector<Chord> &Xi,
                               const Polygon &poly, const std::vector<bool> &reflex_flags,
                               std::size_t a, std::size_t b);

    // The skull, with |Xi| and the number of A and M entries. The counts are
    // unset for a convex P, where the DP is skipped.
    struct Solution {
        Polygon skull;
        std::optional<std::size_t> n_xi, a_entries, m_entries;
    };

    Solution convex_skull(const Polygon &polygon);
}  // namespace potato_peeling::chang_and_yap

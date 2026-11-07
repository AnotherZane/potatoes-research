#pragma once

// Refinement: local ascent from any convex Q inside P (for example a Cabello or
// Aronov answer) until no local move improves its area.
//
// Chang & Yap 1984, Lemma 4: every edge of a maximal convex Q is one of three
// kinds:
//
//     (alpha) part of an edge of P;
//     (beta)  part of a chord whose MIDPOINT is a reflex corner of P;
//     (gamma) part of an extremal chord, through two corners of P.
//
// An edge touching nothing is none of these, so it can be PUSHED out. An edge
// touching one reflex corner off-centre violates (beta), so it can be ROTATED
// until the corner is its midpoint (the Butterfly Lemma). SPLIT and ADJOIN add
// edges and vertices (see Options). A fixed point satisfies Lemma 4 on every
// edge it can reach.
//
// Every accepted state is certified convex and contained by cabello::certify in
// EPECK, so the result's area is a lower bound on A*(P). It is a local optimum
// only.

#include <cstddef>
#include <string>

#include "potato_peeling_core/types.hpp"

namespace potato_peeling::refinement {

    struct Options {
        // Sweeps over all edges. The ascent usually converges in 2-6.
        unsigned max_sweeps{12};
        // Rotation search: rot_steps angles spanning +/- rot_span radians about
        // the edge's current normal, each followed by a fresh push.
        unsigned rot_steps{9};
        double rot_span{0.30};
        // Bisection depth for the push (in double; see push()).
        unsigned push_rounds{40};
        // SPLIT: replace one edge with two, at +/- split_span radians about its
        // normal, each pushed independently; split_steps trial angles, 0 disables
        // it. Push and rotate keep the number of edges fixed, so without split an
        // optimum with more edges than the start is unreachable. One edge pinned
        // by a reflex corner is blocked along its whole length; two edges either
        // side of the corner can each advance past it. Off by default because it
        // costs 20-100x the time.
        unsigned split_steps{0};
        double split_span{0.25};
        // Bisection depth for a split trial. Coarse, since the next sweep's push
        // refines the offset.
        unsigned split_push_rounds{12};
        // ADJOIN: absorb a corner of P into Q, replacing Q with hull(Q u {v}).
        // Rounds of sweeping ALL of P's corners; 0 disables it.
        //
        // Goodman 1981, Proposition 2: if X in P sees every point of a convex Q,
        // then hull(Q u {X}) is the union of the segments Xq, a convex subset of
        // P, larger when X is outside Q. Adjoin puts vertices on P's corners
        // (Lemma 4 type gamma); split works at free angles (type beta).
        unsigned adjoin_rounds{0};
        // Cap on the half-plane count, as a multiple of the start's.
        unsigned max_growth{4};
    };

    struct Result {
        // Area of the refined polygon. A genuine lower bound on A*(P).
        FT area{0};
        // The refined polygon, CCW. Convex and contained, certified in EPECK.
        Polygon polygon;
        // Area of the polygon handed in, for the gain.
        FT start_area{0};
        unsigned sweeps{0};
        // certify calls, the dominant cost.
        std::size_t evaluations{0};
        // Half-plane counts at the start and end.
        std::size_t halfplanes_start{0};
        std::size_t halfplanes_end{0};
        unsigned splits_accepted{0};
        unsigned adjoins_accepted{0};
        bool ok{false};
        // Non-empty iff ok is false.
        std::string message;
    };

    // `start` must be a convex subset of `p`; if it is not, the result carries
    // ok = false and a message rather than throwing. Orientation is normalised
    // internally.
    Result refine(const Polygon &p, const Polygon &start, const Options &opts = {});

}  // namespace potato_peeling::refinement

#pragma once

#include <vector>

#include "potato_peeling_core/algorithms/aronov/aronov.hpp"
#include "potato_peeling_core/algorithms/aronov/mesh.hpp"

namespace potato_peeling::aronov::sweep {

    // The paper's O(n^2) algorithm.
    //
    // §2.1 (p355-356): one left-to-right sweep. Every edge crossing the sweep-line
    // carries Ea/Eb lists of the partners it is still potentially promising with,
    // plus an inverted list EI so the pairs owned by other edges can be unlinked in
    // O(1). Pairs are born when the later edge begins, die when an exterior gap
    // opens between the two, and are reported when the earlier edge ends.
    //
    // Returned sorted by min X(e, e').
    std::vector<EdgePair> promising_pairs(const Mesh &mesh);

    // §2.2 (p356-357): the Q-values, via opp() and prefix maxima instead of scanning
    // every candidate. Uses Mesh::bottom_prefix_len / top_suffix_begin, whose
    // run-contiguity is what makes Observation 1's lookup legitimate.
    QTable q_table(const Mesh &mesh, const std::vector<EdgePair> &promising);

    // promising_pairs, then q_table, then the best pair reconstructed.
    MeshedSkull solve(const Mesh &mesh);
}  // namespace potato_peeling::aronov::sweep

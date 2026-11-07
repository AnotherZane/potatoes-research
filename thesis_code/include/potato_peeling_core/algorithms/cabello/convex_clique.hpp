#pragma once

// phi(H, s): a maximum-area convex clique of H with s as its highest point.
// Lemma 8 [C17 p9-p10].
//
// The paper's O(|V|^2) relies on the valid predecessors of x_i forming a
// contiguous arc in the circular order of edges at x_i [C17 p10], and gets those
// orders from a dual line arrangement. We instead sort the arrivals at each x_i
// and answer each arc with two binary searches and a sparse-table range maximum:
// O(|V|^2 log|V|).

#include <cstddef>

#include "potato_peeling_core/algorithms/cabello/graph.hpp"

namespace potato_peeling::cabello {

    Clique phi(const PointGraph &h, std::size_t s);

}  // namespace potato_peeling::cabello

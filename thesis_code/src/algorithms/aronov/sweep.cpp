#include "potato_peeling_core/algorithms/aronov/sweep.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>

#include "potato_peeling_core/algorithms/aronov/predicates.hpp"
#include "potato_peeling_core/types.hpp"

namespace potato_peeling::aronov::sweep {
    namespace {
        constexpr std::size_t NIL = static_cast<std::size_t>(-1);

        // Collapse a lazy-exact value to its exact representation before storing it.
        // FT is an expression DAG and Q(e,e') = area(T(e,e')) + Q(sub), so otherwise
        // every stored entry keeps the whole arithmetic history it was built from.
        FT settled(FT v) {
            CGAL::exact(v);
            return v;
        }

        // ── The sweep-line status ─────────────────────────────────────────────
        //
        // A plain doubly-linked list, bottom to top, of the edges crossing the
        // sweep-line. No balanced tree is needed: every edge incident to a vertex v
        // occupies a CONTIGUOUS run of the status just left of v (any edge strictly
        // between two of them at x -> v.x would have to pass through v, hence be
        // incident to v as well), so the edges ending at v can be spliced out and
        // the edges starting at v spliced into the same slot in O(deg(v)).
        //
        // The gap between two consecutive status edges is exactly one cell of the
        // arrangement — a triangle or the exterior — so the gap above `e` is interior
        // iff face_above(e) != NO_FACE. That is a STATIC property of the edge, which
        // is what makes "walk until you hit a part of the boundary of P" (p355) an
        // O(1) test per step with no interval bookkeeping at all.
        struct Status {
            std::vector<std::size_t> up, down;
            std::size_t bottom = NIL;

            explicit Status(const std::size_t n) : up(n, NIL), down(n, NIL) {}

            void insert_above(const std::size_t anchor, const std::size_t e) {
                if (anchor == NIL) {
                    up[e] = bottom;
                    down[e] = NIL;
                    if (bottom != NIL) down[bottom] = e;
                    bottom = e;
                } else {
                    const std::size_t nxt = up[anchor];
                    up[anchor] = e;
                    down[e] = anchor;
                    up[e] = nxt;
                    if (nxt != NIL) down[nxt] = e;
                }
            }

            void erase(const std::size_t e) {
                const std::size_t p = down[e], n = up[e];
                if (p != NIL)
                    up[p] = n;
                else
                    bottom = n;
                if (n != NIL) down[n] = p;
                up[e] = down[e] = NIL;
            }
        };

        // ── Potentially-promising pair lists ──────────────────────────────────
        //
        // Ea(e) holds the partners ABOVE e, ordered upward, for the pairs e owns;
        // Eb(e) the partners BELOW e, ordered downward. Both are filled in one walk
        // when e begins, so insertion order IS the order along the sweep-line, and
        // removals from the middle preserve it. EI(e) is the inverted list of the
        // pairs where e is the partner rather than the owner, so "The End" can
        // unlink them in O(1) each (p355).
        struct PairLists {
            struct Node {
                std::size_t upper, lower;
                std::size_t owner, partner;
                bool in_above{};
                std::size_t prev = NIL, next = NIL;    // within the owner's list
                std::size_t iprev = NIL, inext = NIL;  // within EI(partner)
                bool alive = true;
            };

            std::vector<Node> arena;
            std::vector<std::size_t> ea_head, ea_tail, eb_head, eb_tail, ei_head;

            explicit PairLists(const std::size_t n)
                : ea_head(n, NIL), ea_tail(n, NIL), eb_head(n, NIL), eb_tail(n, NIL),
                  ei_head(n, NIL) {}

            void add(const std::size_t upper, const std::size_t lower, const std::size_t owner,
                     const bool in_above) {
                const std::size_t id = arena.size();
                Node node;
                node.upper = upper;
                node.lower = lower;
                node.owner = owner;
                node.partner = in_above ? upper : lower;
                node.in_above = in_above;

                auto &tail = in_above ? ea_tail[owner] : eb_tail[owner];
                auto &head = in_above ? ea_head[owner] : eb_head[owner];
                node.prev = tail;
                if (tail != NIL)
                    arena[tail].next = id;
                else
                    head = id;
                tail = id;

                node.inext = ei_head[node.partner];
                if (node.inext != NIL) arena[node.inext].iprev = id;
                ei_head[node.partner] = id;

                arena.push_back(node);
            }

            void unlink(const std::size_t id) {
                Node &node = arena[id];
                if (!node.alive) return;
                node.alive = false;

                auto &tail = node.in_above ? ea_tail[node.owner] : eb_tail[node.owner];
                auto &head = node.in_above ? ea_head[node.owner] : eb_head[node.owner];
                if (node.prev != NIL)
                    arena[node.prev].next = node.next;
                else
                    head = node.next;
                if (node.next != NIL)
                    arena[node.next].prev = node.prev;
                else
                    tail = node.prev;

                if (node.iprev != NIL)
                    arena[node.iprev].inext = node.inext;
                else
                    ei_head[node.partner] = node.inext;
                if (node.inext != NIL) arena[node.inext].iprev = node.iprev;
            }
        };
    }  // namespace

    std::vector<EdgePair> promising_pairs(const Mesh &mesh) {
        const std::size_t ne = mesh.num_edges();
        Status status(ne);
        PairLists lists(ne);

        std::vector<EdgePair> reported;
        std::vector<bool> is_new(ne, false);

        // Stamp arrays for the exterior kill, so the reachable sets never need to be
        // cleared between events.
        std::vector<std::size_t> below_stamp(ne, 0), above_stamp(ne, 0);
        std::size_t stamp = 0;

        auto report = [&](const std::size_t upper, const std::size_t lower) {
            // A pair whose overlap is a single point has no
            // trapezoid, and a pair failing the p354 conditions can never carry a Q.
            if (!overlap(mesh, upper, lower)) return;
            if (!pair_precondition(mesh, upper, lower)) return;
            reported.emplace_back(upper, lower);
        };

        // Kill every surviving pair that straddles the exterior gap between `lo` and
        // the edge above it. Candidates reach only as far as the next exterior gap in
        // each direction: anything beyond one was already killed when THAT gap
        // formed, so it is no longer in any list.
        auto kill_across_gap = [&](const std::size_t lo, const std::size_t hi) {
            ++stamp;
            std::vector<std::size_t> below, above;

            for (std::size_t cur = lo; cur != NIL;) {
                below_stamp[cur] = stamp;
                below.push_back(cur);
                const std::size_t nxt = status.down[cur];
                if (nxt == NIL || mesh.face_above(nxt) == NO_FACE) break;
                cur = nxt;
            }
            for (std::size_t cur = hi; cur != NIL;) {
                above_stamp[cur] = stamp;
                above.push_back(cur);
                if (mesh.face_above(cur) == NO_FACE) break;
                cur = status.up[cur];
            }

            // Partners below the gap form a suffix of Eb (ordered downward), and
            // partners above it a suffix of Ea (ordered upward), so scanning from the
            // tail costs O(1) whenever there is nothing to discard (p356).
            for (const std::size_t e : above)
                while (lists.eb_tail[e] != NIL
                       && below_stamp[lists.arena[lists.eb_tail[e]].partner] == stamp)
                    lists.unlink(lists.eb_tail[e]);

            for (const std::size_t e : below)
                while (lists.ea_tail[e] != NIL
                       && above_stamp[lists.arena[lists.ea_tail[e]].partner] == stamp)
                    lists.unlink(lists.ea_tail[e]);
        };

        // Vertices are stored in x-order, which is the sweep order.
        for (std::size_t v = 0; v < mesh.num_vertices(); ++v) {
            const auto &el = mesh.vertex_left_edges[v];  // topmost first
            const auto &er = mesh.vertex_right_edges[v];

            // (1) The End (p355): every surviving pair using an edge that stops here
            // is now complete, since X(e, e') closes at the first of the two right
            // endpoints. Report and unlink, then splice the edge out of the status.
            std::size_t anchor = NIL;
            if (!el.empty()) {
                anchor = status.down[el.back()];  // el.back() is the lowest
                for (const std::size_t e : el) {
                    for (std::size_t id = lists.ea_head[e]; id != NIL;) {
                        const std::size_t nxt = lists.arena[id].next;
                        report(lists.arena[id].upper, lists.arena[id].lower);
                        lists.unlink(id);
                        id = nxt;
                    }
                    for (std::size_t id = lists.eb_head[e]; id != NIL;) {
                        const std::size_t nxt = lists.arena[id].next;
                        report(lists.arena[id].upper, lists.arena[id].lower);
                        lists.unlink(id);
                        id = nxt;
                    }
                    for (std::size_t id = lists.ei_head[e]; id != NIL;) {
                        const std::size_t nxt = lists.arena[id].inext;
                        report(lists.arena[id].upper, lists.arena[id].lower);
                        lists.unlink(id);
                        id = nxt;
                    }
                    status.erase(e);
                }
            } else {
                // No edge ends here, so the slot must be located. Nothing in the
                // status passes through v: an edge that did would be incident to v,
                // and E_L(v) is empty while E_R(v) is not yet inserted.
                const FT vx = mesh.point(v).x();
                const FT vy = mesh.point(v).y();
                for (std::size_t cur = status.bottom; cur != NIL; cur = status.up[cur]) {
                    if (edge_y_at(mesh, cur, vx) > vy) break;
                    anchor = cur;
                }
            }

            // (2) A New Beginning (p355). Insert every edge starting here, bottom to
            // top, into the vacated slot.
            std::size_t at = anchor;
            for (std::size_t i = er.size(); i-- > 0;) {
                status.insert_above(at, er[i]);
                at = er[i];
                is_new[er[i]] = true;
            }

            // Then walk out from each new edge, in the order they cross the
            // sweep-line just right of v. Walking up adds nothing for edges that are
            // also new here — the pair is added instead by that edge's downward walk —
            // so each new/new pair is created exactly once, while the walk still
            // passes through them.
            for (std::size_t i = er.size(); i-- > 0;) {
                const std::size_t e = er[i];

                for (std::size_t cur = e; mesh.face_above(cur) != NO_FACE;) {
                    const std::size_t nxt = status.up[cur];
                    if (nxt == NIL) break;
                    if (!is_new[nxt]) lists.add(nxt, e, e, true);
                    cur = nxt;
                }
                for (std::size_t cur = e;;) {
                    const std::size_t nxt = status.down[cur];
                    if (nxt == NIL || mesh.face_above(nxt) == NO_FACE) break;
                    lists.add(e, nxt, e, false);
                    cur = nxt;
                }
            }
            for (const std::size_t e : er) is_new[e] = false;

            // (3) Exterior Parts (p356). Only the adjacencies touched at v can have
            // changed, and those are exactly the ones whose lower edge is the anchor
            // or one of the edges just inserted. A gap that was already exterior has
            // no straddling pairs left, so re-running the scan on it costs O(1).
            auto check_gap = [&](const std::size_t lo) {
                if (lo == NIL || mesh.face_above(lo) != NO_FACE) return;
                const std::size_t hi = status.up[lo];
                if (hi != NIL) kill_across_gap(lo, hi);
            };
            check_gap(anchor);
            for (const std::size_t e : er) check_gap(e);
        }

        std::sort(reported.begin(), reported.end(), [&](const EdgePair &a, const EdgePair &b) {
            const FT a0 = mesh.left_point(a.first).x(), a1 = mesh.left_point(a.second).x();
            const FT b0 = mesh.left_point(b.first).x(), b1 = mesh.left_point(b.second).x();
            return ((a0 > a1) ? a0 : a1) < ((b0 > b1) ? b0 : b1);
        });
        return reported;
    }
    QTable q_table(const Mesh &mesh, const std::vector<EdgePair> &promising) {
        QTable q;
        q.reserve(promising.size());

        // Group the pairs by the vertex where their strip begins, i.e. the further
        // right of the two left endpoints. Vertex ids are x-sorted, so that is simply
        // the larger id, and processing vertices in id order visits every sub-pair
        // before the pair that needs it (p357).
        std::vector<std::vector<EdgePair>> by_vertex(mesh.num_vertices());
        for (const auto &pair : promising)
            by_vertex[std::max(mesh.left(pair.first), mesh.left(pair.second))].push_back(pair);

        // Running maximum over a slice of E_L(v), with the argmax kept for the
        // backpointer. `run` absent means -infinity.
        struct Prefix {
            std::vector<std::optional<FT>> best;
            std::vector<std::size_t> arg;
        };

        std::vector<EdgePair> case_a, case_b, base;
        boost::unordered_map<std::size_t, Prefix> tables;

        for (std::size_t v = 0; v < mesh.num_vertices(); ++v) {
            if (by_vertex[v].empty()) continue;
            const auto &el = mesh.vertex_left_edges[v];  // topmost first

            case_a.clear();
            case_b.clear();
            base.clear();
            for (const auto &pair : by_vertex[v]) {
                const bool top_here = mesh.left(pair.first) == v;
                const bool bot_here = mesh.left(pair.second) == v;
                if (top_here && bot_here)
                    base.push_back(pair);
                else if (bot_here)
                    case_a.push_back(pair);  // extend the bottom chain
                else
                    case_b.push_back(pair);  // extend the top chain
            }

            // Base case (p356): the chains meet at v, so the strip is everything.
            for (const auto &pair : base)
                q.emplace(pair, QEntry{settled(trapezoid_area(mesh, pair.first, pair.second)),
                                       std::nullopt});

            // Case A — e crosses the sweep-line above v, e' starts at v. The valid
            // continuations are a leading run of E_L(v), so a PREFIX maximum keyed on
            // the top edge answers every pair at once. This is Observation 1 (p357).
            tables.clear();
            for (const auto &pair : case_a) {
                auto [it, inserted] = tables.try_emplace(pair.first);
                if (!inserted) continue;
                Prefix &pre = it->second;
                pre.best.assign(el.size(), std::nullopt);
                pre.arg.assign(el.size(), 0);
                std::optional<FT> run;
                std::size_t run_arg = 0;
                for (std::size_t k = 0; k < el.size(); ++k) {
                    if (const auto sub = q.find({pair.first, el[k]});
                        sub != q.end() && (!run || sub->second.area > *run)) {
                        run = sub->second.area;
                        run_arg = k;
                    }
                    pre.best[k] = run;
                    pre.arg[k] = run_arg;
                }
            }
            for (const auto &pair : case_a) {
                const std::size_t len = mesh.bottom_prefix_len[pair.second];
                if (len == 0) continue;
                const Prefix &pre = tables.at(pair.first);
                if (!pre.best[len - 1]) continue;
                q.emplace(pair, QEntry{settled(trapezoid_area(mesh, pair.first, pair.second)
                                               + *pre.best[len - 1]),
                                       EdgePair{pair.first, el[pre.arg[len - 1]]}});
            }

            // Case B — the mirror: e starts at v, e' crosses below it, and the valid
            // continuations are a trailing run, so a SUFFIX maximum keyed on the
            // bottom edge is what opp() indexes into.
            tables.clear();
            for (const auto &pair : case_b) {
                auto [it, inserted] = tables.try_emplace(pair.second);
                if (!inserted) continue;
                Prefix &suf = it->second;
                suf.best.assign(el.size(), std::nullopt);
                suf.arg.assign(el.size(), 0);
                std::optional<FT> run;
                std::size_t run_arg = 0;
                for (std::size_t k = el.size(); k-- > 0;) {
                    if (const auto sub = q.find({el[k], pair.second});
                        sub != q.end() && (!run || sub->second.area > *run)) {
                        run = sub->second.area;
                        run_arg = k;
                    }
                    suf.best[k] = run;
                    suf.arg[k] = run_arg;
                }
            }
            for (const auto &pair : case_b) {
                const std::size_t begin = mesh.top_suffix_begin[pair.first];
                if (begin >= el.size()) continue;
                const Prefix &suf = tables.at(pair.second);
                if (!suf.best[begin]) continue;
                q.emplace(pair, QEntry{settled(trapezoid_area(mesh, pair.first, pair.second)
                                               + *suf.best[begin]),
                                       EdgePair{el[suf.arg[begin]], pair.second}});
            }
        }

        return q;
    }

    MeshedSkull solve(const Mesh &mesh) {
        const std::vector<EdgePair> promising = promising_pairs(mesh);
        const QTable q = q_table(mesh, promising);

        const std::optional<EdgePair> winner = best_pair(mesh, q);
        if (!winner) return {};
        return reconstruct(mesh, q, *winner);
    }
}  // namespace potato_peeling::aronov::sweep

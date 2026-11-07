#include "potato_peeling_core/refinement/refinement.hpp"
#include "potato_peeling_core/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#include "potato_peeling_core/algorithms/cabello/certify.hpp"

namespace potato_peeling::refinement {

    namespace {
        // a*x + b*y <= c, outward normal (a,b).
        struct HalfPlane {
            FT a, b, c;
        };

        std::vector<HalfPlane> halfplanes_of(const Polygon &q) {
            std::vector<HalfPlane> hp;
            const std::size_t n = q.size();
            hp.reserve(n);
            for (std::size_t i = 0; i < n; ++i) {
                const Point &p = q[i];
                const Point &r = q[(i + 1) % n];
                const FT ex = r.x() - p.x(), ey = r.y() - p.y();
                const FT a = ey, b = -ex;  // outward for CCW
                hp.push_back({a, b, a * p.x() + b * p.y()});
            }
            return hp;
        }

        // Clip a CONVEX path by one half-plane. Exact.
        std::vector<Point> clip_convex(const std::vector<Point> &pts, const HalfPlane &h) {
            std::vector<Point> out;
            out.reserve(pts.size() + 2);
            const std::size_t n = pts.size();
            for (std::size_t i = 0; i < n; ++i) {
                const Point &p = pts[i];
                const Point &q = pts[(i + 1) % n];
                const FT sp = h.a * p.x() + h.b * p.y() - h.c;
                const FT sq = h.a * q.x() + h.b * q.y() - h.c;
                if (sp <= FT(0)) out.push_back(p);
                if ((sp < FT(0) && sq > FT(0)) || (sq < FT(0) && sp > FT(0))) {
                    const FT t = sp / (sp - sq);
                    out.emplace_back(p.x() + t * (q.x() - p.x()), p.y() + t * (q.y() - p.y()));
                }
            }
            return out;
        }

        // Intersect all the half-planes, starting from a bounding box.
        std::vector<Point> realise(const std::vector<HalfPlane> &hp,
                                   const std::vector<Point> &box) {
            std::vector<Point> pts = box;
            for (const auto &h : hp) {
                pts = clip_convex(pts, h);
                if (pts.size() < 3) return {};
            }
            std::vector<Point> out;
            out.reserve(pts.size());
            for (const auto &p : pts)
                if (out.empty() || p != out.back()) out.push_back(p);
            if (out.size() > 1 && out.front() == out.back()) out.pop_back();
            return out;
        }

        // area(Q) if Q is a convex subset of P, else nullopt. Exact, via certify.
        class Evaluator {
        public:
            explicit Evaluator(const Polygon &p) : p_(p) {}

            std::optional<std::pair<FT, std::vector<Point>>>
            operator()(const std::vector<HalfPlane> &hp, const std::vector<Point> &box) {
                std::vector<Point> pts = realise(hp, box);
                if (pts.size() < 3) return std::nullopt;
                ++calls_;
                const cabello::Certificate cert = cabello::certify(p_, pts);
                if (!cert.convex || !cert.contained) return std::nullopt;
                return std::make_pair(cert.area, std::move(pts));
            }

            // Certify a point set directly, returning its hull: adjoin works with
            // vertices rather than half-planes.
            std::optional<std::pair<FT, Polygon>> points(const std::vector<Point> &cand) {
                ++calls_;
                const cabello::Certificate cert = cabello::certify(p_, cand);
                if (!cert.convex || !cert.contained) return std::nullopt;
                return std::make_pair(cert.area, cert.polygon);
            }

            [[nodiscard]] std::size_t calls() const { return calls_; }

        private:
            const Polygon &p_;
            std::size_t calls_{0};
        };

        struct PushResult {
            bool ok{false};
            FT c;
            FT area;
            std::vector<Point> pts;
        };

        // Largest feasible offset for half-plane i, by bisection.
        //
        // Q grows with c, so the feasible offsets form an interval and its end is
        // the optimum for this coordinate. The offset is bisected in double; only
        // the evaluation at it is exact, which keeps the lazy-exact DAGs shallow.
        PushResult push(Evaluator &ev, std::vector<HalfPlane> hp, const std::size_t i,
                        const std::vector<Point> &box, const double lo0, const double hi0,
                        const unsigned rounds) {
            PushResult best;
            double lo = lo0, hi = hi0;

            hp[i].c = FT(hi);
            if (auto v = ev(hp, box)) {  // unconstrained in this direction
                best = {true, FT(hi), v->first, std::move(v->second)};
                return best;
            }
            for (unsigned r = 0; r < rounds; ++r) {
                const double mid = (lo + hi) / 2.0;
                hp[i].c = FT(mid);
                if (auto v = ev(hp, box)) {
                    lo = mid;
                    best = {true, FT(mid), v->first, std::move(v->second)};
                } else {
                    hi = mid;
                }
            }
            return best;
        }

        // Drop half-planes that constrain nothing: removing one leaves the realised
        // area unchanged.
        std::size_t prune(std::vector<HalfPlane> &hp, const std::vector<Point> &box) {
            const FT area = path_area(realise(hp, box));
            std::size_t dropped = 0;
            for (std::size_t i = hp.size(); i-- > 0;) {
                if (hp.size() <= 3) break;
                std::vector<HalfPlane> without;
                without.reserve(hp.size() - 1);
                for (std::size_t k = 0; k < hp.size(); ++k)
                    if (k != i) without.push_back(hp[k]);
                if (path_area(realise(without, box)) == area) {
                    hp.swap(without);
                    ++dropped;
                }
            }
            return dropped;
        }

        // The point on half-plane i's line to rotate about: the midpoint of Q's
        // edge on it. If no edge of Q lies on it, the foot of the perpendicular
        // from the centroid.
        Point rotation_anchor(const std::vector<Point> &pts, const HalfPlane &h) {
            std::vector<Point> on;
            for (const auto &p : pts)
                if (h.a * p.x() + h.b * p.y() == h.c) on.push_back(p);
            if (on.size() >= 2) return CGAL::midpoint(on.front(), on.back());
            if (on.size() == 1) return on.front();

            FT cx{0}, cy{0};
            for (const auto &p : pts) {
                cx += p.x();
                cy += p.y();
            }
            const FT k(static_cast<int>(pts.size()));
            cx /= k;
            cy /= k;
            const FT nn = h.a * h.a + h.b * h.b;
            if (nn == FT(0)) return {cx, cy};
            const FT t = (h.c - h.a * cx - h.b * cy) / nn;
            return {cx + h.a * t, cy + h.b * t};
        }
    }  // namespace

    Result refine(const Polygon &p, const Polygon &start, const Options &opts) {
        Result res;
        const Polygon poly = counterclockwise(p);
        const Polygon q0 = counterclockwise(start);
        if (q0.size() < 3) {
            res.message = "start polygon has fewer than 3 vertices";
            return res;
        }

        const auto bb = poly.bbox();
        const double pad = (bb.xmax() - bb.xmin()) + (bb.ymax() - bb.ymin()) + 1.0;
        const std::vector<Point> box
            = {Point(bb.xmin() - pad, bb.ymin() - pad), Point(bb.xmax() + pad, bb.ymin() - pad),
               Point(bb.xmax() + pad, bb.ymax() + pad), Point(bb.xmin() - pad, bb.ymax() + pad)};

        Evaluator ev(poly);
        std::vector<HalfPlane> hp = halfplanes_of(q0);

        auto init = ev(hp, box);
        if (!init) {
            res.message = "start polygon is not a convex subset of P";
            res.evaluations = ev.calls();
            return res;
        }
        FT area = init->first;
        std::vector<Point> pts = std::move(init->second);
        res.start_area = area;
        res.halfplanes_start = hp.size();

        // A convex P is its own answer. The start has still been validated above.
        if (poly.is_convex()) {
            res.area = CGAL::abs(poly.area());
            res.polygon = poly;
            res.halfplanes_end = poly.size();
            res.sweeps = 0;
            res.evaluations = ev.calls();
            res.ok = true;
            return res;
        }
        const std::size_t hp_cap = hp.size() * std::max(1u, opts.max_growth);

        const double c_hi_span = pad * pad * 8.0;  // generous "no constraint" offset

        for (unsigned sweep = 0; sweep < opts.max_sweeps; ++sweep) {
            res.sweeps = sweep + 1;
            bool improved = false;
            std::vector<bool> pinned_edge;
            pinned_edge.reserve(hp.size());

            for (std::size_t i = 0; i < hp.size(); ++i) {
                const double norm
                    = std::sqrt(CGAL::to_double(hp[i].a * hp[i].a + hp[i].b * hp[i].b));
                const double c_now = CGAL::to_double(hp[i].c);

                // ---- PUSH: exact bisection on the offset ---------------------
                bool pinned = true;
                if (const PushResult pr
                    = push(ev, hp, i, box, c_now, c_now + norm * c_hi_span, opts.push_rounds);
                    pr.ok && pr.area > area) {
                    hp[i].c = pr.c;
                    area = pr.area;
                    pts = pr.pts;
                    improved = true;
                    pinned = false;
                }
                pinned_edge.push_back(pinned);

                // ---- ROTATE: 1-D search on the normal, each trial re-pushed ---
                const double phi = std::atan2(CGAL::to_double(hp[i].b), CGAL::to_double(hp[i].a));
                const Point anchor = rotation_anchor(pts, hp[i]);
                const int half = static_cast<int>(opts.rot_steps) / 2;
                std::optional<HalfPlane> best_hp;
                FT best_area = area;
                std::vector<Point> best_pts;

                for (int k = -half; k <= half; ++k) {
                    if (k == 0 || half == 0) continue;
                    const double th = phi + opts.rot_span * k / half;
                    const FT na(std::cos(th)), nb(std::sin(th));
                    if (na == FT(0) && nb == FT(0)) continue;

                    std::vector<HalfPlane> trial = hp;
                    const FT c0 = na * anchor.x() + nb * anchor.y();
                    trial[i] = {na, nb, c0};
                    const double c0d = CGAL::to_double(c0);
                    if (const PushResult pr
                        = push(ev, trial, i, box, c0d - pad, c0d + pad, opts.push_rounds);
                        pr.ok && pr.area > best_area) {
                        best_area = pr.area;
                        best_hp = HalfPlane{na, nb, pr.c};
                        best_pts = pr.pts;
                    }
                }
                if (best_hp) {
                    hp[i] = *best_hp;
                    area = best_area;
                    pts = best_pts;
                    improved = true;
                }
            }

            // ---- SPLIT: replace one edge with two, each pushed ---------------
            //
            // See Options::split_steps.
            if (opts.split_steps && hp.size() < hp_cap) {
                const std::size_t n_before = hp.size();
                // `i` walks the current hp, `orig` the edges before this pass
                // (which pinned_edge indexes). They diverge once a split is
                // accepted, so `i` is advanced explicitly.
                std::size_t i = 0;
                for (std::size_t orig = 0; orig < n_before && hp.size() < hp_cap; ++orig, ++i) {
                    // Only pinned edges are worth splitting; an edge the push just
                    // advanced was not at its limit.
                    if (orig < pinned_edge.size() && !pinned_edge[orig]) continue;
                    const double phi
                        = std::atan2(CGAL::to_double(hp[i].b), CGAL::to_double(hp[i].a));
                    const Point anchor = rotation_anchor(pts, hp[i]);
                    std::optional<std::vector<HalfPlane>> best_hp;
                    FT best_area = area;
                    std::vector<Point> best_pts;

                    for (unsigned k = 1; k <= opts.split_steps; ++k) {
                        const double d = opts.split_span * k / opts.split_steps;
                        std::vector<HalfPlane> trial;
                        trial.reserve(hp.size() + 1);
                        for (std::size_t j = 0; j < hp.size(); ++j) {
                            if (j != i) {
                                trial.push_back(hp[j]);
                                continue;
                            }
                            for (const double th : {phi - d, phi + d}) {
                                const FT na(std::cos(th)), nb(std::sin(th));
                                trial.push_back({na, nb, na * anchor.x() + nb * anchor.y()});
                            }
                        }
                        // push BOTH new half-planes; each can advance past the
                        // contact point that pinned the single edge
                        bool feasible = true;
                        FT got{0};
                        std::vector<Point> got_pts;
                        for (const std::size_t idx : {i, i + 1}) {
                            const double c0 = CGAL::to_double(trial[idx].c);
                            // Coarse; the next sweep's push refines it.
                            const PushResult pr = push(ev, trial, idx, box, c0 - pad, c0 + pad,
                                                       opts.split_push_rounds);
                            if (!pr.ok) {
                                feasible = false;
                                break;
                            }
                            trial[idx].c = pr.c;
                            got = pr.area;
                            got_pts = pr.pts;
                        }
                        if (feasible && got > best_area) {
                            best_area = got;
                            best_hp = trial;
                            best_pts = got_pts;
                        }
                    }
                    if (best_hp) {
                        hp = *best_hp;
                        area = best_area;
                        pts = best_pts;
                        improved = true;
                        ++res.splits_accepted;
                        ++i;  // one edge became two; step over the second
                    }
                }
            }

            // Splitting can leave redundant half-planes.
            prune(hp, box);

            // ---- ADJOIN: absorb a corner of P, adding a VERTEX ---------------
            //
            // See Options::adjoin_rounds. The candidates are P's corners, one
            // certify call each. Runs after split, because it rebuilds the
            // half-plane list and would invalidate pinned_edge.
            if (opts.adjoin_rounds) {
                for (unsigned round = 0; round < opts.adjoin_rounds; ++round) {
                    bool gained = false;
                    for (std::size_t vi = 0; vi < poly.size(); ++vi) {
                        if (hp.size() >= hp_cap) break;
                        const Point &v = poly[vi];
                        // already a vertex of Q: the hull would be unchanged
                        if (std::find(pts.begin(), pts.end(), v) != pts.end()) continue;
                        std::vector<Point> cand = pts;
                        cand.push_back(v);
                        auto got = ev.points(cand);
                        if (!got || got->first <= area) continue;
                        area = got->first;
                        pts.assign(got->second.vertices_begin(), got->second.vertices_end());
                        hp = halfplanes_of(got->second);
                        ++res.adjoins_accepted;
                        gained = true;
                        improved = true;
                    }
                    if (!gained) break;
                }
            }

            if (!improved) break;
        }

        res.halfplanes_end = hp.size();
        res.area = area;
        res.polygon.clear();
        for (const auto &v : pts) res.polygon.push_back(v);
        if (res.polygon.is_clockwise_oriented()) res.polygon.reverse_orientation();
        res.evaluations = ev.calls();
        res.ok = true;
        return res;
    }

}  // namespace potato_peeling::refinement

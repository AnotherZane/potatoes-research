#pragma once

// Memory instrumentation for the Chang & Yap solve, used by solve_memory.py.
//
// Memory, not time, bounds the exact solver. Snapshots at phase boundaries
// separate the two candidate costs: the per-entry Lazy_exact_nt DAGs (scaling
// with entries stored) and the matrices' bucket arrays (scaling with their
// reservation). `VmHWM` from /proc/self/status is the kernel's high-water mark of
// resident memory, which is what the OOM killer acts on; `VmRSS` is the current
// value.
//
// It never influences a result. The report is a process-wide static, so it is
// for one solve per process.

#include <cstddef>
#include <string>
#include <vector>

namespace potato_peeling::chang_and_yap::memprobe {

    // One phase boundary of the solve.
    struct Snapshot {
        std::string phase;
        std::size_t rss_bytes;       // VmRSS at this point
        std::size_t peak_rss_bytes;  // VmHWM so far
        double seconds;              // since the solve began

        // Structural state. 0 when the phase precedes the structure's creation.
        std::size_t xi = 0;         // |Xi|, the chord count
        std::size_t a_stored = 0;   // A entries actually present
        std::size_t a_buckets = 0;  // A's bucket array
        std::size_t m_stored = 0;   // M entries actually present
        std::size_t m_buckets = 0;  // M's bucket array
        std::size_t crit_l0 = 0;    // CritParamCache: l0 keys
        std::size_t crit_l1 = 0;    // CritParamCache: l1 keys
        std::size_t crit_ft = 0;    // CritParamCache: total FT values held

        // Heap attribution (glibc mallinfo2). RSS alone cannot tell live data
        // from memory freed but retained by the allocator.
        //   live  = uordblks, bytes actually allocated and not yet freed
        //   freed = fordblks, free chunks still held inside the arena
        //   arena = arena + hblkhd, total obtained from the OS
        // live / arena is the useful ratio: near 1 means the peak IS live data
        // and shrinking structures is the only lever; well below 1 means the
        // peak is retention or fragmentation and the lever is allocation
        // strategy. mallinfo2 aggregates over all arenas.
        std::size_t heap_live = 0;
        std::size_t heap_freed = 0;
        std::size_t heap_arena = 0;
    };

    struct Report {
        std::size_t n = 0;  // corners after pre_checks
        std::vector<Snapshot> snapshots;
        std::size_t peak_rss_bytes = 0;      // VmHWM at the end of the solve
        std::size_t baseline_rss_bytes = 0;  // VmRSS before any solve work

        [[nodiscard]] bool empty() const { return snapshots.empty(); }
    };

    // Current VmRSS in bytes; 0 if /proc is unavailable.
    std::size_t rss_bytes();

    // Peak VmHWM in bytes; 0 if /proc is unavailable. This is what the OOM
    // killer reacts to, so it is the figure a memory ceiling is stated in.
    std::size_t peak_rss_bytes();

    // Enable or disable collection; disabled by default.
    void set_enabled(bool on);

    bool enabled();

    // Begin a fresh report; records the baseline RSS.
    void begin(std::size_t n);

    // Record a phase boundary, optionally with structural counts from the caller.
    void mark(const std::string &phase);

    void mark(const std::string &phase, const Snapshot &structural);

    // Record CritParamCache's size while it is alive (it is local to
    // compute_A_matrix_all_stages). Counts rather than bytes: freed memory does
    // not show up as an RSS drop.
    void note_crit(std::size_t l0_keys, std::size_t l1_keys, std::size_t ft_values);

    // The most recent completed report.
    const Report &last();

}  // namespace potato_peeling::chang_and_yap::memprobe

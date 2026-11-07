#include "potato_peeling_core/algorithms/chang_and_yap/memprobe.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <utility>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace potato_peeling::chang_and_yap::memprobe {
    namespace {
        bool g_enabled = false;
        Report g_report;
        std::chrono::steady_clock::time_point g_start;
        // Latched by note_crit and copied into every subsequent snapshot, so the
        // cache's size survives its own destruction.
        std::size_t g_crit_l0 = 0, g_crit_l1 = 0, g_crit_ft = 0;

        // Parse one "VmXXX:  N kB" line out of /proc/self/status. Returns bytes,
        // or 0 when /proc is unavailable (non-Linux, or a restricted container),
        // which callers treat as "unknown" rather than "zero used".
        std::size_t proc_status_kib(const char *key) {
            std::FILE *f = std::fopen("/proc/self/status", "r");
            if (!f) return 0;
            char line[256];
            const std::size_t klen = std::strlen(key);
            std::size_t out = 0;
            while (std::fgets(line, sizeof(line), f)) {
                if (std::strncmp(line, key, klen) == 0) {
                    unsigned long long kib = 0;
                    if (std::sscanf(line + klen, " %llu", &kib) == 1)
                        out = static_cast<std::size_t>(kib) * 1024u;
                    break;
                }
            }
            std::fclose(f);
            return out;
        }
    }  // namespace

    std::size_t rss_bytes() {
        return proc_status_kib("VmRSS:");
    }

    std::size_t peak_rss_bytes() {
        return proc_status_kib("VmHWM:");
    }

    void set_enabled(const bool on) {
        g_enabled = on;
    }

    bool enabled() {
        return g_enabled;
    }

    void begin(const std::size_t n) {
        if (!g_enabled) return;
        g_report = Report{};
        g_report.n = n;
        g_report.baseline_rss_bytes = rss_bytes();
        g_crit_l0 = g_crit_l1 = g_crit_ft = 0;
        g_start = std::chrono::steady_clock::now();
    }

    void mark(const std::string &phase) {
        mark(phase, Snapshot{});
    }

    void mark(const std::string &phase, const Snapshot &structural) {
        if (!g_enabled) return;
        Snapshot s = structural;
        s.phase = phase;
        s.crit_l0 = g_crit_l0;
        s.crit_l1 = g_crit_l1;
        s.crit_ft = g_crit_ft;
#if defined(__GLIBC__) && defined(__GLIBC_PREREQ)
#if __GLIBC_PREREQ(2, 33)
        {
            const struct mallinfo2 mi = mallinfo2();
            s.heap_live = mi.uordblks;
            s.heap_freed = mi.fordblks;
            s.heap_arena = mi.arena + mi.hblkhd;
        }
#endif
#endif
        s.rss_bytes = rss_bytes();
        s.peak_rss_bytes = peak_rss_bytes();
        s.seconds
            = std::chrono::duration<double>(std::chrono::steady_clock::now() - g_start).count();
        g_report.peak_rss_bytes = s.peak_rss_bytes;
        g_report.snapshots.push_back(std::move(s));
    }

    void note_crit(const std::size_t l0_keys, const std::size_t l1_keys,
                   const std::size_t ft_values) {
        if (!g_enabled) return;
        g_crit_l0 = l0_keys;
        g_crit_l1 = l1_keys;
        g_crit_ft = ft_values;
    }

    const Report &last() {
        return g_report;
    }

}  // namespace potato_peeling::chang_and_yap::memprobe

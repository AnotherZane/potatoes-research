#pragma once

// Seeded, position-addressed random substreams.
//
// There is no shared generator drawn from in sequence. Every draw belongs to a
// substream identified by (seed, purpose, i, j), so what one edge of one
// iteration sees does not depend on earlier draws, and a run replays identically
// for a given seed whatever the number of threads.

#include <cstdint>
#include <random>

namespace potato_peeling::cabello {

    // Purpose tags. Substreams differing only in purpose are independent.
    enum class Stream : std::uint64_t {
        SampleR = 1,        // line 5: the r points sampled in P
        SampleRab = 2,      // line 9: R_ab, sampled in Gamma
        SampleSab = 3,      // line 10: S_ab, sampled in Gamma
        EdgeSubsample = 4,  // choosing which edges of G(P,R) to search
    };

    inline std::uint64_t splitmix64(std::uint64_t x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

    // The substream for (seed, purpose, i, j): `i` is the iteration and `j` the
    // edge index.
    inline std::mt19937_64 substream(const std::uint64_t seed, const Stream purpose,
                                     const std::uint64_t i = 0, const std::uint64_t j = 0) {
        std::uint64_t h = splitmix64(seed);
        h = splitmix64(h ^ (static_cast<std::uint64_t>(purpose) * 0xD1B54A32D192ED03ULL));
        h = splitmix64(h ^ (i * 0xA24BAED4963EE407ULL));
        h = splitmix64(h ^ (j * 0x9FB21C651E98DF25ULL));
        const std::uint64_t h2 = splitmix64(h);

        std::seed_seq seq{static_cast<std::uint32_t>(h), static_cast<std::uint32_t>(h >> 32),
                          static_cast<std::uint32_t>(h2), static_cast<std::uint32_t>(h2 >> 32)};
        return std::mt19937_64(seq);
    }

    // Uniform on [0, 1) from the top 53 bits of one draw. Not
    // std::generate_canonical, which can return 1.0 in libstdc++ and differs
    // between standard libraries.
    inline double canonical(std::mt19937_64 &gen) {
        return static_cast<double>(gen() >> 11) * 0x1.0p-53;
    }

}  // namespace potato_peeling::cabello

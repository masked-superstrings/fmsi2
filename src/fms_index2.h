#pragma once

// FMSI 2: the masked Burrows-Wheeler transform in a cache-friendly
// layout, searched in batches with prefetching.
//
// Differences from the original index (fms_index.h):
//  - the BWT is stored as 2-bit symbols with in-line rank counts (one cache
//    line per rank instead of a wavelet tree over three bitvectors),
//  - the SA-transformed mask is either a plain bitvector with in-line rank
//    counts or the original RRR vector,
//  - the SA intervals of all short patterns are tabulated when the index is loaded,
//  - k-mers are searched in lockstep, each step prefetching the cache lines of the next one.

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include <sdsl/bit_vectors.hpp>

#include "fms_index.h"
#include "rank2.h"

namespace fmsi2 {

enum class mask_kind : uint32_t {
    plain = 0,
    rrr = 1,
};

struct index {
    uint32_t k = 0;
    /// Length of the BWT, i.e., the length of the superstring plus one.
    uint64_t size = 0;
    uint64_t dollar_position = 0;
    uint64_t counts[4] = {0, 0, 0, 0};
    quad_rank bwt;
    mask_kind mask_type = mask_kind::plain;
    bit_rank mask;
    sdsl::rrr_vector<RRR_BLOCK_SIZE> mask_rrr;
    sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE> mask_rrr_rank;

    /// SA intervals of all patterns of length prefix_length; computed at load time, not stored.
    uint32_t prefix_length = 0;
    std::vector<uint64_t> prefix_table;

    /// The number of set bits of the mask, i.e., of hashes; computed when the index is built or loaded.
    uint64_t mask_ones = 0;

    /// Bytes of the index in memory: the BWT, the mask and the table of patterns.
    size_t bytes() const { return bwt.bytes() + mask_bytes() + prefix_table.size() * 8; }

    index() = default;
    index(const index &) = delete;
    index &operator=(const index &) = delete;

    inline bool mask_get(size_t i) const {
        return mask_type == mask_kind::plain ? mask[i] : (bool)mask_rrr[i];
    }

    inline size_t mask_rank(size_t i) const {
        return mask_type == mask_kind::plain ? mask.rank(i) : mask_rrr_rank(i);
    }

    inline void mask_prefetch(size_t i) const {
        if (mask_type == mask_kind::plain) mask.prefetch(i);
    }

    size_t mask_bytes() const {
        return mask_type == mask_kind::plain ? mask.bytes() : sdsl::size_in_bytes(mask_rrr);
    }
};

/// Go from range (i,j) for pattern P to range for c+P.
inline void update_range(const index &ix, size_t &i, size_t &j, byte c) {
    size_t count = ix.counts[c];
    size_t rank_i, rank_j;
    ix.bwt.rank_pair(i, j, c, rank_i, rank_j);
    // The dollar is stored as A in the BWT.
    i = count + rank_i - (c == 0 && i > ix.dollar_position);
    j = count + rank_j - (c == 0 && j > ix.dollar_position);
}

/// Build the index from the MBWT, which is consumed.
void build(index &ix, mbwt &m, int k, mask_kind mask_type) {
    ix.k = k;
    ix.size = m.bwt.size();
    ix.dollar_position = m.dollar_position;

    ix.bwt = quad_rank(m.bwt.size());
    uint64_t totals[4] = {0, 0, 0, 0};
    for (size_t i = 0; i < m.bwt.size(); ++i) {
        ix.bwt.set(i, m.bwt[i]);
        totals[m.bwt[i]]++;
    }
    ix.bwt.finalize();
    m.bwt = std::vector<byte>();
    // The count of A includes the dollar, which stands for the smallest suffix.
    ix.counts[0] = 1;
    ix.counts[1] = totals[0];
    ix.counts[2] = totals[0] + totals[1];
    ix.counts[3] = totals[0] + totals[1] + totals[2];

    ix.mask_type = mask_type;
    if (mask_type == mask_kind::plain) {
        ix.mask = bit_rank(m.sa_transformed_mask.size());
        for (size_t i = 0; i < m.sa_transformed_mask.size(); ++i) {
            if (m.sa_transformed_mask[i]) ix.mask.set(i);
        }
        ix.mask.finalize();
    } else {
        ix.mask_rrr = sdsl::rrr_vector<RRR_BLOCK_SIZE>(m.sa_transformed_mask);
        ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
    }
    m.sa_transformed_mask.resize(0);
    ix.mask_ones = ix.mask_rank(ix.size);
}

/// Tabulate the SA intervals of all patterns of the given length. An entry is indexed by the
/// characters in the order of the backward search, the first one being the most significant.
void build_prefix_table(index &ix, unsigned length) {
    length = std::min(length, ix.k);
    std::vector<uint64_t> current = {0, ix.size};
    for (unsigned depth = 0; depth < length; ++depth) {
        std::vector<uint64_t> next(current.size() * 4);
        for (size_t key = 0; key < current.size() / 2; ++key) {
            for (byte c = 0; c < 4; ++c) {
                size_t i = current[2 * key], j = current[2 * key + 1];
                if (i != j) update_range(ix, i, j, c);
                next[2 * (4 * key + c)] = i;
                next[2 * (4 * key + c) + 1] = j;
            }
        }
        current = std::move(next);
    }
    ix.prefix_length = length;
    ix.prefix_table = std::move(current);
}

struct lookup_stats {
    uint64_t kmers = 0;
    uint64_t found = 0;
    /// Number of steps of the backward search (LF-mappings of an interval).
    uint64_t steps = 0;
    /// Number of searches started from scratch.
    uint64_t searches = 0;
    double seconds = 0;

    void add(const lookup_stats &other) {
        kmers += other.kmers;
        found += other.found;
        steps += other.steps;
        searches += other.searches;
    }
};

/// Look up the k-mers starting at the given offsets of the sequence, [batch_size] of them at a time.
/// For each k-mer, out[offset] is set to its presence (0/1), or to its hash (-1 if absent) if [hashes].
/// A k-mer is represented if the k-mer or its reverse complement has an occurrence set in the mask.
template <bool hashes, typename out_t>
void lookup_batch(const index &ix, const char *sequence, const uint32_t *offsets, size_t count, out_t *out,
                  unsigned batch_size, lookup_stats &stats) {
    const unsigned k = ix.k;
    const unsigned p = ix.prefix_length;
    const uint64_t *table = ix.prefix_table.data();

    enum : uint8_t { IDLE, START, SEARCH, MASK };
    struct slot {
        size_t i, j;
        const char *kmer;
        uint64_t key;
        uint32_t offset;
        /// Number of characters of the k-mer already searched.
        uint32_t t;
        uint8_t state;
        /// Whether the reverse complement is being searched, and whether it is the second strand tried.
        uint8_t reverse, second;
    };
    std::vector<slot> slots(std::max(1u, batch_size));
    size_t next = 0, active = slots.size();
    // Saturating counter predicting the strand of the next k-mer from the previous hits.
    int score = 0;
    uint64_t found = 0, steps = 0, searches = 0;

    auto begin_strand = [&](slot &s) {
        ++searches;
        if (p == 0) {
            s.i = 0;
            s.j = ix.size;
            s.t = 0;
            s.state = SEARCH;
            return;
        }
        uint64_t key = 0;
        if (s.reverse) {
            for (unsigned t = 0; t < p; ++t) key = (key << 2) | (3 - nucleotideToInt[(uint8_t)s.kmer[t]]);
        } else {
            for (unsigned t = 0; t < p; ++t) key = (key << 2) | nucleotideToInt[(uint8_t)s.kmer[k - 1 - t]];
        }
        s.key = key;
        s.state = START;
        __builtin_prefetch(table + 2 * key);
    };
    auto load = [&](slot &s) {
        if (next == count) {
            s.state = IDLE;
            --active;
            return;
        }
        s.offset = offsets[next++];
        s.kmer = sequence + s.offset;
        s.reverse = score < 0;
        s.second = 0;
        begin_strand(s);
    };
    auto strand_failed = [&](slot &s) {
        if (s.second) {
            out[s.offset] = hashes ? -1 : 0;
            load(s);
        } else {
            s.second = 1;
            s.reverse ^= 1;
            begin_strand(s);
        }
    };
    auto advance = [&](slot &s) {
        if (s.i == s.j) {
            strand_failed(s);
        } else if (s.t == k) {
            s.state = MASK;
            ix.mask_prefetch(s.i);
        } else {
            s.state = SEARCH;
            ix.bwt.prefetch(s.i);
            ix.bwt.prefetch(s.j);
        }
    };

    for (slot &s : slots) load(s);
    while (active > 0) {
        for (slot &s : slots) {
            switch (s.state) {
            case START:
                s.i = table[2 * s.key];
                s.j = table[2 * s.key + 1];
                s.t = p;
                advance(s);
                break;
            case SEARCH: {
                byte c = s.reverse ? 3 - nucleotideToInt[(uint8_t)s.kmer[s.t]]
                                   : nucleotideToInt[(uint8_t)s.kmer[k - 1 - s.t]];
                update_range(ix, s.i, s.j, c);
                ++s.t;
                ++steps;
                advance(s);
                break;
            }
            case MASK: {
                bool represented = false;
                for (size_t q = s.i; q < s.j && !represented; ++q) represented = ix.mask_get(q);
                if (represented) {
                    out[s.offset] = hashes ? (out_t)ix.mask_rank(s.i) : (out_t)1;
                    ++found;
                    score = std::max(-7, std::min(7, score + (s.reverse ? -1 : 1)));
                    load(s);
                } else {
                    strand_failed(s);
                }
                break;
            }
            default:
                break;
            }
        }
    }

    stats.kmers += count;
    stats.found += found;
    stats.steps += steps;
    stats.searches += searches;
}

} // namespace fmsi2

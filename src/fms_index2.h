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
#include "index_file.h"
#include "rank2.h"
#include "version.h"

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

    /// The type of the mask as the index file states it: "min-one" (every k-mer of the set has one set
    /// occurrence, so that lookup is a minimal perfect hash), "max-one" (all occurrences of the k-mers of the
    /// set are set), "min-one,max-one", "other", or "" if unknown; and how it is known ("computed", "declared").
    std::string mask_type_name, mask_type_source;
    /// The FASTA header of the masked superstring that was indexed, if the index file has it.
    std::string superstring_header;

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

/// The type of a mask. FMSI answers membership queries the fastest on a max-one mask (all occurrences of the
/// k-mers of the set are set, so the first bit of an interval decides) and dictionary queries on a min-one mask
/// (one set occurrence per k-mer, so the rank of the interval is a minimal perfect hash). FMSI 2 tests any mask
/// so the type matters for the hashes only: with a mask other than min-one they stay unique
/// per set occurrence, but a k-mer set several times has several, and their range exceeds the number of k-mers.
///
/// The type is computed from the rows that start with the same k-mer (mbwt::same_kmer). A k-mer and its reverse
/// complement are not related to each other there: a mask found min-one may still set a k-mer once in each
/// orientation (masks of KmerCamel in its bidirectional mode never do). The type that the producer of the masked
/// superstring states in its FASTA header (mask=min-one) or the user gives is therefore preferred, and what is
/// computed checks it.
struct mask_statistics {
    /// "min-one", "max-one", "min-one,max-one" (each k-mer of the set occurs once), "other" or "unknown".
    std::string type = "unknown";
    /// "computed-one-strand", "header" (of the masked superstring), "declared" (by the user), or empty.
    std::string source;
    /// Whether the type was computed, and what was found: no k-mer with several set occurrences; no k-mer of
    /// the set with an unset occurrence.
    bool computed = false, single = false, uniform = false;
    /// Number of different k-mers of the text with a set occurrence (computed only).
    uint64_t kmers = 0;
    /// Set if a stated type contradicts what was computed.
    std::string warning;
};

/// The type of the mask computed from the MBWT (needs mbwt::same_kmer; "unknown" without it).
inline mask_statistics classify_mask(const mbwt &m) {
    mask_statistics result;
    if (m.bwt.empty() || m.same_kmer.size() != m.bwt.size()) return result;
    bool single = true, uniform = true;
    uint64_t ones = 0, length = 0;
    const size_t rows = m.bwt.size();
    for (size_t i = 0; i < rows; ++i) {
        ones += m.sa_transformed_mask[i];
        ++length;
        if (!m.same_kmer[i]) {
            if (ones > 0) {
                ++result.kmers;
                single = single && ones == 1;
                uniform = uniform && ones == length;
            }
            ones = length = 0;
        }
    }
    result.computed = true;
    result.single = single;
    result.uniform = uniform;
    result.type = single && uniform ? "min-one,max-one" : single ? "min-one" : uniform ? "max-one" : "other";
    result.source = "computed-one-strand";
    return result;
}

/// The type of the mask given what was computed and what was stated: by the user ([declared]: "min-one",
/// "max-one", or "auto" for nothing) or in the header of the masked superstring ([header_type], may be empty).
/// A stated type is kept unless the computed one contradicts it.
inline mask_statistics decide_mask_type(mask_statistics computed, const std::string &declared, const std::string &header_type) {
    const bool by_user = declared == "min-one" || declared == "max-one";
    const bool by_header = !by_user && (header_type == "min-one" || header_type == "max-one");
    if (!by_user && !by_header) return computed;
    const std::string stated = by_user ? declared : header_type;
    if (computed.computed && ((stated == "min-one" && !computed.single) || (stated == "max-one" && !computed.uniform))) {
        computed.warning = "the mask is stated to be " + stated + (by_user ? " (-M)" : " (in the header of the masked superstring)") +
                           ", but it is not: " + (stated == "min-one" ? "a k-mer has several set occurrences" : "a k-mer of the set has an unset occurrence") +
                           "; its type is recorded as " + computed.type;
        return computed;
    }
    computed.type = stated;
    computed.source = by_user ? "declared" : "header";
    return computed;
}

/// Build the index from the MBWT, which is consumed.
void build(index &ix, mbwt &m, int k, mask_kind mask_type) {
    ix.k = k;
    ix.size = m.bwt.size();
    ix.dollar_position = m.dollar_position;
    const mask_statistics statistics = classify_mask(m);
    ix.mask_type_name = statistics.type == "unknown" ? "" : statistics.type;
    ix.mask_type_source = statistics.source;

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

/// The layout of an index in memory. It is chosen when the index is loaded, as the file holds the BWT and the
/// mask in a form of its own; the file recommends one.
struct layout {
    /// Encoding of the mask (a mask_kind); -1 stands for the recommended one.
    int mask = -1;
};

/// The layout that a file recommends unless fmsi2 index is told otherwise (-m): the one that answers the
/// fastest (the plain mask).
inline layout default_layout() {
    layout result;
    result.mask = (int)mask_kind::plain;
    return result;
}

inline const char *mask_name(mask_kind kind) {
    return kind == mask_kind::plain ? "plain" : "rrr";
}

inline bool parse_mask_name(const std::string &name, mask_kind &kind) {
    if (name == "plain") {
        kind = mask_kind::plain;
    } else if (name == "rrr") {
        kind = mask_kind::rrr;
    } else {
        return false;
    }
    return true;
}

/// Format 2 (index_file.h): one file with tagged sections. The sections of an index:
///   META  text, key=value per line (below)
///   BWT   encoding 1: the characters of the BWT, 2 bits each, 4 per byte (character i at bits 2 (i % 4) of
///         byte i / 4), the sentinel stored as A; not compressed, as it does not compress
///   MASK  encoding 1: the mask in the order of the suffix array, one bit per row (bit i at bit i % 8 of byte
///         i / 8); compressed
/// The rank counts of the BWT and of the mask and the table of the first search steps are not stored: they are
/// computed when the file is loaded, in the layout asked for.
///
/// Keys of META of an index: kind (index), creator, k, rows (the length of the BWT), dollar (the row of the
/// sentinel), counts (the four cumulative character counts), mask_ones, mask_order (sa), layout (the recommended
/// layout in memory: the mask encoding). Unknown keys are ignored.
static const uint32_t ENCODING_BWT_PACKED = 1;
static const uint32_t ENCODING_MASK_SA_BITS = 1;
static const char *const INDEX_EXTENSION = ".fmsi";

/// The mask as a bitvector, whatever its encoding.
inline sdsl::bit_vector mask_bits(const index &ix) {
    sdsl::bit_vector bits(ix.size, 0);
    for (size_t i = 0; i < ix.size; ++i) {
        if (ix.mask_get(i)) bits[i] = 1;
    }
    return bits;
}

/// What a file says about its index besides the BWT and the mask.
struct description {
    uint32_t k = 0;
    uint64_t rows = 0;
    uint64_t dollar_position = 0;
    uint64_t counts[4] = {0, 0, 0, 0};
    uint64_t mask_ones = 0;
    layout recommended;
    /// The type of the mask and how it is known (mask_statistics), and the FASTA header of the masked superstring.
    std::string mask_type, mask_type_source, superstring_header;
};

inline std::string creator() { return std::string("fmsi2 ") + VERSION; }

/// Write an index file. bwt(begin, end, buffer) gives the bytes [begin, end) of the packed BWT, mask(begin, end,
/// buffer) those of the mask; both are called from several threads.
inline void write_index_file(const std::string &path, const description &d,
                             const std::function<void(uint64_t, uint64_t, uint8_t *)> &bwt,
                             const std::function<void(uint64_t, uint64_t, uint8_t *)> &mask, unsigned threads = 1) {
    file::metadata meta;
    meta.set("kind", "index");
    meta.set("creator", creator());
    meta.set("k", d.k);
    meta.set("rows", d.rows);
    meta.set("dollar", d.dollar_position);
    meta.set("counts", std::to_string(d.counts[0]) + "," + std::to_string(d.counts[1]) + "," + std::to_string(d.counts[2]) +
                           "," + std::to_string(d.counts[3]));
    meta.set("mask_ones", d.mask_ones);
    meta.set("mask_type", d.mask_type.empty() ? "unknown" : d.mask_type);
    if (!d.mask_type_source.empty()) meta.set("mask_type_source", d.mask_type_source);
    if (!d.superstring_header.empty()) meta.set("superstring_header", d.superstring_header);
    meta.set("mask_order", "sa");
    meta.set("layout", mask_name((mask_kind)d.recommended.mask));
    const std::string text = meta.text();

    file::writer out(path);
    out.add("META", 1, file::REQUIRED, meta.items().size(), text.data(), text.size());
    out.add_stream("BWT", ENCODING_BWT_PACKED, file::REQUIRED, d.rows, (d.rows + 3) / 4, false, threads, bwt);
    out.add_stream("MASK", ENCODING_MASK_SA_BITS, file::REQUIRED, d.rows, (d.rows + 7) / 8, true, threads, mask);
    out.finish();
}

/// Bytes [begin, end) of a bit vector of sdsl, bit i at bit i % 8 of byte i / 8.
inline void bit_vector_bytes(const sdsl::bit_vector &bits, uint64_t begin, uint64_t end, uint8_t *buffer) {
    // The words of the vector cover its bits; a byte past them is zero.
    const uint64_t available = ((bits.size() + 63) / 64) * 8;
    const uint64_t stop = std::min(end, available);
    if (begin < stop) memcpy(buffer, reinterpret_cast<const uint8_t *>(bits.data()) + begin, stop - begin);
    if (stop < end) memset(buffer + (std::max(begin, stop) - begin), 0, end - std::max(begin, stop));
}

/// Write the index of an MBWT to a file, without building the index in memory. Returns its description.
inline description save_mbwt(const mbwt &m, int k, layout recommended, const mask_statistics &statistics,
                             const std::string &path, unsigned threads = 1, const std::string &superstring_header = "") {
    description d;
    d.k = k;
    d.rows = m.bwt.size();
    d.dollar_position = m.dollar_position;
    d.recommended = recommended;
    d.mask_type = statistics.type;
    d.mask_type_source = statistics.source;
    // One line of the metadata.
    d.superstring_header = superstring_header;
    std::replace(d.superstring_header.begin(), d.superstring_header.end(), '\n', ' ');
    uint64_t totals[4] = {0, 0, 0, 0};
    for (size_t i = 0; i < m.bwt.size(); ++i) totals[m.bwt[i]]++;
    // The count of A includes the sentinel, which stands for the smallest suffix.
    d.counts[0] = 1;
    d.counts[1] = totals[0];
    d.counts[2] = totals[0] + totals[1];
    d.counts[3] = totals[0] + totals[1] + totals[2];
    const uint64_t *words = m.sa_transformed_mask.data();
    for (size_t w = 0; w < (m.sa_transformed_mask.size() + 63) / 64; ++w) d.mask_ones += __builtin_popcountll(words[w]);
    write_index_file(
        path, d,
        [&](uint64_t begin, uint64_t end, uint8_t *buffer) {
            for (uint64_t b = begin; b < end; ++b) {
                uint8_t packed = 0;
                for (uint64_t i = 4 * b; i < std::min<uint64_t>(4 * b + 4, d.rows); ++i) packed |= m.bwt[i] << (2 * (i & 3));
                buffer[b - begin] = packed;
            }
        },
        [&](uint64_t begin, uint64_t end, uint8_t *buffer) { bit_vector_bytes(m.sa_transformed_mask, begin, end, buffer); }, threads);
    return d;
}

inline description describe(const index &ix) {
    description d;
    d.k = ix.k;
    d.rows = ix.size;
    d.dollar_position = ix.dollar_position;
    for (int c = 0; c < 4; ++c) d.counts[c] = ix.counts[c];
    d.mask_ones = ix.mask_ones;
    d.recommended.mask = (int)ix.mask_type;
    d.mask_type = ix.mask_type_name;
    d.mask_type_source = ix.mask_type_source;
    d.superstring_header = ix.superstring_header;
    return d;
}

/// Write an index that is in memory to a file; the file recommends its layout, or the given one.
inline void save(const index &ix, const std::string &path, unsigned threads = 1, const layout *recommended = nullptr) {
    description d = describe(ix);
    if (recommended != nullptr) d.recommended = *recommended;
    write_index_file(
        path, d, [&](uint64_t begin, uint64_t end, uint8_t *buffer) { ix.bwt.get_packed(begin, buffer, end - begin); },
        [&](uint64_t begin, uint64_t end, uint8_t *buffer) {
            for (uint64_t b = begin; b < end; ++b) {
                uint8_t byte = 0;
                for (uint64_t i = 8 * b; i < std::min<uint64_t>(8 * b + 8, ix.size); ++i) byte |= (uint8_t)ix.mask_get(i) << (i & 7);
                buffer[b - begin] = byte;
            }
        },
        threads);
}

/// The bits of the pieces of a section, 64 at a time: sink(first bit, word, count).
template <typename sink_t>
inline void for_each_word(uint64_t offset, const uint8_t *data, size_t bytes, uint64_t total_bits, sink_t sink) {
    if (offset % 8 != 0) throw std::logic_error("pieces of a section must start at whole words");
    for (size_t b = 0; b < bytes && 8 * (offset + b) < total_bits; b += 8) {
        uint64_t word = 0;
        memcpy(&word, data + b, std::min<size_t>(8, bytes - b));
        const uint64_t first = 8 * (offset + b);
        sink(first, word, (unsigned)std::min<uint64_t>(64, total_bits - first));
    }
}

/// The description of an index file.
inline description read_description(const file::reader &in) {
    const file::metadata meta = in.meta();
    if (meta.get("kind") != "index") throw std::invalid_argument("file " + in.path() + " is not an index of FMSI 2");
    description d;
    d.k = (uint32_t)meta.number("k");
    d.rows = meta.number("rows");
    d.dollar_position = meta.number("dollar");
    d.mask_ones = meta.number("mask_ones");
    d.mask_type = meta.get("mask_type", "unknown");
    d.mask_type_source = meta.get("mask_type_source");
    d.superstring_header = meta.get("superstring_header");
    std::istringstream counts(meta.get("counts"));
    std::string item;
    for (int c = 0; c < 4; ++c) {
        if (!std::getline(counts, item, ',')) throw std::invalid_argument("file " + in.path() + " has no character counts");
        d.counts[c] = std::stoull(item);
    }
    d.recommended = default_layout();
    mask_kind kind;
    if (parse_mask_name(meta.get("layout"), kind)) d.recommended.mask = (int)kind;
    if (meta.get("mask_order", "sa") != "sa") {
        throw std::invalid_argument("file " + in.path() + " stores its mask in an order that this version of FMSI does not know");
    }
    return d;
}

/// Load an index in the given layout (by default the one the file recommends). The BWT is read into its blocks
/// and the mask decoded into its encoding piece by piece, so that the peak memory is that of the loaded index;
/// the sections are read on the given threads.
inline void load(index &ix, const std::string &path, unsigned threads = 1, layout wanted = layout()) {
    file::reader in(path);
    in.check_known({"META", "BWT", "MASK"});
    const description d = read_description(in);
    ix.k = d.k;
    ix.size = d.rows;
    ix.dollar_position = d.dollar_position;
    for (int c = 0; c < 4; ++c) ix.counts[c] = d.counts[c];
    ix.mask_type_name = d.mask_type == "unknown" ? "" : d.mask_type;
    ix.mask_type_source = d.mask_type_source;
    ix.superstring_header = d.superstring_header;
    ix.mask_type = (mask_kind)(wanted.mask >= 0 ? wanted.mask : d.recommended.mask);

    const file::section_entry &bwt = in.get("BWT");
    if (bwt.encoding != ENCODING_BWT_PACKED || bwt.items != d.rows || bwt.raw_bytes != (d.rows + 3) / 4) {
        throw std::invalid_argument("file " + path + " stores its BWT in a form that this version of FMSI does not know");
    }
    ix.bwt = quad_rank(d.rows);
    in.read(bwt, threads, [&](uint64_t offset, const uint8_t *data, size_t bytes) { ix.bwt.set_packed(offset, data, bytes); });
    ix.bwt.finalize();

    const file::section_entry &mask = in.get("MASK");
    if (mask.encoding != ENCODING_MASK_SA_BITS || mask.items != d.rows || mask.raw_bytes != (d.rows + 7) / 8) {
        throw std::invalid_argument("file " + path + " stores its mask in a form that this version of FMSI does not know");
    }
    if (ix.mask_type == mask_kind::plain) {
        ix.mask = bit_rank(d.rows);
        in.read(mask, threads, [&](uint64_t offset, const uint8_t *data, size_t bytes) {
            for_each_word(offset, data, bytes, d.rows, [&](uint64_t first, uint64_t word, unsigned count) { ix.mask.set_bits(first, word, count); });
        });
        ix.mask.finalize();
    } else {
        sdsl::bit_vector bits(d.rows, 0);
        in.read(mask, threads, [&](uint64_t offset, const uint8_t *data, size_t bytes) {
            for_each_word(offset, data, bytes, d.rows, [&](uint64_t first, uint64_t word, unsigned count) { bits.set_int(first, word, count); });
        });
        ix.mask_rrr = sdsl::rrr_vector<RRR_BLOCK_SIZE>(bits);
        ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
    }
    ix.mask_ones = ix.mask_rank(ix.size);
    if (ix.mask_ones != d.mask_ones) throw std::invalid_argument("file " + path + " is damaged: the mask does not match its description");
}

/// The file of an index given as the file itself, or as the name without its extension (a prefix, as the
/// commands of FMSI took it). Empty if there is none.
inline std::string find_index_file(const std::string &name) {
    auto is_file = [](const std::string &path) { return file::file_version(path) != 0; };
    if (is_file(name)) return name;
    if (is_file(name + INDEX_EXTENSION)) return name + INDEX_EXTENSION;
    return "";
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
    /// Number of searches started from scratch, and those of them caused by an unset mask bit of a long match.
    uint64_t searches = 0;
    uint64_t rechecks = 0;
    double seconds = 0;

    void add(const lookup_stats &other) {
        kmers += other.kmers;
        found += other.found;
        steps += other.steps;
        searches += other.searches;
        rechecks += other.rechecks;
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

/// A stretch of a query sequence without invalid characters, at least k characters long.
struct segment {
    uint32_t offset, length;
};

/// Look up all k-mers of the segments by streaming, [batch_size] segments at a time; the results
/// are as in lookup_batch (out must be initialized to 0, or to -1 if [hashes]).
///
/// A segment is searched on one strand and the k-mers not found there on the other one. On a strand,
/// the k-mers (windows) are resolved from the right to the left: the match of a window is extended by one
/// character to the match of the next one, so the match can get longer than k. A longer match has a
/// subset of the occurrences of the k-mer, hence a set mask bit proves the k-mer represented, while an
/// unset one requires searching the k-mer itself. This needs no kLCP array.
///
/// When a window cannot be found, the search stopped on a substring that does not occur in the
/// superstring, so that all windows containing it are absent on this strand. To find such substrings
/// cheaply, the unresolved windows containing the position where the search stopped are searched
/// starting from the lowest one, whose search stops the earliest.
template <bool hashes, typename out_t>
void lookup_streaming(const index &ix, const char *sequence, const segment *segments, size_t count, out_t *out,
                      unsigned batch_size, lookup_stats &stats) {
    const int k = ix.k;
    const int p = ix.prefix_length;
    const uint64_t *table = ix.prefix_table.data();
    // Values of out while a segment is processed; a non-negative hash or 1 means found.
    const out_t UNRESOLVED = hashes ? -1 : 0;
    const out_t ABSENT[2] = {(out_t)(hashes ? -2 : 2), (out_t)(hashes ? -3 : 3)};
    const int NO_BLOCK = -1;

    enum : uint8_t { IDLE, START, SEARCH, CHECK };
    struct slot {
        size_t i, j;
        const char *sequence;
        /// Results of the windows of the segment.
        out_t *out;
        uint64_t key;
        int32_t length, windows;
        /// The window searched or matched, in the coordinates of the strand.
        int32_t x;
        /// Number of characters of the window already searched.
        int32_t t;
        /// The highest unresolved window if the lowest one of a block is searched, NO_BLOCK otherwise.
        int32_t block_top;
        /// Windows above this one are searched from the top only.
        int32_t top_only_above;
        int32_t hits[2];
        uint8_t state;
        /// Whether the reverse complement is searched, and whether it is the second strand,
        /// and whether the matched string has exactly k characters.
        uint8_t reverse, second, exact;
    };
    std::vector<slot> slots(std::max(1u, batch_size));
    size_t next = 0, active = slots.size();
    int score = 0;
    uint64_t found = 0, steps = 0, kmers = 0, searches = 0, rechecks = 0;

    // Character of the strand at position v, and the result of its window y.
    auto symbol = [&](const slot &s, int v) -> byte {
        return s.reverse ? 3 - nucleotideToInt[(uint8_t)s.sequence[s.length - 1 - v]]
                         : nucleotideToInt[(uint8_t)s.sequence[v]];
    };
    auto result = [&](const slot &s, int y) -> out_t & { return s.out[s.reverse ? s.windows - 1 - y : y]; };
    auto unresolved = [&](const slot &s, int y) {
        out_t value = result(s, y);
        return value == UNRESOLVED || (s.second && value == ABSENT[0]);
    };

    std::function<void(slot &)> load;

    auto begin_search = [&](slot &s, int y) {
        ++searches;
        s.x = y;
        s.t = 0;
        if (p == 0) {
            s.i = 0;
            s.j = ix.size;
            s.state = SEARCH;
            return;
        }
        uint64_t key = 0;
        for (int t = 0; t < p; ++t) key = (key << 2) | symbol(s, y + k - 1 - t);
        s.key = key;
        s.state = START;
        __builtin_prefetch(table + 2 * key);
    };
    std::function<void(slot &)> strand_done;
    // Search the highest unresolved window not above [from].
    auto next_top = [&](slot &s, int from) {
        while (from >= 0 && !unresolved(s, from)) --from;
        s.block_top = NO_BLOCK;
        if (from < 0) {
            strand_done(s);
        } else {
            begin_search(s, from);
        }
    };
    strand_done = [&](slot &s) {
        if (!s.second) {
            s.second = 1;
            s.reverse ^= 1;
            s.top_only_above = INT32_MAX;
            next_top(s, s.windows - 1);
        } else {
            for (int y = 0; y < s.windows; ++y) {
                if (s.out[y] == ABSENT[0] || s.out[y] == ABSENT[1]) s.out[y] = UNRESOLVED;
            }
            if (s.hits[0] != s.hits[1]) score = std::max(-7, std::min(7, score + (s.hits[0] > s.hits[1] ? 1 : -1)));
            load(s);
        }
    };
    load = [&](slot &s) {
        if (next == count) {
            s.state = IDLE;
            --active;
            return;
        }
        const segment &g = segments[next++];
        s.sequence = sequence + g.offset;
        s.out = out + g.offset;
        s.length = g.length;
        s.windows = g.length - k + 1;
        kmers += s.windows;
        s.reverse = score < 0;
        s.second = 0;
        s.hits[0] = s.hits[1] = 0;
        s.top_only_above = INT32_MAX;
        next_top(s, s.windows - 1);
    };
    // The search of window x stopped: the substring from position a to the end of the window does not occur.
    auto search_failed = [&](slot &s) {
        int a = s.x + k - s.t;
        if (s.block_top == NO_BLOCK) {
            result(s, s.x) = ABSENT[s.second];
            // The unresolved windows containing position a.
            int lowest = std::max(0, a - k + 1);
            int top = s.x - 1;
            while (top >= lowest && !unresolved(s, top)) --top;
            if (top < lowest) {
                next_top(s, top);
                return;
            }
            int bottom = lowest;
            while (!unresolved(s, bottom)) ++bottom;
            if (bottom == top || top > s.top_only_above) {
                begin_search(s, top);
            } else {
                s.block_top = top;
                begin_search(s, bottom);
            }
        } else {
            // All windows from x to a contain the substring.
            for (int y = s.x; y <= std::min(a, s.block_top); ++y) {
                if (unresolved(s, y)) result(s, y) = ABSENT[s.second];
            }
            int bottom = a + 1;
            while (bottom <= s.block_top && !unresolved(s, bottom)) ++bottom;
            if (bottom > s.block_top) {
                next_top(s, s.x - 1);
            } else if (bottom == s.block_top) {
                s.block_top = NO_BLOCK;
                begin_search(s, bottom);
            } else {
                begin_search(s, bottom);
            }
        }
    };
    auto to_check = [&](slot &s) {
        s.state = CHECK;
        ix.mask_prefetch(s.i);
        ix.bwt.prefetch(s.i);
        ix.bwt.prefetch(s.j);
    };
    auto after_search_step = [&](slot &s) {
        if (s.i == s.j) {
            search_failed(s);
        } else if (s.t == k) {
            s.exact = 1;
            to_check(s);
        } else {
            s.state = SEARCH;
            ix.bwt.prefetch(s.i);
            ix.bwt.prefetch(s.j);
        }
    };
    // The interval is the one of a string starting with window x: resolve the window and go to the next one.
    auto check = [&](slot &s) {
        size_t q = s.i;
        while (q < s.j && !ix.mask_get(q)) ++q;
        if (q < s.j) {
            result(s, s.x) = hashes ? (out_t)ix.mask_rank(q) : (out_t)1;
            ++found;
            ++s.hits[s.reverse];
        } else if (!s.exact) {
            // Other occurrences of the k-mer may be set in the mask.
            ++rechecks;
            s.block_top = NO_BLOCK;
            begin_search(s, s.x);
            return;
        } else {
            result(s, s.x) = ABSENT[s.second];
        }
        if (s.block_top != NO_BLOCK) {
            // The lowest window of a block occurs, so nothing is known about the windows above it.
            int top = s.block_top;
            s.top_only_above = s.x;
            s.block_top = NO_BLOCK;
            begin_search(s, top);
            return;
        }
        if (s.x == 0 || !unresolved(s, s.x - 1)) {
            next_top(s, s.x - 1);
            return;
        }
        --s.x;
        update_range(ix, s.i, s.j, symbol(s, s.x));
        ++steps;
        if (s.i == s.j) {
            begin_search(s, s.x);
        } else {
            s.exact = 0;
            to_check(s);
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
                after_search_step(s);
                break;
            case SEARCH:
                update_range(ix, s.i, s.j, symbol(s, s.x + k - 1 - s.t));
                ++s.t;
                ++steps;
                after_search_step(s);
                break;
            case CHECK:
                check(s);
                break;
            default:
                break;
            }
        }
    }

    stats.kmers += kmers;
    stats.found += found;
    stats.steps += steps;
    stats.searches += searches;
    stats.rechecks += rechecks;
}

} // namespace fmsi2

#pragma once

#include <random>
#include <vector>

#include "../src/rank2.h"
#include "gtest/gtest.h"

namespace {

TEST(RANK2, QUAD_RANK) {
    std::mt19937_64 rng(1);
    for (size_t size : {0, 1, 223, 224, 225, 4000, 57344, 57345, 130000}) {
        std::vector<uint8_t> symbols(size);
        fmsi2::quad_rank ranks(size);
        for (size_t i = 0; i < size; ++i) {
            symbols[i] = rng() % 4;
            ranks.set(i, symbols[i]);
        }
        ranks.finalize();
        size_t counts[4] = {0, 0, 0, 0};
        for (size_t i = 0; i <= size; ++i) {
            for (int c = 0; c < 4; ++c) ASSERT_EQ(ranks.rank(i, c), counts[c]) << "size " << size << " position " << i;
            for (size_t distance : {0, 1, 2, 30, 250}) {
                if (i + distance > size) continue;
                size_t rank_i, rank_j;
                uint8_t c = rng() % 4;
                ranks.rank_pair(i, i + distance, c, rank_i, rank_j);
                ASSERT_EQ(rank_i, ranks.rank(i, c));
                ASSERT_EQ(rank_j, ranks.rank(i + distance, c));
            }
            if (i < size) {
                ASSERT_EQ(ranks[i], symbols[i]);
                counts[symbols[i]]++;
            }
        }
    }
}

TEST(RANK2, BIT_RANK) {
    std::mt19937_64 rng(2);
    for (size_t size : {0, 1, 495, 496, 497, 5000, 63488, 63489, 140000}) {
        std::vector<bool> bits(size);
        fmsi2::bit_rank ranks(size);
        for (size_t i = 0; i < size; ++i) {
            bits[i] = rng() % 3 == 0;
            if (bits[i]) ranks.set(i);
        }
        ranks.finalize();
        size_t ones = 0;
        for (size_t i = 0; i <= size; ++i) {
            ASSERT_EQ(ranks.rank(i), ones) << "size " << size << " position " << i;
            if (i < size) {
                ASSERT_EQ(ranks[i], bits[i]);
                ones += bits[i];
            }
        }
    }
}

/// The packed bytes of a sequence of 2-bit symbols, written and read in pieces.
TEST(RANK2, PACKED_SYMBOLS_IN_PIECES) {
    std::mt19937_64 rng(11);
    for (size_t size : {(size_t)1, (size_t)223, (size_t)224, (size_t)225, (size_t)100003}) {
        std::vector<uint8_t> symbols(size), bytes((size + 3) / 4, 0);
        for (size_t i = 0; i < size; ++i) {
            symbols[i] = rng() % 4;
            bytes[i / 4] |= symbols[i] << (2 * (i % 4));
        }
        fmsi2::quad_rank sequence(size);
        for (size_t first = 0; first < bytes.size();) {
            const size_t count = std::min<size_t>(bytes.size() - first, 1 + rng() % 300);
            sequence.set_packed(first, bytes.data() + first, count);
            first += count;
        }
        sequence.finalize();
        size_t ranks[4] = {0, 0, 0, 0};
        for (size_t i = 0; i < size; ++i) {
            for (uint8_t c = 0; c < 4; ++c) ASSERT_EQ(sequence.rank(i, c), ranks[c]);
            ranks[symbols[i]]++;
        }
        std::vector<uint8_t> back(bytes.size(), 0xFF);
        for (size_t first = 0; first < bytes.size();) {
            const size_t count = std::min<size_t>(bytes.size() - first, 1 + rng() % 300);
            sequence.get_packed(first, back.data() + first, count);
            first += count;
        }
        ASSERT_EQ(back, bytes);
    }
}

} // namespace

#ifndef COLORGRAM_UTILS_HPP
#define COLORGRAM_UTILS_HPP

#include "config.h"
#include <iostream>
#include <fstream>
#include <sdsl/sfstream.hpp>
#include <sdsl/structure_tree.hpp>

// using spp::sparse_hash_map;
using namespace std;

static inline uint8_t symbol_to_id(const char c) {
    switch (c) {
        case '$':
            return 0;
        case 'A':
        case 'a':
            return 1;
        case 'C':
        case 'c':
            return 2;
        case 'G':
        case 'g':
            return 3;
        case 'T':
        case 't':
            return 4;
        default:
            return 255;
    }
}

static inline uint8_t symbol_to_bits(const char c) {
    switch (c) {
        case '$':
            return 0b000;
        case 'A':
            return 0b001;
        case 'C':
            return 0b011;
        case 'G':
            return 0b101;
        case 'T':
            return 0b111;
        default:
            return 255;
    }
}

static inline char bits_to_char(const uint8_t s) {
    switch (s) {
        case 0b000:
            return '$';
        case 0b001:
            return 'A';
        case 0b011:
            return 'C';
        case 0b101:
            return 'G';
        case 0b111:
            return 'T';
        default:
            return 0;
    }
}

static inline uint8_t bits_to_id(const uint8_t s) {
    switch (s) {
        case 0b000:
            return 0;
        case 0b001:
            return 1;
        case 0b011:
            return 2;
        case 0b101:
            return 3;
        case 0b111:
            return 4;
        default:
            return 255;
    }
}

static inline uint8_t id_to_bits(const uint8_t i) {
    switch (i) {
        case 0:
            return 0b000;
        case 1:
            return 0b001;
        case 2:
            return 0b011;
        case 3:
            return 0b101;
        case 4:
            return 0b111;
        default:
            return 255;
    }
}

static bool file_exists(const string& fname) {
    ifstream f(fname);
    if (!f.is_open()) {
        return false;
    }
    f.close();
    return true;
}

template<uint16_t KMERBITS>
struct compare_bit_vector {
    const std::bitset<KMERBITS> bs;
    const std::bitset<KMERBITS> mask = std::bitset<KMERBITS>(string(LOGSIGMA, '1'));
    int kmer_bits;

    explicit compare_bit_vector(int pkmer_bits) : kmer_bits(pkmer_bits) {}

    bool operator()(const std::bitset<KMERBITS>& a, const std::bitset<KMERBITS>& b) const {
        for (int i = kmer_bits - 1; i >= 0; --i) {
            if (a[i] ^ b[i]) {
                return b[i];
            }
        }
        return false;
    }

    std::bitset<KMERBITS> min_value() const { return bs; }

    std::bitset<KMERBITS> max_value() const { return mask; }
};


template<uint16_t KMERBITS>
struct compare_symbol_bit_vector {
    const std::bitset<KMERBITS> bs;
    const std::bitset<KMERBITS> mask = std::bitset<KMERBITS>(string(LOGSIGMA, '1'));
    int kmer_bits;

    explicit compare_symbol_bit_vector(int pkmer_bits) : kmer_bits(pkmer_bits) {}

    bool operator()(const std::bitset<KMERBITS>& a, const std::bitset<KMERBITS>& b) const {
        for (int i = kmer_bits - LOGSIGMA - 1; i >= 0; --i) {
            if (a[i] ^ b[i]) {
                return b[i];
            }
        }
        // return false;
        if (a[kmer_bits - 1] ^ b[kmer_bits - 1]) {
            return b[kmer_bits - 1];
        }
        if (a[kmer_bits - 2] ^ b[kmer_bits - 2]) {
            return b[kmer_bits - 2];
        }
        if (a[kmer_bits - 3] ^ b[kmer_bits - 3]) {
            return b[kmer_bits - 3];
        }
        return false;
    }

    std::bitset<KMERBITS> min_value() const { return bs; }

    std::bitset<KMERBITS> max_value() const { return mask; }
};


#endif //COLORGRAM_UTILS_HPP

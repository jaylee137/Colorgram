#ifndef COLORGRAM_DBG_H
#define COLORGRAM_DBG_H

#include <iostream>
#include <array>
#include <string>
#include <bitset>
#include <stxxl/vector>
#include <sdsl/bit_vectors.hpp>
#include <sparsepp/spp.h>
#include <stxxl.h>
#include <stack>

#include "utils.hpp"
#include "config.h"
#include "succinct_dbg.h"

using namespace std;
using namespace sdsl;
using namespace spp;


template<uint16_t KMERBITS>
class ColoredDeBrujinGraph {
public:
    ColoredDeBrujinGraph(uint32_t pk, uint8_t pcologram_type, const string& pkmer_db_fname, const string& pbegin_db_fname,
                  const string& pend_db_fname) : k(pk), kmer_bits(LOGSIGMA * pk), colorgram_type(pcologram_type) {
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            bitset<KMERBITS> sid = symbol_to_bits(base[i]);
            sid <<= kmer_bits - LOGSIGMA;
            shifted_sids[i] = sid;
        }

        build_graph(pkmer_db_fname, pbegin_db_fname, pend_db_fname);
    }

    ~ColoredDeBrujinGraph() { delete sdbg; }

    void build_colored_graph(uint32_t color, const string& dna_str);

    SuccinctDeBruijnGraph* get_sdbg() { return sdbg; }

    void print_stats();

    void sort_color_table();

    bool save_graph(const string& output_fname);

private:
    void build_graph(const string& kmer_db_fname, const string& begin_db_fname, const string& end_db_fname);

    inline bool compare_nodes_begin(const bitset<KMERBITS>& a, const bitset<KMERBITS>& b);

    inline bool compare_nodes_end(const bitset<KMERBITS>& a, const bitset<KMERBITS>& b);

    string kmer_to_str(bitset<KMERBITS> kmer_str, uint32_t k);

    inline void add_color(size_t& kmer_color_hash, size_t color);

    inline size_t add_color_class(const bitset<MAXCOLORS>& bitvector);

    inline size_t gen_next_hash_id(size_t hashv) const {
        // mt19937_64 rand_rng(hashv);
        // return rand_rng();
        return hashv + 1;
    }

    array<bitset<KMERBITS>, SIGMA + 1> shifted_sids;

    uint32_t k;
    uint32_t C = 0;
    uint32_t kmer_bits;
    uint32_t colorgram_type;

    SuccinctDeBruijnGraph *sdbg = nullptr;

    sparse_hash_map<uint64_t, uint64_t> cids;

    typedef typename stxxl::VECTOR_GENERATOR<color_class_t>::result color_vector_type;
    color_vector_type color_table;
    size_t_vector_type label_hash_vector;
    size_t_vector_type gaps;
    hash<bitset<MAXCOLORS>> hash_color_class;

    size_t set_bits = 0;
};


#endif //COLORGRAM_DBG_H

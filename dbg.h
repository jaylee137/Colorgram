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
    ColoredDeBrujinGraph(uint32_t pk, uint8_t pcologram_type, const string& pkmer_db_fname,
                         const string& pbegin_db_fname, const string& pend_db_fname,
                         const string& pout_fname) : k(pk), kmer_bits(LOGSIGMA * pk), colorgram_type(pcologram_type),
                                                     out_fname(pout_fname) {
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            bitset<KMERBITS> sid = symbol_to_bits(base[i]);
            sid <<= kmer_bits - LOGSIGMA;
            shifted_sids[i] = sid;
        }

        build_graph(pkmer_db_fname, pbegin_db_fname, pend_db_fname);
    }

    ~ColoredDeBrujinGraph() { delete sdbg; }

    void build_label_vector(size_t color, const string& dna_str, sparse_hash_map<uint64_t, uint64_t>& H,
                            sparse_hash_map<uint64_t, uint8_t>& visited);

    void sort_label_vector();

    void build_color_table(size_t color, const string& dna_str);

    SuccinctDeBruijnGraph *get_sdbg() { return sdbg; }

    void print_stats();

    void create_succinct_structures();

    bool save_graph();

private:
    void build_graph(const string& kmer_db_fname, const string& begin_db_fname, const string& end_db_fname);

    inline bool compare_nodes_begin(const bitset<KMERBITS>& a, const bitset<KMERBITS>& b);

    inline bool compare_nodes_end(const bitset<KMERBITS>& a, const bitset<KMERBITS>& b);

    string kmer_to_str(bitset<KMERBITS> kmer_str, uint32_t k);

    inline void add_color(size_t color, size_t r);

    // inline size_t add_color_class(const bitset<MAXCOLORS>& bitvector);

    inline size_t gen_next_hash_id(size_t hashv) const {
        // mt19937_64 rand_rng(hashv);
        // return rand_rng();
        return hashv + 1;
    }

    void save_multiplicities(const vector<vector<size_t>>& buckets);

    array<bitset<KMERBITS>, SIGMA + 1> shifted_sids;

    uint32_t k;
    size_t C = 0;
    uint32_t kmer_bits;
    uint32_t colorgram_type;
    string out_fname;

    SuccinctDeBruijnGraph *sdbg = nullptr;

    // the first element of the pair stores the label and the second the number of set bits
    // typedef typename stxxl::VECTOR_GENERATOR<pair<size_t, size_t>>::result label_vector_type;
    vector<pair<size_t, size_t>> X; // label_vector_type
    size_t max_label = 0;
    size_t label_vector_length = 0;

    // color_vector_class color_table;
    vector<sd_vector_builder *> color_table;
    vector<sd_vector<> *> color_table_succinct;
    vector<size_t> color_table_row_capacities;
    size_t color_table_set_bits = 0;
    // vector<size_t> label_hash_vector;
    // size_t_vector_type gaps;
    // hash<bitset<MAXCOLORS>> hash_color_class;
};


#endif //COLORGRAM_DBG_H

#ifndef COLORGRAM_SUCCINCT_DBG_H
#define COLORGRAM_SUCCINCT_DBG_H

#include <iostream>
#include <array>
#include <string>
#include <sparsepp/spp.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "utils.hpp"
#include "config.h"

using namespace std;
using namespace sdsl;
using namespace spp;


class SuccinctDeBruijnGraph {
public:
    SuccinctDeBruijnGraph(size_t pn, uint32_t pk, uint8_t ptype, const array<size_t, SIGMA>& pT,
                          const bit_vector& pBL, const bit_vector& pBF,
                          const string& tmp_fname) : n(pn), k(pk), T(pT), BL(sd_vector<>(pBL)), BF(sd_vector<>(pBF)),
                                                     BL_rank(&BL), BL_select(&BL), BF_rank(&BF), BF_select(&BF),
                                                     SBV(calculate_SBV(pBL, ptype)), SBV_rank(&SBV), SBV_select(&SBV),
                                                     label_vect_size(SBV_rank.rank(SBV.size())) {

        // number of nodes
        v = BL_rank.rank(BL.size());

        // construct the edges wavelet tree
        construct(edges, tmp_fname, 1);

        // construct the edges static vector for faster scanning purposes...
        // construct_edges_static();

        // calculate T_F that stores the starting positions of the symbols in F
        calc_F_L_node_cnt();
    }

    SuccinctDeBruijnGraph(const string& input_fname) {
        load(input_fname);

        BL_rank = sd_vector<>::rank_1_type(&BL);
        BL_select = sd_vector<>::select_1_type(&BL);
        BF_rank = sd_vector<>::rank_1_type(&BF);
        BF_select = sd_vector<>::select_1_type(&BF);
        X_select = sd_vector<>::select_1_type(&X);
        CT_rank = sd_vector<>::rank_1_type(&CT);
        CT_select = sd_vector<>::select_1_type(&CT);

        SBV = sd_vector<>(calculate_SBV());
        SBV_rank = sd_vector<>::rank_1_type(&SBV);
        SBV_select = sd_vector<>::select_1_type(&SBV);
        label_vect_size = SBV_rank.rank(SBV.size());

        calc_F_L_node_cnt();
    }

    inline uint8_t indegree(size_t index);

    size_t forward(size_t index, uint8_t c) const;

    inline size_t backward(size_t index);

    bitset<MAXCOLORS> get_color_class(size_t index);

    size_t get_next_symbol_index(size_t index, uint8_t c) const;

    size_t get_label_index(size_t index) const;

    size_t get_label_vect_size() const;

    size_t get_num_of_edges() const { return n; }

    size_t get_num_of_nodes() const { return v; }

    void set_C(uint32_t pC) { C = pC; }

    void set_X(sd_vector_builder *& vector_builder) {
        // X = select_support_mcl2(label_hash_vector, label_permutation, cids);
        X = sd_vector<>(*vector_builder);
        X_select = sd_vector<>::select_1_type(&X);
    }

    void set_CT(sd_vector_builder *& vector_builder) {
        CT = sd_vector<>(*vector_builder);
        CT_rank = sd_vector<>::rank_1_type(&CT);
        CT_select = sd_vector<>::select_1_type(&CT);
    }

    bool save(const string& output_fname) const;

    void print_stats(ostream& out);

private:
    // void construct_edges_static();

    bit_vector calculate_SBV(const bit_vector& pBL, uint8_t ptype = 0);

    bit_vector calculate_SBV(uint8_t ptype = 0);

    void calc_F_L_node_cnt() {
        // calculate F_node_cnt and L_node_cnt
        size_t fsum = 0;
        for (uint8_t i = 0; i < SIGMA; ++i) {
            T_F[i] = fsum;
            fsum += edges.rank(edges.size(), id_to_bits(i + 1));
            F_node_cnt[i] = BF_rank.rank(T_F[i]);
            L_node_cnt[i] = BL_rank.rank(T[i]);
        }
    }

    void update_color_class(size_t index, bitset<MAXCOLORS>& color_class);

    size_t save_dbg(ostream& out, structure_tree_node *v = NULL, string name = "") const;

    size_t save_label_vect(ostream& out, structure_tree_node *v = NULL, string name = "") const;

    size_t save_color_table(ostream& out, structure_tree_node *v = NULL, string name = "") const;

    size_t save_storage_vect(ostream& out, structure_tree_node *v = NULL, string name = "") const;

    void load_dbg(istream& in);

    void load_label_vect(istream& in);

    void load_color_table(istream& in);

    void load_storage_vect(istream& in);

    bool load(const string& input_fname);


    size_t n;
    size_t v;
    uint32_t k;
    array<size_t, SIGMA> T{};
    array<size_t, SIGMA> T_F{};
    sd_vector<> BL;
    sd_vector<> BF;
    sd_vector<>::rank_1_type BL_rank;
    sd_vector<>::select_1_type BL_select;
    sd_vector<>::rank_1_type BF_rank;
    sd_vector<>::select_1_type BF_select;
    array<size_t, SIGMA> L_node_cnt{};
    array<size_t, SIGMA> F_node_cnt{};
    typedef wt_huff<rrr_vector<63>> wt_t;
    wt_t edges;
    typedef typename stxxl::VECTOR_GENERATOR<uint8_t>::result uint8_t_vector_type;
    uint8_t_vector_type edges_static;
    sd_vector<> SBV;
    sd_vector<>::rank_1_type SBV_rank;
    sd_vector<>::select_1_type SBV_select;
    size_t label_vect_size;


    sd_vector<> X;
    sd_vector<>::select_1_type X_select;
    uint32_t C;
    sd_vector<> CT;
    sd_vector<>::rank_1_type CT_rank;
    sd_vector<>::select_1_type CT_select;
};


#endif //COLORGRAM_SUCCINCT_DBG_H

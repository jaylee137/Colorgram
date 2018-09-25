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

        // calculate T_F that stores the starting positions of the symbols in F
        // calculate F_node_cnt and L_node_cnt
        size_t fsum = 0;
        for (uint8_t i = 0; i < SIGMA; ++i) {
            T_F[i] = fsum;
            fsum += edges.rank(edges.size(), id_to_bits(i + 1));
            F_node_cnt[i] = BF_rank.rank(T_F[i]);
            L_node_cnt[i] = BL_rank.rank(T[i]);
        }
    }

    uint8_t indegree(size_t index);

    size_t forward(size_t index, uint8_t c) const;

    // size_t backward(size_t index, uint8_t c);

    // size_t backward(size_t index);

    size_t get_next_symbol_index(size_t index, uint8_t c) const;

    size_t get_label_index(size_t index) const;

    size_t get_label_vect_size() const;

    size_t get_num_of_edges() const { return n; }

    size_t get_num_of_nodes() const { return v; }

    void set_X(bit_vector*& bv) { X = select_support_mcl<1, 1>(bv); }

    void set_CT(sd_vector_builder*& vector_builder) {
        CT = sd_vector<>(*vector_builder);
        CT_rank = sd_vector<>::rank_1_type(&CT);
        CT_select = sd_vector<>::select_1_type(&CT);
    }

    size_t serialize_dbg(ostream& out, structure_tree_node *v = NULL, string name = "") const;

    size_t serialize_label_vect(ostream& out, structure_tree_node *v = NULL, string name = "") const;

    size_t serialize_color_table(ostream& out, structure_tree_node *v = NULL, string name = "") const;

    size_t serialize_storage_vect(ostream& out, structure_tree_node *v = NULL, string name = "") const;

private:
    bit_vector calculate_SBV(const bit_vector& pBL, uint8_t ptype);

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
    sd_vector<> SBV;
    sd_vector<>::rank_1_type SBV_rank;
    sd_vector<>::select_1_type SBV_select;
    size_t label_vect_size;


    select_support_mcl<1, 1> X;
    sd_vector<> CT;
    sd_vector<>::rank_1_type CT_rank;
    sd_vector<>::select_1_type CT_select;
};


#endif //COLORGRAM_SUCCINCT_DBG_H

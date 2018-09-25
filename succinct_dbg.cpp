#include "succinct_dbg.h"


bit_vector SuccinctDeBruijnGraph::calculate_SBV(const bit_vector& pBL, uint8_t ptype = 0) {
    bit_vector rv(pBL.size());
    bool prev_zero = false;
    for (size_t i = 0; i < pBL.size(); ++i) {
        if (prev_zero) {
            rv[i] = true;
            prev_zero = false;
        }
        else if (pBL[i] == 0) {
            rv[i] = true;
            prev_zero = true;
        }
        else if (ptype == 1 && edges[i] != 0 && indegree(forward(i, edges[i])) > 0) {

        }
        else if (ptype == 2) {

        }
        else {
            rv[i] = false;
        }
    }
    return rv;
}


uint8_t SuccinctDeBruijnGraph::indegree(size_t index) {
    return 0;
}


// returns the node index (first outgoing edge index) when we read symbol c on the given position
size_t SuccinctDeBruijnGraph::forward(size_t index, uint8_t c) const {
    assert(c != 0);

    // uint8_t ac = edges[index];
    auto cid = (uint8_t) (bits_to_id(c) - 1);
    size_t j = edges.rank(index, c);
    size_t l = T_F[cid] + j;
    size_t node_cnt = BF_rank.rank(l) - F_node_cnt[cid];
    return BL_select.select(L_node_cnt[cid] + node_cnt) + 1;
}


size_t SuccinctDeBruijnGraph::get_next_symbol_index(size_t index, uint8_t c) const {
    assert(c != 0);
    return edges.select(edges.rank(index, c) + 1, c);
    // for (size_t i = index; i < edges.size(); ++i) {
    //     if (edges[i] == c) {
    //         return i;
    //     }
    // }
    //
    // // in normal circumstances this could never happen...
    // return edges.size();
}


size_t SuccinctDeBruijnGraph::get_label_index(size_t index) const {
    if (SBV[index] == 0) {
        return label_vect_size;
    }

    return SBV_rank.rank(index);
}


// getter for label_vect_size
size_t SuccinctDeBruijnGraph::get_label_vect_size() const {
    return label_vect_size;
}


size_t SuccinctDeBruijnGraph::serialize_dbg(ostream& out, structure_tree_node *v, string name) const {
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += write_member(n, out, child, "n");
    written_bytes += write_member(k, out, child, "k");
    for (size_t i = 0; i < SIGMA; ++i) {
        written_bytes += write_member(T[i], out, child, "T[" + to_string(i) +"]");
    }
    written_bytes += BL.serialize(out, child, "BL");
    written_bytes += BL_rank.serialize(out, child, "BL_rank");
    written_bytes += BL_select.serialize(out, child, "BL_select");
    written_bytes += BF.serialize(out, child, "BF");
    written_bytes += BF_rank.serialize(out, child, "BF_rank");
    written_bytes += BF_select.serialize(out, child, "BF_select");
    written_bytes += edges.serialize(out, child, "edges");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}


size_t SuccinctDeBruijnGraph::serialize_label_vect(ostream& out, structure_tree_node *v, string name) const {
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += X.serialize(out, child, "X");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}


size_t SuccinctDeBruijnGraph::serialize_color_table(ostream& out, structure_tree_node *v, string name) const {
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += CT.serialize(out, child, "CT");
    written_bytes += CT_rank.serialize(out, child, "CT_rank");
    written_bytes += CT_select.serialize(out, child, "CT_select");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}


size_t SuccinctDeBruijnGraph::serialize_storage_vect(ostream& out, structure_tree_node *v, string name) const {
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += SBV.serialize(out, child, "SBV");
    written_bytes += SBV_rank.serialize(out, child, "SBV_rank");
    written_bytes += SBV_select.serialize(out, child, "SBV_select");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}


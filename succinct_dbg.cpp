#include "succinct_dbg.h"
#include <ctime>


void SuccinctDeBruijnGraph::calc_F_L_node_cnt() {
    // calculate F_node_cnt and L_node_cnt
    size_t fsum = 0;
    for (uint8_t i = 0; i < SIGMA; ++i) {
        T_F[i] = fsum;
        fsum += edges.rank(edges.size(), id_to_bits((uint8_t)(i + 1)));
        F_node_cnt[i] = BF_rank.rank(T_F[i]);
        L_node_cnt[i] = BL_rank.rank(T[i]);
    }
}


bit_vector SuccinctDeBruijnGraph::calculate_SBV(const bit_vector& pBL, uint8_t ptype) {
    cerr << "Generating Storage Bit Vector (SBV)..." << endl;
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
        // else if (pBL[i] == 1 && indegree(i) > 1) {
        //     rv[i] = true;
        // }
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


bit_vector SuccinctDeBruijnGraph::calculate_SBV(uint8_t ptype) {
    cerr << "Generating Storage Bit Vector (SBV)..." << endl;
    cerr << "logn: " << (size_t(log2(this->n)) + 1) << endl;
    bit_vector rv(BL.size());
    sparse_hash_map<size_t, uint8_t> visited;
    size_t goodcnt = 0;
    auto sample = [this, &visited, &rv, &goodcnt](size_t i) {
        const size_t logn = max(size_t(log2(n)) + 1, (size_t)20);
        size_t cnt = 0;
        do {
            visited[i] = 1;
            if (++cnt >= logn) {
                cnt = 0;
                rv[i] = true;
                goodcnt++;
            }
            uint8_t c = edges[i];
            if (c == 0) {
                break;
            }
            i = forward(i, c);
        } while (visited.find(i) == visited.end() && BL[i] == 1 && indegree(i) == 1);
    };

    bool prev_zero = false;
    for (size_t i = 0; i < BL.size(); ++i) {
        if (prev_zero) {
            // sample(i);
            rv[i] = true;
            prev_zero = false;
            goodcnt++;
        }
        else if (BL[i] == 0) {
            // sample(i);
            rv[i] = true;
            prev_zero = true;
            goodcnt++;
        }
        else if (BL[i] == 1 && indegree(i) > 1) {
            // rv[i] = true;
            Ep1++;
        }
        else if (ptype == 1 && edges[i] != 0 && indegree(forward(i, edges[i])) > 0) {

        }
        else if (ptype == 2) {

        }
        else {
            rv[i] = false;
        }
    }
    cerr << "Label vector size:" << goodcnt << endl;
    return rv;
}


/// Returns the in-degree of the given node corresponding to the edge_index
/// \param index - edge index of W(G)
/// \return in-degree
uint8_t SuccinctDeBruijnGraph::indegree(size_t index) {
    size_t _v = BL_rank.rank(index);
    if (_v == 0) {
        return 0;
    }
    size_t q = BF_select.select(_v);
    size_t r = (_v > 1) ? q - BF_select.select(_v - 1) : q + 1;
    return (uint8_t) (r);
}


/// Returns the out-degree of the given node corresponding to the edge_index
/// \param index - edge index of W(G)
/// \return out-degree
uint8_t SuccinctDeBruijnGraph::outdegree(size_t index) {
    size_t r = BL_rank.rank(index);
    return (r == 0) ? (uint8_t) (BL_select.select(1) + 1) : (uint8_t) (BL_select.select(r + 1) - BL_select.select(r));
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


/// Updates a given color_class bit vector with the given row of the Color Table
/// \param index - row index of the Color Table
/// \param color_class
/// \param num_of_colors - increments each newly introduced color
/// \return number of operations
size_t SuccinctDeBruijnGraph::update_color_class(size_t index, bitset<MAXCOLORS>& color_class, size_t& num_of_colors) {
    size_t number_of_operations = 0;
    if (C > 100) {
        // for bigger number of colors use select operations to navigate through the rows of the color table...
        size_t start_index = index * C;
        size_t i = CT_rank.rank(start_index) + 1;
        index = CT_select.select(i);
        number_of_operations = 1;
        for (size_t mc = start_index + C; index < mc; index = CT_select.select(++i)) {
            ++number_of_operations;
            if (!color_class[index - start_index]) {
                ++num_of_colors;
                color_class[index - start_index] = true;
            }
            // if we hit the last set bit in CT => break (otherwise the last select operation would fail...)
            if (index >= CT_last_set_bit_position) {
                break;
            }
        }
    }
    else {
        // for small number of colors => scan through the current line of the color table
        for (size_t i = index * C, mc = (index + 1) * C, j = 0; i < mc; ++i, ++j) {
            if (CT[i]) {
                // if we introduce a new color to the color class => increment the number of used colors
                if (!color_class[j]) {
                    ++num_of_colors;
                    color_class[j] = true;
                }
            }
        }
        number_of_operations = C;
    }
    // for (size_t i, ri = CT_rank.rank(index * C + 1), mc = (index + 1) * C;
    //      i < mc && i < CT.size() - 1; i = CT_select.select(++ri)) {
    //     CT.
    //     i = CT_select.select()
    //     color_class[i] = 1;
    // }
    return number_of_operations;
}


/// Jumps backward in W(G) - !!!Note that here the given index corresponds to BF!!!
/// \param index - edge index in W(G)
/// \return the edge index of the node in W(G) that preceeds the given index
inline size_t SuccinctDeBruijnGraph::backward(size_t index) {
    // find character c
    uint8_t cid = 0;
    for (cid = SIGMA - 1; cid > 0; --cid) {
        if (T_F[cid] <= index) {
            break;
        }
    }
    return edges.select(index - T_F[cid] + 1, id_to_bits((uint8_t)(cid + 1)));
}


/// Returns the color class of the given edge
/// \param index - edge index in W(G)
/// \return color class
tuple<size_t, size_t, size_t> SuccinctDeBruijnGraph::get_color_class(bitset<MAXCOLORS>& color_class, size_t index) {
    // bitset<MAXCOLORS> color_class;
    // mark the visited nodes...
    sparse_hash_map<size_t, uint8_t> visited;
    // size_t r = BL_rank.rank(index);
    // size_t startnode_index = (r == 0) ? 0 : BL_select.select(r) + 1;
    // the number of colors appearing in the current color class
    size_t num_of_operations = 0;
    size_t num_of_colors = 0;
    stack<size_t> s;
    s.push(index);
    // visited[index] = 1;
    size_t tree_size = 0, merge_cnt = 0, leaf_cnt = 0;
    while (!s.empty()) {
        size_t edge_index = s.top();
        ++tree_size;
        s.pop();
        // size_t ci = get_label(edge_index);
        // num_of_operations += update_color_class(ci, color_class, num_of_colors);
        // break;

        // if the actual edge was saved
        if (SBV[edge_index] == 1) {
            ++leaf_cnt;
            // get the label index
            size_t xi = SBV_rank.rank(edge_index);
            // get the row index of the color table
            size_t ci = get_label(xi);
            num_of_operations += update_color_class(ci, color_class, num_of_colors);
            // if all the colors appearing in the color class => return
            if (num_of_colors >= C) {
                // exit from the while loop and return
                break;
            }
        }
        else {
            // get the incoming edges
            size_t _v = BL_rank.rank(edge_index);
            // if we hit the first node => the color class contain all the colors...
            if (_v == 0) {
                ++leaf_cnt;
                // color_class = bitset<MAXCOLORS>(string(C, '1'));
                for (size_t i = 0; i < C; ++i) {
                    color_class.set(i);
                }

                // exit from the while loop and return
                break;
            }
            size_t q = BF_select.select(_v);
            size_t r = (_v > 1) ? q - BF_select.select(_v - 1) : q + 1;
            if (r > 1) {
                ++merge_cnt;
            }
            for (size_t _q = q - r + 1; _q <= q; ++_q) {
                // step backward - push to the stack each incoming edges
                if (r == 1 || visited.find(_q) == visited.end()) {
                    s.push(backward(_q));
                    if (r > 1) {
                        visited[_q] = 1;
                    }
                }
            }
        }
    }
    return tuple<size_t, size_t, size_t>(tree_size, num_of_operations, leaf_cnt);
}


/// Returns the next edge tagged by the given character
/// \param index
/// \param c - character, according to the bit order (not char type)
/// \return the index of the next character
size_t SuccinctDeBruijnGraph::get_next_symbol_index(size_t index, uint8_t c) const {
    // static wt_t::const_iterator start_it = edges.begin();
    // static uint8_t first_char = symbol_to_bits('A');

    if (c == 0) {
        return index;
    }
    return edges.select(edges.rank(index, c) + 1, c);
    // for (size_t i = index; i < edges.size(); ++i) {
    //     if (edges[i] == c) {
    //         return i;
    //     }
    // }
    // return edges.size();

}


size_t SuccinctDeBruijnGraph::get_label_index(size_t index) const {
    if (SBV[index] == 0) {                 // UPDATE using sparse_hash
        return label_vect_size;
    }

    return SBV_rank.rank(index);
}


size_t SuccinctDeBruijnGraph::save_dbg(ostream& out, structure_tree_node *v, string name) const {
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += write_member(n, out, child, "n");
    written_bytes += write_member(this->v, out, child, "v");
    written_bytes += write_member(k, out, child, "k");
    for (size_t i = 0; i < SIGMA; ++i) {
        written_bytes += write_member(T[i], out, child, "T[" + to_string(i) + "]");
    }
    written_bytes += BL.serialize(out, child, "BL");
    written_bytes += BF.serialize(out, child, "BF");
    written_bytes += edges.serialize(out, child, "edges");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}


void SuccinctDeBruijnGraph::load_dbg(istream& in) {
    read_member(n, in);
    read_member(v, in);
    read_member(k, in);
    for (size_t i = 0; i < SIGMA; ++i) {
        read_member(T[i], in);
    }
    BL.load(in);
    BF.load(in);
    edges.load(in);
}


size_t SuccinctDeBruijnGraph::save_label_vect(ostream& out, structure_tree_node *v, string name) const {
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += X.serialize(out, child, "X");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}


void SuccinctDeBruijnGraph::load_label_vect(istream& in) {
    X.load(in);
}


size_t SuccinctDeBruijnGraph::save_color_table(ostream& out, structure_tree_node *v, string name) const {
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += write_member(C, out, child, "C");
    written_bytes += CT.serialize(out, child, "CT");
    // written_bytes += CT_rank.serialize(out, child, "CT_rank");
    // written_bytes += CT_select.serialize(out, child, "CT_select");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}


void SuccinctDeBruijnGraph::load_color_table(istream& in) {
    read_member(C, in);
    CT.load(in);
}


size_t SuccinctDeBruijnGraph::save_storage_vect(ostream& out, structure_tree_node *v, string name) const {
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += SBV.serialize(out, child, "SBV");
    // written_bytes += SBV_rank.serialize(out, child, "SBV_rank");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}


void SuccinctDeBruijnGraph::load_storage_vect(istream& in) {
    SBV.load(in);
}


bool SuccinctDeBruijnGraph::save(const string& output_fname) const {
    auto save_to_bin_file = [](string fname, auto save_fn) {
        sdsl::osfstream out(fname, std::ios::binary | std::ios::trunc | std::ios::out);
        if (!out) {
            if (sdsl::util::verbose) {
                std::cerr << "ERROR: store_to_file not successful for: '" << fname << "'" << std::endl;
            }
            return false;
        }
        save_fn(out);
        out.close();
        if (util::verbose) {
            std::cerr << "INFO: store_to_file: `" << fname << "`" << std::endl;
        }

        return true;
    };

    cerr << "Saving DBG..." << endl;
    if (!save_to_bin_file(output_fname + ".dbg", [this](ostream& out) { save_dbg(out); }))
        return false;

    cerr << "Saving Label Vector..." << endl;
    if (!save_to_bin_file(output_fname + ".x", [this](ostream& out) { save_label_vect(out); }))
        return false;

    cerr << "Saving Color Table..." << endl;
    if (!save_to_bin_file(output_fname + ".ct", [this](ostream& out) { save_color_table(out); }))
        return false;

    // cerr << "Saving Storage Vector..." << endl;
    // if (!save_to_bin_file(out_fname + ".sbv", [this](ostream& out) { save_storage_vect(out); }))
    //     return false;

    return true;
}


bool SuccinctDeBruijnGraph::load(const string& input_fname) {
    auto load_from_bin_file = [](string fname, auto load_fn) {
        isfstream in(fname, std::ios::binary | std::ios::in);
        if (!in) {
            if (util::verbose) {
                std::cerr << "Could not load file `" << fname << "`" << std::endl;
            }
            return false;
        }
        load_fn(in);
        in.close();
        if (util::verbose) {
            std::cerr << "Load file `" << fname << "`" << std::endl;
        }
        return true;
    };

    cerr << "Loading DBG..." << endl;
    if (!load_from_bin_file(input_fname + ".dbg", [this](istream& in) { load_dbg(in); }))
        return false;

    cerr << "Loading Label Vector..." << endl;
    if (!load_from_bin_file(input_fname + ".x", [this](istream& in) { load_label_vect(in); }))
        return false;

    cerr << "Loading Color Table..." << endl;
    if (!load_from_bin_file(input_fname + ".ct", [this](istream& in) { load_color_table(in); }))
        return false;

    cerr << "Loading Storage Vector..." << endl;
    if (!load_from_bin_file(input_fname + ".sbv", [this](istream& in) { load_storage_vect(in); }))
        return false;

    return true;
}


void SuccinctDeBruijnGraph::print_stats(ostream& out) {
    // ../datasets/human_transcript_orig/human_transcript_smpl
    // ../datasets/hx500k/hx500k_smpl
    // ../datasets/ecoli569/ecoli569_smpl
    // ../datasets/plant/plant_smpl

    out << "Number of edges:               " << n << endl;
    out << "Number of nodes:               " << v << endl;
    out << "k-mer size:                    " << k << endl;
    out << "Number of colors:              " << C << endl;
    out << "#Rows of Color Table (CT):     " << CT.size() / C << endl;
    out << "Size of original CT in bits:   " << CT.size() << endl;
    out << "Size of CT in bytes:           " << size_in_bytes(CT) << endl;
    out << "Length of Label Vector:        " << label_vect_size << endl;
    out << "Size of Label Vector in bytes: " << size_in_bytes(X) << endl;
    out << "# E_+^1:                       " << Ep1 << endl;
    out << std::fixed;
    out << "Ratio of 1s / (size of CT):    " << std::setprecision(6) << CT_rank.rank(CT.size()) / (double) CT.size()
        << endl;
    out << "AVG 1s in CT (row):            " << std::setprecision(6) << CT_rank.rank(CT.size()) / (double) CT.size() * C
        << endl;
    out << "AVG 1s in CT (column):         " << std::setprecision(6) << CT_rank.rank(CT.size()) / (double) C << endl;
    out << endl;

    size_t experiment_cnt = 100000;
    vector<bitset<MAXCOLORS>> color_classes(experiment_cnt);
    cerr << "experiment started..." << endl;
    size_t sum_tree_size = 0, sum_merge_cnt = 0, sum_leaf_cnt = 0;

    std::mt19937_64 rng(0);
    std::uniform_int_distribution<uint64_t> distribution(0, edges.size() - 1);
    auto dice = bind(distribution, rng);

    clock_t begin = clock();
    // struct timespec start, finish;
    // uint64_t elapsed;
    // clock_gettime(CLOCK_MONOTONIC, &start);


    for (size_t i = 0; i < experiment_cnt; ++i) {
        size_t index = dice(); // rand() % edges.size();
        auto r = get_color_class(color_classes[i], index);
        // cout << get<1>(r) << " " << get<2>(r) << endl;
        sum_tree_size += get<0>(r);
        sum_merge_cnt += get<1>(r);
        sum_leaf_cnt += get<2>(r);
        if (i % 10000 == 0) {
            out << i << " ";
        }
    }

    clock_t end = clock();
    double elapsed = double(end - begin) / CLOCKS_PER_SEC;
    // clock_gettime(CLOCK_MONOTONIC, &finish);
    // elapsed = (finish.tv_sec - start.tv_sec) * 1000000000;
    // elapsed += (finish.tv_nsec - start.tv_nsec);
    out << "Elapsed secs:                  " << elapsed << endl;
    out << "AVG time / get_color_class()   " << elapsed / (double) experiment_cnt << endl;
    out << "AVG tree size                  " << sum_tree_size / (double) experiment_cnt << endl;
    out << "AVG merge cnt                  " << sum_merge_cnt / (double) experiment_cnt << endl;
    out << "AVG leaf cnt                   " << sum_leaf_cnt / (double) experiment_cnt << endl;
}

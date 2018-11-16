#include "dbg.h"

#include <numeric>
#include <fstream>
#include <kmc_api/kmc_file.h>


template<uint16_t KMERBITS>
inline bool ColoredDeBrujinGraph<KMERBITS>::compare_nodes_begin(const bitset<KMERBITS>& a, const bitset<KMERBITS>& b) {
    for (int i = kmer_bits - LOGSIGMA - 1; i >= 0; --i) {
        if (a[i] ^ b[i]) {
            return false;
        }
    }
    return true;
}

template<uint16_t KMERBITS>
inline bool ColoredDeBrujinGraph<KMERBITS>::compare_nodes_end(const bitset<KMERBITS>& a, const bitset<KMERBITS>& b) {
    for (uint32_t i = kmer_bits; i >= LOGSIGMA; --i) {
        if (a[i] ^ b[i]) {
            return false;
        }
    }
    return true;
}

template<uint16_t KMERBITS>
string ColoredDeBrujinGraph<KMERBITS>::kmer_to_str(bitset<KMERBITS> kmer_str, uint32_t k) {
    static const bitset<KMERBITS> mask(string(LOGSIGMA, '1'));
    stringstream ss;
    for (uint32_t i = 0; i < k; ++i) {
        uint8_t ac = (uint8_t) (kmer_str & mask).to_ulong();
        ss << bits_to_char(ac);
        kmer_str >>= LOGSIGMA;
    }
    return ss.str();
}

template<uint16_t KMERBITS>
void ColoredDeBrujinGraph<KMERBITS>::build_graph(const string& kmer_db_fname, const string& begin_db_fname,
                                                 const string& end_db_fname) {
    uint32 _mode;
    uint32 _counter_size;
    uint32 _lut_prefix_length;
    uint32 _signature_len;
    uint32 _min_count;
    uint64 _max_count;
    uint64 _total_kmers;

    uint64 counter = 0;

    cerr << "Processing begin (k-1)mers... " << endl;

    CKMCFile kmer_data_base_prefixes;
    uint32_t _k = k - 1;
    if (!kmer_data_base_prefixes.OpenForListing(begin_db_fname)) {
        cerr << "Could not open database... " << begin_db_fname << endl;
        exit(1);
    }

    kmer_data_base_prefixes.Info(_k, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count,
                                 _total_kmers);

    cerr << "Total number of begin kmers: " << _total_kmers << endl;

    CKmerAPI _kmer_object(_k);

    typedef typename stxxl::VECTOR_GENERATOR<bitset<KMERBITS>>::result kmer_vector_type;
    // append $$...$ to the beginning of the starting kmcs...
    size_t ii = 0;

    array<size_t, SIGMA + 1> last_symbol_cnt{};
    cerr << "Loading begin (k-1)mers..." << endl;
    kmer_vector_type kmers_begin;
    kmer_vector_type kmers_begin2;
    while (kmer_data_base_prefixes.ReadNextKmer(_kmer_object, counter)) {
        bitset<KMERBITS> akmer;
        for (uint32_t j = 0; j < _k; ++j) {
            uint8_t sid = symbol_to_id(_kmer_object.get_asci_symbol(j));
            akmer >>= LOGSIGMA;
            akmer |= shifted_sids[sid];
        }
        akmer >>= LOGSIGMA;
        kmers_begin.push_back(akmer);

        if (++ii % 1000000 == 0) {
            cerr << ii << " ";
        }
    }
    kmer_data_base_prefixes.Close();

    kmer_vector_type kmers;


    array<size_t, SIGMA> T{};
    const bitset<KMERBITS> mask_kmer(string(kmer_bits, '1'));
    const bitset<KMERBITS> mask_last_char(string(LOGSIGMA, '1'));
    for (size_t i = 1; i <= _k; ++i) {
        cerr << i << endl;
        kmer_vector_type& km1 = kmers_begin;
        kmer_vector_type& km2 = kmers_begin2;
        ii = 0;
        for (auto akmer : km1) {
            akmer <<= LOGSIGMA;
            akmer &= mask_kmer;


            km2.push_back(akmer);

            if (++ii % 1000000 == 0) {
                cerr << ii << " ";
            }
        }

        km1.clear();

        cerr << "Sorting begin (k-1)mers... " << " kmer sizes: " << km2.size() << endl;
        stxxl::sort(km2.begin(), km2.end(), compare_bit_vector<KMERBITS>(kmer_bits), 64 * 1024ULL * 1024ULL * 1024ULL);

        cerr << "Filtering begin (k-1)mers..." << endl;
        kmers.push_back(km2[0]);
        km1.push_back(km2[0]);

        // save the last but one symbols of the begin kmers
        if (i < _k) {
            // cerr << kmer_to_str(km1.back(), k) << endl;
            auto ac = (uint8_t) (km1.back() >> (LOGSIGMA * (_k - 1)) & mask_last_char).to_ulong();
            uint8_t aid = bits_to_id(ac) - 1;
            T[aid]++;
        }
        // update the last character
        auto ac = (uint8_t) (km1.back() >> (LOGSIGMA * _k) & mask_last_char).to_ulong();
        uint8_t aid = bits_to_id(ac);
        last_symbol_cnt[aid]++;
        for (size_t j = 1; j < km2.size(); ++j) {
            if (km2[j] != km2[j - 1]) {
                kmers.push_back(km2[j]);
                km1.push_back(km2[j]);

                if (i < _k) {
                    // cerr << kmer_to_str(km1.back(), k) << endl;
                    auto ac = (uint8_t) (km1.back() >> (LOGSIGMA * (_k - 1)) & mask_last_char).to_ulong();
                    uint8_t aid = bits_to_id(ac) - 1;
                    T[aid]++;
                }
                // update the last character
                auto ac = (uint8_t) (km1.back() >> (LOGSIGMA * _k) & mask_last_char).to_ulong();
                uint8_t aid = bits_to_id(ac);
                last_symbol_cnt[aid]++;
            }
        }

        km2.clear();
    }

    kmers_begin.clear();
    kmers_begin2.clear();

    cerr << "Processing end (k-1)mers... " << endl;
    CKMCFile kmer_data_base_suffixes;
    if (!kmer_data_base_suffixes.OpenForListing(end_db_fname)) {
        cerr << "Could not open database... " << end_db_fname << endl;
        exit(1);
    }

    kmer_data_base_suffixes.Info(_k, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count,
                                 _total_kmers);

    cerr << "Total number of end kmers: " << _total_kmers << endl;

    while (kmer_data_base_suffixes.ReadNextKmer(_kmer_object, counter)) {
        bitset<KMERBITS> akmer;
        for (uint32_t j = 0; j < _k; ++j) {
            uint8_t sid = symbol_to_id(_kmer_object.get_asci_symbol(j));
            if (j == _k - 1) {
                T[sid - 1]++;
            }
            akmer >>= LOGSIGMA;
            akmer |= shifted_sids[sid];
        }
        // add an extra $ to the end
        akmer >>= LOGSIGMA;
        kmers.push_back(akmer);
    }

    kmer_data_base_suffixes.Close();


    cerr << "Processing kmers... " << endl;
    CKMCFile kmer_data_base;
    if (!kmer_data_base.OpenForListing(kmer_db_fname)) {
        cerr << "Could not open database... " << kmer_db_fname << endl;
        exit(1);
    }

    kmer_data_base.Info(k, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count,
                        _total_kmers);

    cerr << "Total number of kmers: " << _total_kmers << endl;

    CKmerAPI kmer_object(k);
    ii = 0;
    while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
        bitset<KMERBITS> akmer;
        for (uint32_t j = 0; j < k; ++j) {
            uint8_t sid = symbol_to_id(kmer_object.get_asci_symbol(j));
            if (j == k - 1) {
                last_symbol_cnt[sid]++;
            }
            else if (j == k - 2) {
                T[sid - 1]++;
            }
            akmer >>= LOGSIGMA;
            akmer |= shifted_sids[sid];
        }
        kmers.push_back(akmer);

        if (++ii % 1000000 == 0) {
            cerr << ii << "/" << _total_kmers << " ";
        }
    }
    cerr << endl;

    kmer_data_base.Close();

    cerr << "Number of edges: " << kmers.size() << endl;

    cerr << "Sorting kmers... " << endl;

    stxxl::sort(kmers.begin(), kmers.end(), compare_symbol_bit_vector<KMERBITS>(kmer_bits), 64 * 1024UL * 1024UL *
                                                                                            1024UL);

    cerr << "Creating BL, BF, edges..." << endl;
    bit_vector BL(kmers.size());

    // array<size_t, SIGMA> T { };
    array<size_t, SIGMA> char_indexes{};
    array<size_t, SIGMA> prev_indexes{};
    prev_indexes.fill(-1ULL);
    size_t border_sum = 0;
    for (size_t i = 0; i < SIGMA; ++i) {
        char_indexes[i] = border_sum;
        border_sum += last_symbol_cnt[i + 1];
        // T[i] = border_sum;
        // cout << last_symbol_cnt[i + 1] << " " << border_sum << endl;
    }

    // create temporary file for edges - wavelet tree
    string tmp_fname = "edges.temp";
    stxxl::syscall_file tmp_edge_file(tmp_fname,
                                      stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT |
                                      stxxl::file::TRUNC);
    tmp_edge_file.set_size(kmers.size());
    auto W = new stxxl::vector<uint8_t>(&tmp_edge_file);
    auto w_bufwriter = new stxxl::vector<uint8_t>::bufwriter_type(*W);

    // create static vector for the edge tags
    // uint8_t_vector_type edges_static(kmers.size());

    // javitani - BF merete: BL merete - $$$X edgek szama...
    bit_vector BF(border_sum);


    for (size_t i = 0; i < kmers.size(); ++i) {
        if (i == kmers.size() - 1 || !compare_nodes_begin(kmers[i], kmers[i + 1])) {
            BL[i] = 1;
        }
        else {
            BL[i] = 0;
        }

        // cout << kmer_to_str(kmers[i], k) << " " << BL[i] << endl;

        // get the last character
        auto ac = (uint8_t) (kmers[i] >> (LOGSIGMA * (k - 1)) & mask_last_char).to_ulong();
        // write the last char to the W vector (file)
        *w_bufwriter << ac;
        // edges_static[i] = ac;
        uint8_t aid = bits_to_id(ac) - 1;
        if (aid < SIGMA) {
            BF[char_indexes[aid]] = 1;
            if (prev_indexes[aid] != -1ULL && compare_nodes_end(kmers[i], kmers[prev_indexes[aid]])) {
                BF[char_indexes[aid] - 1] = 0;
            }
            prev_indexes[aid] = i;
            char_indexes[aid]++;
        }

        if ((i + 1) % 1000000 == 0) {
            cerr << (i + 1) << "/" << kmers.size() << " ";
        }
    }

    cerr << endl;


    delete w_bufwriter;
    delete W;

    // correct the starting positions of vector T
    size_t pv = kmers.size() - accumulate(T.begin(), T.end(), 0ULL);
    for (size_t i = 0; i < T.size(); ++i) {
        size_t a = T[i];
        T[i] = pv;
        pv += a;
    }

    sdbg = new SuccinctDeBruijnGraph(kmers.size(), k, colorgram_type, T, BL, BF, tmp_fname);
    tmp_edge_file.close_remove();
    kmers.clear();
    // edges_static.clear();

    // reset the label vector
    // label_hash_vector.resize(sdbg->get_label_vect_size());
    // fill(label_hash_vector.begin(), label_hash_vector.end(), 0);
    X.resize(sdbg->get_label_vect_size(), pair<size_t, size_t>(0, 0));
}


template<uint16_t KMERBITS>
void ColoredDeBrujinGraph<KMERBITS>::build_label_vector(size_t color, const string& dna_str,
                                                        sparse_hash_map<uint64_t, uint64_t>& H,
                                                        sparse_hash_map<uint64_t, uint8_t>& visited) {
    size_t index = 0;
    for (char c : dna_str) {
        uint8_t ac = symbol_to_bits(c);
        if (ac != 255) {
            index = sdbg->get_next_symbol_index(index, ac);
            size_t r = sdbg->get_label_index(index);
            if (r < X.size() && visited.find(r) == visited.end()) {
                visited[r] = 1;
                // add_color(label_hash_vector[r], color);
                if (H.find(X[r].first) == H.end()) {
                    H[X[r].first] = ++max_label;
                    X[r] = pair<size_t, size_t>(max_label, X[r].second + 1);
                }
                else {
                    X[r] = pair<size_t, size_t>(H[X[r].first], X[r].second + 1);
                }
            }
            index = sdbg->forward(index, ac);
        }
    }
    index = sdbg->get_next_symbol_index(index, 0);
    // process the $-tagged edge... (no need to search for the index of the edge - it must be the = with index)
    size_t r = sdbg->get_label_index(index);
    if (r < X.size() && visited.find(r) == visited.end()) {
        // add_color(label_hash_vector[r], color);
        // X[r] = (H.find(X[r].first) == H.end()) ? pair<size_t, size_t>(++max_label, X[r].second + 1)
        //                                        : pair<size_t, size_t>(H[X[r].first], X[r].second + 1);
        visited[r] = 1;
        if (H.find(X[r].first) == H.end()) {
            H[X[r].first] = ++max_label;
            X[r] = pair<size_t, size_t>(max_label, X[r].second + 1);
        }
        else {
            X[r] = pair<size_t, size_t>(H[X[r].first], X[r].second + 1);
        }
    }

    // update the number of colors
    if (color >= C) {
        C = color + 1;
        sdbg->set_C(C);
    }
}


/// Sorts label vector by the frequencies and redistributes the labels accordingly, initializes the color table
/// running time and space complexity: O(|E^+|)
template<uint16_t KMERBITS>
void ColoredDeBrujinGraph<KMERBITS>::sort_label_vector() {
    cerr << "Sorting label vector..." << endl;

    sparse_hash_map<size_t, size_t> freq;
    for (auto x : X) {
        freq[x.first]++;

        if (x.second == 0) {
            cerr << "Error the number of set bits is... " << x.first << " " << x.second << endl;
        }
    }
    vector<vector<size_t>> buckets(X.size() + 1);
    for (auto it = freq.begin(); it != freq.end(); ++it) {
        buckets[it->second].push_back(it->first);
    }

    size_t label = 0;
    for (size_t i = X.size() + 1; i > 0; --i) {
        for (auto j : buckets[i - 1]) {
            freq[j] = label++;
        }
    }

#define SAVE_MULTIPLICITIES
#ifdef SAVE_MULTIPLICITIES
    save_multiplicities(buckets);
#endif
    buckets.clear();

    // initialize color table
    color_table.resize(label, nullptr);
    color_table_succinct.resize(label, nullptr);
    color_table_row_capacities.resize(label, 0);
    label_vector_length = 0;
    for (auto it = X.begin(); it != X.end(); ++it) {
        size_t set_bits = it->second; // miert 0??? nem lehetne...
        size_t l = freq[it->first];
        it->first = l;
        label_vector_length += l + 1;
        if (color_table[l] == nullptr) {
            color_table[l] = new sd_vector_builder(C, set_bits);
            color_table_row_capacities[l] = set_bits;
            color_table_set_bits += set_bits;
        }
    }
}


/// Saves the For debugging purpose
/// \param buckets
template<uint16_t KMERBITS>
void ColoredDeBrujinGraph<KMERBITS>::save_multiplicities(const vector<vector<size_t>>& buckets) {
    cerr << "Saving multiplicities..." << endl;
    string fname = out_fname + ".m";
    sdsl::osfstream out(fname, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        if (sdsl::util::verbose) {
            std::cerr << "ERROR: store_to_file not successful for: '" << fname << "'" << std::endl;
        }
    }
    structure_tree_node *child = structure_tree::add_child(NULL, "", util::class_name(*this));
    size_t written_bytes = 0;
    for (size_t i = X.size(); i > 0; --i) {
        size_t m = i - 1;
        for (size_t j = 0; j < buckets[m].size(); ++j) {
            written_bytes += write_member(m, out, child);
        }
    }
    structure_tree::add_size(child, written_bytes);
    out.close();
}


template<uint16_t KMERBITS>
void ColoredDeBrujinGraph<KMERBITS>::build_color_table(size_t color, const string& dna_str) {
    size_t index = 0;
    for (char c : dna_str) {
        uint8_t ac = symbol_to_bits(c);
        if (ac != 255) {
            index = sdbg->get_next_symbol_index(index, ac);
            size_t r = sdbg->get_label_index(index);
            if (r < X.size()) {
                add_color(color, r);
            }
            index = sdbg->forward(index, ac);
        }
    }
    // process the $-tagged edge... (no need to search for the index of the edge - it must be the = with index)
    size_t r = sdbg->get_label_index(index);
    if (r < X.size()) {
        add_color(color, r);
    }
}


/// Set the color-th position of the X[r]th row of the color table
/// \tparam KMERBITS
/// \param color
/// \param r
template<uint16_t KMERBITS>
void ColoredDeBrujinGraph<KMERBITS>::add_color(size_t color, size_t r) {
    size_t l = X[r].first;
    // if the row is not yet finished
    if (color_table_succinct[l] == nullptr) {
        auto& row = color_table[l];
        // if we haven't added the color-th bit to the actual row
        if (row->tail() <= color) {
            row->set(color);
            // if we have finished the actual row => make it as an sd vector
            if (row->items() >= X[r].second) {
                color_table_succinct[l] = new sd_vector<>(*row);
                delete row;
                row = nullptr;
            }
        }
    }
}


/// Creates succinct label vector, succinct color table and passes to the succinct colored de Bruijn graph
/// \tparam KMERBITS
template<uint16_t KMERBITS>
void ColoredDeBrujinGraph<KMERBITS>::create_succinct_structures() {
    cerr << "Creating succinct label vector..." << endl;
    auto label_vector_builder = new sd_vector_builder(label_vector_length, X.size());
    size_t index = 0;
    for (auto x : X) {
        index += x.first;
        label_vector_builder->set(index);
        ++index;
    }
    sdbg->set_X(label_vector_builder);
    delete label_vector_builder;
    X.clear();

    cerr << "Creating succinct color table..." << endl;
    auto ct_vector_builder = new sd_vector_builder(color_table_succinct.size() * C, color_table_set_bits);
    cerr << "ct size: " << color_table_succinct.size() * C << " " << color_table_set_bits << endl;
    index = 0;
    for (size_t i = 0; i < color_table_succinct.size(); ++i) {
        if (color_table_succinct[i] == nullptr) {
            cerr << "Error, the following row is NULL: " << i << endl;
        }
    }
    cerr << endl;
    for (size_t i = 0; i < color_table_succinct.size(); ++i, index += C) {
        auto X_select = sd_vector<>::select_1_type(color_table_succinct[i]);
        if (color_table_succinct[i] == nullptr) {
            cerr << "Error..." << i << " capacity: " << color_table_row_capacities[i] << endl;
        }
        for (uint32_t j = 1; j <= color_table_row_capacities[i]; ++j) {
            size_t set_index = index + X_select.select(j);
            ct_vector_builder->set(set_index);
        }
        delete color_table_succinct[i];
    }
    color_table.clear();
    color_table_succinct.clear();
    color_table_row_capacities.clear();
    sdbg->set_CT(ct_vector_builder);
    delete ct_vector_builder;
}


template<uint16_t KMERBITS>
void ColoredDeBrujinGraph<KMERBITS>::print_stats() {
    cerr << "Number of edges:           " << sdbg->get_num_of_edges() << endl;
    cerr << "Number of nodes:           " << sdbg->get_num_of_nodes() << endl;
    cerr << "k-mer size:                " << k << endl;
    cerr << "Number of colors:          " << C << endl;
    cerr << "#Rows of Color Table (CT): " << color_table.size() << endl;
    cerr << "Size of CT in bits:        " << color_table.size() * C << endl;
    cerr << "# of set bits in CT:       " << color_table_set_bits << endl;
    // cerr << "# of gaps in CT:           " << gaps.size() << endl;
    cerr << "Size of Label Vector (X):  " << X.size() << endl;
    cerr << endl;
}

// saves the graph
template<uint16_t KMERBITS>
bool ColoredDeBrujinGraph<KMERBITS>::save_graph() {
    return sdbg->save(out_fname);
}

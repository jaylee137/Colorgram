#ifndef COLORGRAM_BUILDER_HPP
#define COLORGRAM_BUILDER_HPP

#include <iostream>
#include "dbg_wrapper.hpp"

using namespace std;


void create_color_information(DBGWrapper *& dbg, const string& kmer_list_fname) {
    // check whether the kmer list files do exist
    ifstream kmer_lst_file(kmer_list_fname);
    string fname;
    vector<string> fnames;
    while (getline(kmer_lst_file, fname)) {
        if (file_exists(fname)) {
            fnames.push_back(fname);
        }
        else {
            cerr << "Unable to open file '" << fname << "'" << endl;
            exit(1);
        }
    }
    kmer_lst_file.close();

    auto build = [&fnames](auto fn) {
        // add colors for the dbg
        for (size_t i = 0; i < fnames.size(); ++i) {
            ifstream f(fnames[i]);
            if (fnames.size() < 1000000 || i % 1000000 == 0) {
                cerr << i << " ";
            }
            sparse_hash_map<uint64_t, uint64_t> H;
            sparse_hash_map<uint64_t, uint8_t> visited;
            if (f.is_open()) {
                string line;
                // the first line is just a header...
                while (getline(f, line)) {
                    // the second line contains the actual (flatten) DNA data
                    getline(f, line);

                    fn(i, line, H, visited);
                }

                f.close();
            }
            else {
                cerr << "Unable to open file '" << fnames[i] << "'" << endl;
                exit(1);
            }
        }
    };

    // first build the label vector
    build([&dbg](size_t i, const string& line, sparse_hash_map<uint64_t, uint64_t>& H,
                 sparse_hash_map<uint64_t, uint8_t>& visited) { dbg->build_label_vector(i, line, H, visited); });

    // then sort it, create the proper labels
    dbg->sort_label_vector();

    // then create color table
    build([&dbg](size_t i, const string& line, sparse_hash_map<uint64_t, uint64_t>& H,
                 sparse_hash_map<uint64_t, uint8_t>& visited) { dbg->build_color_table(i, line); });

    dbg->print_stats();

    // create the succinct label vector and color table
    dbg->create_succinct_structures();
}

void create_color_information_single_file_mode(DBGWrapper *& dbg, const string& fname) {
    auto build = [&fname](auto fn) {
        ifstream f(fname);

        if (f.is_open()) {
            string line;
            size_t i = 0;
            // the first line is just a header...
            while (getline(f, line)) {
                if (i < 1000000 || i % 1000000 == 0) {
                    cerr << i << " ";
                }
                // the second line contains the actual (flatten) DNA data
                getline(f, line);

                sparse_hash_map<uint64_t, uint64_t> H;
                sparse_hash_map<uint64_t, uint8_t> visited;

                fn(i++, line, H, visited); //dbg->build_label_vector(i++, line);
            }

            f.close();
        }
        else {
            cerr << "Unable to open file '" << fname << "'" << endl;
            exit(1);
        }
    };

    // first build the label vector
    build([&dbg](size_t i, const string& line, sparse_hash_map<uint64_t, uint64_t>& H,
                 sparse_hash_map<uint64_t, uint8_t>& visited) { dbg->build_label_vector(i, line, H, visited); });

    // then sort it, create the proper labels
    dbg->sort_label_vector();

    // then create color table
    build([&dbg](size_t i, const string& line, sparse_hash_map<uint64_t, uint64_t>& H,
                 sparse_hash_map<uint64_t, uint8_t>& visited) { dbg->build_color_table(i, line); });

    dbg->print_stats();

    // create the succinct label vector and color table
    dbg->create_succinct_structures();
}

DBGWrapper *build_graph(uint32_t k, uint8_t cologram_type, const string& kmer_list_fname, const string& kmer_db_fname,
                        const string& begin_db_fname, const string& end_db_fname, const string& out_db_fname,
                        bool single_mode) {
    // check that the kmer list file exists
    if (!file_exists(kmer_list_fname)) {
        cerr << "Unable to open file '" << kmer_list_fname << "'" << endl;
        exit(1);
    }

    // create the de Bruijn graph
    auto dbg = new DBGWrapper(k, cologram_type, kmer_db_fname, begin_db_fname, end_db_fname, out_db_fname);

    cerr << "Creating color information..." << endl;
    if (single_mode) {
        create_color_information_single_file_mode(dbg, kmer_list_fname);
    }
    else {
        create_color_information(dbg, kmer_list_fname);
    }
    cerr << endl;

    // save the graph: out_db_fname.dbg, out_db_fname.ct, out_db_fname.x
    dbg->save_graph();

    return dbg;
}

#endif //COLORGRAM_BUILDER_HPP

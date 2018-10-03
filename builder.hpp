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

    // add colors for the dbg
    for (uint32_t i = 0; i < (uint32_t) fnames.size(); ++i) {
        ifstream f(fnames[i]);
        cerr << i << " ";
        if (f.is_open()) {
            string line;
            // the first line is just a header...
            while (getline(f, line)) {
                // the second line contains the actual (flatten) DNA data
                getline(f, line);

                dbg->build_colored_graph(i, line);
            }

            f.close();
        }
        else {
            cerr << "Unable to open file '" << fnames[i] << "'" << endl;
            exit(1);
        }
    }
}

void create_color_information_single_file_mode(DBGWrapper *& dbg, const string& fname) {
    ifstream f(fname);

    if (f.is_open()) {
        string line;
        size_t i = 0;
        // the first line is just a header...
        while (getline(f, line)) {
            cerr << i << " ";
            // the second line contains the actual (flatten) DNA data
            getline(f, line);

            dbg->build_colored_graph(i++, line);
        }

        f.close();
    }
    else {
        cerr << "Unable to open file '" << fname << "'" << endl;
        exit(1);
    }
}

DBGWrapper *build_graph(uint32_t k, uint8_t cologram_type, const string& kmer_list_fname, const string& kmer_db_fname,
                        const string& begin_db_fname, const string& end_db_fname, const string& out_db_fname,
                        bool single_mode, bool save_db = true) {
    // check that the kmer list file exists
    if (!file_exists(kmer_list_fname)) {
        cerr << "Unable to open file '" << kmer_list_fname << "'" << endl;
        exit(1);
    }

    // create the de Bruijn graph
    auto dbg = new DBGWrapper(k, cologram_type, kmer_db_fname, begin_db_fname, end_db_fname);

    cerr << "Creating color information..." << endl;
    if (single_mode) {
        create_color_information_single_file_mode(dbg, kmer_list_fname);
    }
    else {
        create_color_information(dbg, kmer_list_fname);
    }
    cerr << endl;

    dbg->print_stats();

    dbg->sort_color_table();

    if (save_db) {
        // save the graph: out_db_fname.dbg, out_db_fname.ct, out_db_fname.x
        dbg->save_graph(out_db_fname);
    }

    return dbg;
}

#endif //COLORGRAM_BUILDER_HPP

#include <iostream>
#include "dbg_wrapper.hpp"

using namespace std;


int main(int argc, char *argv[]) {
    if (argc < 7) {
        cerr
                << "Usage: cologram <kmer> <cologram-type \\in {0, 1, 2}> <kmer list file> <kmer database>"
                   "<begin kmer database> <end kmer database> <output-file>" << endl;
        cerr << "Example usage: ./cologram 32 0 kmer_lst.txt kmers.kmc begin_kmers.kmc end_kmers.kmc out_db_name"
             << endl;
        exit(1);
    }

    uint32_t k = (uint32_t) stoul(argv[1]);
    uint8_t cologram_type = (uint8_t) stoi(argv[2]);
    string kmer_list_fname = argv[3];
    string kmer_db_fname = argv[4];
    string begin_db_fname = argv[5];
    string end_db_fname = argv[6];
    string out_db_fname = argv[7];

    // check that the kmer list file exists
    if (!file_exists(kmer_list_fname)) {
        cerr << "Unable to open file '" << kmer_list_fname << "'" << endl;
        exit(1);
    }

    // chech whether the kmer list files do exist
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

    // create the de Bruijn graph
    DBGWrapper dbg(k, cologram_type, kmer_db_fname, begin_db_fname, end_db_fname);

    cerr << "Creating color information..." << endl;
    // add colors for the dbg
    size_t color_cnt = fnames.size();
    for (size_t i = 0; i < color_cnt; ++i) {
        ifstream f(fnames[i]);
        cerr << i << " ";
        if (f.is_open()) {
            string line;
            // the first line is just a header...
            while(getline(f, line)) {
                // the second line contains the actual (flatten) DNA data
                getline(f, line);

                dbg.build_colored_graph(i, line);
            }

            f.close();
        }
        else {
            cerr << "Unable to open file '" << fnames[i] << "'" << endl;
            exit(1);
        }
    }
    cerr << endl;

    dbg.print_stats();

    // save the graph: out_db_fname.dbg, out_db_fname.ct, out_db_fname.x
    dbg.save_graph(out_db_fname);

    return 0;
}

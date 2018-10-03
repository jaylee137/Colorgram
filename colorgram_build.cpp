#include <iostream>
#include "builder.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 7) {
        cerr << "Usage: cologram-build <kmer> <cologram-type \\in {0, 1, 2}> <kmer list file> <kmer database>"
                "<begin kmer database> <end kmer database> <output-file>" << endl;
        cerr << "Example usage: ./cologram-build 32 0 kmer_lst.txt kmers.kmc begin_kmers.kmc end_kmers.kmc out_db_name"
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

    DBGWrapper *dbg = build_graph(k, cologram_type, kmer_list_fname, kmer_db_fname, begin_db_fname, end_db_fname,
                                  out_db_fname, argc >= 9);

    delete dbg;

    return 0;
}

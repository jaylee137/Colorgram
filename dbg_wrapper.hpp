#ifndef COLORGRAM_DBG_WRAPPER_HPP
#define COLORGRAM_DBG_WRAPPER_HPP

#include "dbg.cpp"

#define KMER8BYTES 64 / LOGSIGMA
#define KMER16BYTES 128 / LOGSIGMA
#define KMER24BYTES 192 / LOGSIGMA
#define KMER32BYTES 256 / LOGSIGMA
#define KMER40BYTES 320 / LOGSIGMA


class DBGWrapper {
public:
    DBGWrapper(uint32_t pk, uint8_t pcologram_type, const string& pkmer_db_fname, const string& pbegin_db_fname,
               const string& pend_db_fname, const string& pout_fname) : k(pk) {
        if (k <= KMER8BYTES) {
            dbg8 = new ColoredDeBrujinGraph<64>(pk, pcologram_type, pkmer_db_fname, pbegin_db_fname, pend_db_fname,
                                                pout_fname);
        }
        else if (k <= KMER16BYTES) {
            dbg16 = new ColoredDeBrujinGraph<128>(pk, pcologram_type, pkmer_db_fname, pbegin_db_fname, pend_db_fname,
                                                  pout_fname);
        }
        else if (k <= KMER24BYTES) {
            dbg24 = new ColoredDeBrujinGraph<192>(pk, pcologram_type, pkmer_db_fname, pbegin_db_fname, pend_db_fname,
                                                  pout_fname);
        }
        else if (k <= KMER32BYTES) {
            dbg32 = new ColoredDeBrujinGraph<256>(pk, pcologram_type, pkmer_db_fname, pbegin_db_fname, pend_db_fname,
                                                  pout_fname);
        }
        else if (k <= KMER40BYTES) {
            dbg40 = new ColoredDeBrujinGraph<320>(pk, pcologram_type, pkmer_db_fname, pbegin_db_fname, pend_db_fname,
                                                  pout_fname);
        }
        else {
            cerr << "Maximal k-mer size is " << KMER40BYTES << "..." << endl;
        }
    }

    ~DBGWrapper() {
        delete dbg8;
        delete dbg16;
        delete dbg24;
        delete dbg32;
        delete dbg40;
    }

    void build_label_vector(size_t color, const string& dna_str, sparse_hash_map<uint64_t, uint64_t>& H,
                            sparse_hash_map<uint64_t, uint8_t>& visited) {
        if (k <= KMER8BYTES) {
            dbg8->build_label_vector(color, dna_str, H, visited);
        }
        else if (k <= KMER16BYTES) {
            dbg16->build_label_vector(color, dna_str, H, visited);
        }
        else if (k <= KMER24BYTES) {
            dbg24->build_label_vector(color, dna_str, H, visited);
        }
        else if (k <= KMER32BYTES) {
            dbg32->build_label_vector(color, dna_str, H, visited);
        }
        else if (k <= KMER40BYTES) {
            dbg40->build_label_vector(color, dna_str, H, visited);
        }
    }

    void sort_label_vector() {
        if (k <= KMER8BYTES) {
            dbg8->sort_label_vector();
        }
        else if (k <= KMER16BYTES) {
            dbg16->sort_label_vector();
        }
        else if (k <= KMER24BYTES) {
            dbg24->sort_label_vector();
        }
        else if (k <= KMER32BYTES) {
            dbg32->sort_label_vector();
        }
        else if (k <= KMER40BYTES) {
            dbg40->sort_label_vector();
        }
    }

    void build_color_table(size_t color, const string& dna_str) {
        if (k <= KMER8BYTES) {
            dbg8->build_color_table(color, dna_str);
        }
        else if (k <= KMER16BYTES) {
            dbg16->build_color_table(color, dna_str);
        }
        else if (k <= KMER24BYTES) {
            dbg24->build_color_table(color, dna_str);
        }
        else if (k <= KMER32BYTES) {
            dbg32->build_color_table(color, dna_str);
        }
        else if (k <= KMER40BYTES) {
            dbg40->build_color_table(color, dna_str);
        }
    }

    SuccinctDeBruijnGraph *get_sdbg() {
        if (k <= KMER8BYTES) {
            return dbg8->get_sdbg();
        }
        else if (k <= KMER16BYTES) {
            return dbg16->get_sdbg();
        }
        else if (k <= KMER24BYTES) {
            return dbg24->get_sdbg();
        }
        else if (k <= KMER32BYTES) {
            return dbg32->get_sdbg();
        }
        else if (k <= KMER40BYTES) {
            return dbg40->get_sdbg();
        }

        return nullptr;
    }

    void print_stats() {
        if (k <= KMER8BYTES) {
            dbg8->print_stats();
        }
        else if (k <= KMER16BYTES) {
            dbg16->print_stats();
        }
        else if (k <= KMER24BYTES) {
            dbg24->print_stats();
        }
        else if (k <= KMER32BYTES) {
            dbg32->print_stats();
        }
        else if (k <= KMER40BYTES) {
            dbg40->print_stats();
        }
    }

    void create_succinct_structures() {
        if (k <= KMER8BYTES) {
            dbg8->create_succinct_structures();
        }
        else if (k <= KMER16BYTES) {
            dbg16->create_succinct_structures();
        }
        else if (k <= KMER24BYTES) {
            dbg24->create_succinct_structures();
        }
        else if (k <= KMER32BYTES) {
            dbg32->create_succinct_structures();
        }
        else if (k <= KMER40BYTES) {
            dbg40->create_succinct_structures();
        }
    }

    bool save_graph() {
        if (k <= KMER8BYTES) {
            return dbg8->save_graph();
        }
        else if (k <= KMER16BYTES) {
            return dbg16->save_graph();
        }
        else if (k <= KMER24BYTES) {
            return dbg24->save_graph();
        }
        else if (k <= KMER32BYTES) {
            return dbg32->save_graph();
        }
        else if (k <= KMER40BYTES) {
            return dbg40->save_graph();
        }

        return false;
    }


private:
    ColoredDeBrujinGraph<64> *dbg8 = nullptr;
    ColoredDeBrujinGraph<128> *dbg16 = nullptr;
    ColoredDeBrujinGraph<192> *dbg24 = nullptr;
    ColoredDeBrujinGraph<256> *dbg32 = nullptr;
    ColoredDeBrujinGraph<320> *dbg40 = nullptr;

    uint32_t k;
};

#endif // COLORGRAM_DBG_WRAPPER_HPP

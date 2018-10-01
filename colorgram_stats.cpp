#include <iostream>
#include "succinct_dbg.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "Usage: cologram-stats <fname>" << endl;
        cerr << "Example usage: ./cologram-stats out" << endl;
        cerr << "Note that there should be out.dbg, out.x, out.ct, out.sbv..." << endl;
        exit(1);
    }

    string fname = argv[1];
    SuccinctDeBruijnGraph sdbg(fname);
    sdbg.print_stats(cerr);
}

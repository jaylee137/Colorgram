#include <gtest/gtest.h>
#include "builder.hpp"

void test_str_forward(SuccinctDeBruijnGraph *sdbg, string str, vector<size_t> fv) {
    size_t index = 0;
    for (size_t i = 0; i < fv.size(); ++i) {
        uint8_t ac = symbol_to_bits(str[i]);
        index = sdbg->get_next_symbol_index(index, ac);
        index = sdbg->forward(index, ac);
        ASSERT_EQ(index, fv[i]);
    }
}

// void test_str_backward(SuccinctDeBruijnGraph *sdbg, string str, vector<size_t> fv) {
//     size_t index = fv[0];
//     for (size_t i = 1; i < fv.size(); ++i) {
//         index = sdbg->backward(index - sdbg->get_start_node_length());
//         uint8_t ac = symbol_to_bits(str[i]);
//         while (sdbg->get_edge(index) != ac) {
//             index++;
//         }
//         // index = sdbg->get_next_symbol_index(index, ac);
//         ASSERT_EQ(index, fv[i]);
//     }
// }

TEST(DBGTest, DBG1) {
    DBGWrapper *dbg = build_graph(4, 0, "../tests/test_lst.txt", "../tests/edges/kmers", "../tests/edges/begin",
                                  "../tests/edges/end", "", false, false);
    auto sdbg = dbg->get_sdbg();

    ASSERT_EQ(25, sdbg->get_num_of_edges());
    ASSERT_EQ(20, sdbg->get_num_of_nodes());
    ASSERT_EQ(2, sdbg->get_start_node_length());
    ASSERT_EQ(3, sdbg->get_C());
    ASSERT_EQ(4, sdbg->get_k());
    ASSERT_EQ(10, sdbg->get_label_vect_size());
    ASSERT_GE(MAXCOLORS, 4);

    // test edge list
    string edge_lst = "CGGTTCTTCCGATA$AAAAACCGTA";
    for (size_t i = 0; i < edge_lst.length(); ++i) {
        ASSERT_EQ(edge_lst[i], bits_to_char(sdbg->get_edge(i)));
    }

    // test in-degree
    vector<uint8_t> in_degrees{0, 0, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1};
    for (size_t i = 0; i < in_degrees.size(); ++i) {
        ASSERT_EQ(in_degrees[i], sdbg->indegree(i));
    }

    // test out-degree
    vector<uint8_t> out_degrees{2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1};
    for (size_t i = 0; i < out_degrees.size(); ++i) {
        ASSERT_EQ(out_degrees[i], sdbg->outdegree(i));
    }

    // test forward
    test_str_forward(sdbg, "GACTTACAGATC", {16, 5, 11, 22, 24, 9, 13, 2, 17, 6, 21, 14});
    test_str_forward(sdbg, "GACTGACATC", {16, 5, 11, 22, 19, 8, 11, 2, 20, 14});
    test_str_forward(sdbg, "CGATCATC", {10, 18, 7, 21, 14, 4, 20, 14});

    // // test backward
    // test_str_backward(sdbg, "CTAGACATTCAG", {14, 21, 6, 17, 2, 13, 9, 24, 22, 11, 5, 16}, {1, 0, 0, 0, 1, 0, 0, 0});

    // test BL
    vector<uint8_t> BL{0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1};
    for (size_t i = 0; i < BL.size(); ++i) {
        ASSERT_EQ((bool) BL[i], sdbg->get_BL(i));
    }

    // test BF
    vector<uint8_t> BF{0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1};
    for (size_t i = 0; i < BF.size(); ++i) {
        ASSERT_EQ((bool) BF[i], sdbg->get_BF(i));
    }

    // test label vector
    vector<size_t> label_vect{3, 2, 1, 0, 0, 2, 4, 3, 0, 1};
    for (size_t i = 0; i < label_vect.size(); ++i) {
        ASSERT_EQ(label_vect[i], sdbg->get_label(i));
    }

    // test SBV
    vector<uint8_t> SBV{1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0};
    for (size_t i = 0; i < SBV.size(); ++i) {
        ASSERT_EQ((bool) SBV[i], sdbg->get_label_index(i) < sdbg->get_label_vect_size());
    }

    // test CT
    vector<bitset<MAXCOLORS>> CT{bitset<MAXCOLORS>("010"), bitset<MAXCOLORS>("001"), bitset<MAXCOLORS>("011"),
                                 bitset<MAXCOLORS>("100"), bitset<MAXCOLORS>("111")};
    vector<size_t> ncv{1, 1, 2, 1, 3};
    for (size_t i = 0; i < CT.size(); ++i) {
        bitset<MAXCOLORS> ac;
        size_t num_of_colors = 0;
        sdbg->update_color_class(i, ac, num_of_colors);
        ASSERT_EQ(ac, CT[i]);
        ASSERT_EQ(ncv[i], num_of_colors);
    }

    // test getting the color classes of the edges
    vector<bitset<MAXCOLORS>> CT_all{bitset<MAXCOLORS>("010"), bitset<MAXCOLORS>("001"), bitset<MAXCOLORS>("011"),
                                     bitset<MAXCOLORS>("100"), bitset<MAXCOLORS>("111"), bitset<MAXCOLORS>("110"),
                                     bitset<MAXCOLORS>("101")};
    vector<size_t> color_class_ids{3, 2, 1, 0, 3, 2, 1, 3, 0, 1, 3, 0, 2, 1, 4, 3, 2, 1, 3, 0, 5, 6, 0, 1, 1};
    for (size_t i = 0; i < color_class_ids.size(); ++i) {
        bitset<MAXCOLORS> acc;
        sdbg->get_color_class(acc, i);

        ASSERT_EQ(CT_all[color_class_ids[i]], acc);
    }

    delete sdbg;
}


TEST(DBGTest, DBG2) {
    DBGWrapper *dbg = build_graph(4, 0, "../tests/test_kmers_all.txt", "../tests/edges2/kmers", "../tests/edges2/begin",
                                  "../tests/edges2/end", "", true, false);
    auto sdbg = dbg->get_sdbg();

    ASSERT_EQ(27, sdbg->get_num_of_edges());
    ASSERT_EQ(22, sdbg->get_num_of_nodes());
    ASSERT_EQ(1, sdbg->get_start_node_length());
    ASSERT_EQ(3, sdbg->get_C());
    ASSERT_EQ(4, sdbg->get_k());
    ASSERT_EQ(10, sdbg->get_label_vect_size());
    ASSERT_GE(MAXCOLORS, 4);

    // test edge list
    string edge_lst = "CGTGC$CGTAAT$CCAAGTACATAACC";
    for (size_t i = 0; i < edge_lst.length(); ++i) {
        ASSERT_EQ(edge_lst[i], bits_to_char(sdbg->get_edge(i)));
    }

    // test in-degree
    vector<uint8_t> in_degrees{0, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1};
    for (size_t i = 0; i < in_degrees.size(); ++i) {
        ASSERT_EQ(in_degrees[i], sdbg->indegree(i));
    }
    //
    // test out-degree
    vector<uint8_t> out_degrees{1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};
    for (size_t i = 0; i < out_degrees.size(); ++i) {
        ASSERT_EQ(out_degrees[i], sdbg->outdegree(i));
    }

    // test BL
    vector<uint8_t> BL{1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1};
    for (size_t i = 0; i < BL.size(); ++i) {
        ASSERT_EQ((bool) BL[i], sdbg->get_BL(i));
    }

    // test BF
    vector<uint8_t> BF{0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1};
    for (size_t i = 0; i < BF.size(); ++i) {
        ASSERT_EQ((bool) BF[i], sdbg->get_BF(i));
    }

    // test CT
    vector<bitset<MAXCOLORS>> CT{bitset<MAXCOLORS>("010"), bitset<MAXCOLORS>("001"), bitset<MAXCOLORS>("100"),
                                 bitset<MAXCOLORS>("011"), bitset<MAXCOLORS>("110")};
    vector<size_t> ncv{1, 1, 1, 2, 2};
    for (size_t i = 0; i < CT.size(); ++i) {
        bitset<MAXCOLORS> ac;
        size_t num_of_colors = 0;
        sdbg->update_color_class(i, ac, num_of_colors);
        ASSERT_EQ(ac, CT[i]);
        ASSERT_EQ(ncv[i], num_of_colors);
    }

    bitset<MAXCOLORS> acc;
    sdbg->get_color_class(acc, 0);
    // get color classes
    ASSERT_EQ(bitset<MAXCOLORS>(string(3, '1')), acc);
}


/// Before running the tests make sure that the kmer databases were created
/// colorgram_tool 4 0 ../tests/test_lst.txt ../tests/edges/kmers ../tests/edges/begin ../tests/edges/end out_tests
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

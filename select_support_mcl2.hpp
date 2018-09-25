#ifndef COLORGRAM_SELECT_SUPPORT_MCL2_H
#define COLORGRAM_SELECT_SUPPORT_MCL2_H

#include <sdsl/select_support_mcl.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sparsepp/spp.h>
#include <stxxl.h>

#include "config.h"

using namespace spp;
using namespace sdsl;

class select_support_mcl2 : public select_support_mcl<1, 1> {
public:
    explicit select_support_mcl2(const int_vector_type& label_hash_vector, const int_vector_type& label_permutation,
                                 sparse_hash_map<uint64_t, uint64_t>& cids) : select_support_mcl<1, 1>(nullptr) {
        set_vector(nullptr);
        initData();

        // calculate the length of the labels
        size_type length = 0;
        for (auto i : label_hash_vector) {
            length += label_permutation[cids[i]] + 1;
        }

        m_logn = bits::hi(length) + 1;
        m_logn2 = m_logn * m_logn;
        m_logn4 = m_logn2 * m_logn2;

        // if (select_support_mcl<t_b,t_pat_len>::m_v==nullptr)
        //     return;
        // Count the number of arguments in the bit vector
        m_arg_cnt = label_hash_vector.size(); // arg_cnt(*v);

        const size_type SUPER_BLOCK_SIZE = 4096;

        if (m_arg_cnt == 0) // if there are no arguments in the vector we are done...
            return;

        size_type sb = (m_arg_cnt + SUPER_BLOCK_SIZE - 1) / SUPER_BLOCK_SIZE; // number of superblocks
        delete[] m_miniblock;
        m_miniblock = new int_vector<0>[sb];

        m_superblock = int_vector<0>(sb, 0, m_logn);


        size_type arg_position[SUPER_BLOCK_SIZE], arg_cnt = 0;
        size_type sb_cnt = 0;

        size_type i = 0;
        for (auto j : label_hash_vector) {
            size_t aid = label_permutation[cids[j]];
            i += aid;


            arg_position[arg_cnt % SUPER_BLOCK_SIZE] = i;
            assert(arg_position[arg_cnt % SUPER_BLOCK_SIZE] == i);
            ++arg_cnt;
            if (arg_cnt % SUPER_BLOCK_SIZE == 0 or arg_cnt == m_arg_cnt) { //
                assert(sb_cnt < sb);
                m_superblock[sb_cnt] = arg_position[0];

                size_type pos_diff = arg_position[(arg_cnt - 1) % SUPER_BLOCK_SIZE] - arg_position[0];
                if (pos_diff > m_logn4) { // longblock
                    if (m_longsuperblock == nullptr)
                        m_longsuperblock = new int_vector<0>[sb]; // create longsuperblock
                    m_longsuperblock[sb_cnt] = int_vector<0>(SUPER_BLOCK_SIZE, 0,
                                                             bits::hi(arg_position[(arg_cnt - 1) % SUPER_BLOCK_SIZE]) + 1);

                    for (size_type j = 0; j <= (arg_cnt - 1) % SUPER_BLOCK_SIZE; ++j)
                        m_longsuperblock[sb_cnt][j] = arg_position[j]; // copy argument positions to longsuperblock
                }
                else { // short block
                    m_miniblock[sb_cnt] = int_vector<0>(64, 0, bits::hi(pos_diff) + 1);
                    for (size_type j = 0; j <= (arg_cnt - 1) % SUPER_BLOCK_SIZE; j += 64) {
                        m_miniblock[sb_cnt][j / 64] = arg_position[j] - arg_position[0];
                    }
                }
                ++sb_cnt;
            }
            i++;
        }
    }

    select_support_mcl2() : select_support_mcl<1, 1>(nullptr) {}
};

#endif //COLORGRAM_SELECT_SUPPORT_MCL2_H

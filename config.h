#ifndef COLORGRAM_CONFIG_H
#define COLORGRAM_CONFIG_H

#include <stxxl.h>
#include <bitset>

#define SIGMA       4u
#define LOGSIGMA    3u

static const char base[SIGMA + 1] = {'$', 'A', 'C', 'G', 'T'};

struct color_class_t {
    color_class_t() = default;

    explicit color_class_t(const std::bitset<MAXCOLORS>& pbitvector) : bitvector(pbitvector) {}

    std::bitset<MAXCOLORS> bitvector;
    size_t cnt = 0;
};

typedef typename stxxl::VECTOR_GENERATOR<size_t>::result size_t_vector_type;
typedef typename stxxl::VECTOR_GENERATOR<uint8_t>::result uint8_t_vector_type;


#endif //COLORGRAM_CONFIG_H

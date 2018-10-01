cmake_minimum_required(VERSION 2.6.2 FATAL_ERROR)
project(Colorgram)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS_DEBUG "-pipe -m64 -std=c++14 -pedantic-errors -W -Wall -Wextra -Wpointer-arith -Wunused -Wwrite-strings -ggdb -O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing -fopenmp")
set(CMAKE_CXX_FLAGS "-pipe -m64 -std=c++14 -pedantic-errors -W -Wall -Wextra -Wpointer-arith -Wunused -Wwrite-strings -O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing -fopenmp")

IF (NOT DEFINED MAXCOLORS)
    SET(MAXCOLORS 569)
ENDIF()

add_definitions(-DMAXCOLORS=${MAXCOLORS} -DNDEBUG)

include_directories(SYSTEM 3rd_party/KMC
                    SYSTEM 3rd_party/sdsl-lite/install/include
                    SYSTEM 3rd_party/stxxl/include
                    SYSTEM 3rd_party/sparsepp
                    SYSTEM 3rd_party/googletest/googletest/include)

cmake_policy(SET CMP0015 NEW)

link_directories(3rd_party/stxxl/build/install/lib
                 3rd_party/sdsl-lite/install/lib
                 3rd_party/googletest/googletest/build/)

add_executable(colorgram-build colorgram_build.cpp succinct_dbg.cpp)
target_link_libraries(colorgram-build
                      sdsl
                      stxxl
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/kmc_file.o
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/kmer_api.o
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/mmer.o)

add_executable(colorgram-stats colorgram_stats.cpp succinct_dbg.cpp)
target_link_libraries(colorgram-stats
                      sdsl
                      stxxl
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/kmc_file.o
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/kmer_api.o
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/mmer.o)

# tests
add_executable(succinct-dbg-tests succinct_dbg_tests.cpp succinct_dbg.cpp)
target_link_libraries(succinct-dbg-tests
                      gtest
                      pthread
                      sdsl
                      stxxl
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/kmc_file.o
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/kmer_api.o
                      ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/KMC/kmc_api/mmer.o)

cmake_minimum_required(VERSION 2.6.2 FATAL_ERROR)
project(Colorgram)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS_DEBUG "-pipe -m64 -std=c++14 -pedantic-errors -W -Wall -Wextra -Wpointer-arith -Wunused -Wwrite-strings -ggdb -O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing -fopenmp")
set(CMAKE_CXX_FLAGS "-pipe -m64 -std=c++14 -pedantic-errors -W -Wall -Wextra -Wpointer-arith -Wunused -Wwrite-strings -O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing -fopenmp")

IF (NOT DEFINED MAXCOLORS)
    SET(MAXCOLORS 100)
ENDIF()

add_definitions(-DMAXCOLORS=${MAXCOLORS})

include_directories(SYSTEM 3rd_party/KMC
                    SYSTEM 3rd_party/sdsl-lite/install/include
                    SYSTEM 3rd_party/stxxl/include
                    SYSTEM 3rd_party/sparsepp)

link_directories(3rd_party/stxxl/build/install/lib)

add_executable(colorgram-build cologram_build.cpp succinct_dbg.cpp)
target_link_libraries(colorgram
                      /home/denes/PhD/singapore/color_graph/Colorgram/3rd_party/sdsl-lite/install/lib/libsdsl.a
                      stxxl
                      /home/denes/PhD/singapore/color_graph/cosmo_old/3rd_party_src/KMC/kmc_api/kmc_file.o
                      /home/denes/PhD/singapore/color_graph/cosmo_old/3rd_party_src/KMC/kmc_api/kmer_api.o
                      /home/denes/PhD/singapore/color_graph/cosmo_old/3rd_party_src/KMC/kmc_api/mmer.o)

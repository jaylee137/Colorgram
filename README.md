# Colorgram
Succinct Colored de Bruijn Graph

## Creating the environment

Colorgram requirese five third party packages:

1. KMC tool - https://github.com/refresh-bio/KMC
2. STXXL - https://github.com/stxxl/stxxl
3. SDSL-Lite - https://github.com/simongog/sdsl-lite
4. Sparsepp - https://github.com/greg7mdp/sparsepp
5. Google Test - https://github.com/google/googletest

You can simply run `./create_environment.sh` that creates the 3rd_party folder, clones and complies everything there automatically.

## Compiling and running Colorgram

To create the succinct colored de Bruijn graph run `colorgram_tool.sh` with the following parameters:

`./colorgram_tool.sh -k=32 <file-list> <kmc-files-dir> <out-filename>`

OR

`./colorgram_tool.sh -k=32 <fasta-file> <kmc-files-dir> <out-filename> -m`

In the first case you have to specify a `<file-list>` that contains file names (with absolute paths) where each file is in FASTA-format.
In this case each line in the `<file-list>` correspond to one color. The `<kmc-files-dir>` should be an arbitrary empty (or non-existent)
directory that is the temporary working directory storing the KMC files output by the KMC tool. Finally the `<out-filename>`
should be the name of your desired output file. Colorgram will create 3 different files: out-filename.dbg (storing the succinct
de Bruijn graph), out-filename.x (storing the label vector) and out-filename.ct (storing the color table).

The second case is when you specify the `-m` argument. Here you have to give one FASTA-format file `<fasta-file>`. The
number of colors will be the the number of reads in `<fasta-file>`.

Note that in both cases the FASTA-files should be valid meaning that each read should contain only `ACTG` characters and
for technical reasons each read should be stored on ONE line. For examples you can check the Testing part.

## Statistics and Bubble Calling
If you want to use the features of Colorgram, you have to manually compile `colorgram-stats` program. By default when you
run the `colorgram_tool.sh` it creates a `build` folder and compiles the programs there. You can simply run
`build/colorgram-stats <out-filename>` where you have to specify the same file name that you used for the building (there
must exist out-filename.dbg, out-filename.x and out-filename.ct files).

Note that Bubble Calling is not yet added.

## K-mer size
Colorgram was tested with k<=63 values - note that KMC tool has a bug that sometimes close to k=64 it generates corrupted
k-mer values.

## Testing
You can simply run

`./colorgram_tool.sh -k=4 tests/test_lst.txt tests/edges tests/out_tests`

AND

`./colorgram_tool.sh -k=4 tests/test_kmers_all.fna tests/edges2 tests/out_tests2 -m`

Then you can cd to `build` and run `./succinct-dbg-tests`.

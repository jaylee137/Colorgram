#!/bin/bash

if test "$#" -lt 4; then
    echo "Too few arguments were given..."
elif test "$#" -eq 4; then
    # example usage: ./colorgra_tool.sh -k=32 file_list.txt kmc-files-dir out-filename

    mkdir -p $3/
    rm $3/*kmc_pre
    rm $3/*kmc_suf
    # cut the beginnings, ends of the lists
    python cut_sides.py "$@"


    #############CALCULATE THE SIDES#############

    # calculate k - 1
    K_1="$((${1/-k=/}-1))"

    # create temp tidrectory for kmc
    mkdir -p kmc_temp
    ls -d -1 $3/* | xargs -l -i ./3rd_party/KMC/bin/kmc -ci0 -fm -k$K_1 -b -m1 -cs1 {} {} kmc_temp

    # put the kmer
    begin_array=($(ls -d -1 $3/*begin))
    end_array=($(ls -d -1 $3/*end))


    if test ${#begin_array[@]} -eq 1; then
        mv ${begin_array[0]}.kmc_pre $3/begin.kmc_pre
        mv ${begin_array[0]}.kmc_suf $3/begin.kmc_suf

        mv ${end_array[0]}.kmc_pre $3/end.kmc_pre
        mv ${end_array[0]}.kmc_suf $3/end.kmc_suf
    else
        # begin edges - union
        #kmc_tools simple out1 out2 intersect out_ins
        mkdir -p $3/begin_tmp
        ./3rd_party/KMC/bin/kmc_tools simple ${begin_array[0]} ${begin_array[1]} union $3/begin_tmp/begin
        for (( i=1; i < ${#begin_array[@]}; i++ )); do
            ./3rd_party/KMC/bin/kmc_tools simple ${begin_array[$i]} $3/begin_tmp/begin union $3/begin
            mv $3/begin.kmc* $3/begin_tmp/
        done
        mv $3/begin_tmp/* $3/

        # end edges - union
        mkdir -p $3/end_tmp
        ./3rd_party/KMC/bin/kmc_tools simple ${end_array[0]} ${end_array[1]} union $3/end_tmp/end
        for (( i=1; i < ${#end_array[@]}; i++ )); do
            ./3rd_party/KMC/bin/kmc_tools simple ${end_array[$i]} $3/end_tmp/end union $3/end
            mv $3/end.kmc* $3/end_tmp/
        done
        mv $3/end_tmp/* $3/
    fi



    # #############CALCULATE THE EDGES#############

    # get value k
    K="$((${1/-k=/}))"
    cat $2 | xargs -l -i ./3rd_party/KMC/bin/kmc -ci0 -fm -k$K -b -m1 -cs1 {} {} kmc_temp

    # move files to the edge folder
    # cat $2 | xargs -l -i bash -c 'ls {}.kmc*' | xargs -l -i mv {} $3/

    reads=( $(cat $2) )
    # echo "${reads[@]}"

    if test ${#begin_array[@]} -eq 1; then
        mv ${reads[0]}.kmc_pre $3/kmers.kmc_pre
        mv ${reads[0]}.kmc_suf $3/kmers.kmc_suf
    else
        # if there are more reads => put the union of the edges into one database
        mkdir -p $3/kmers_tmp
        ./3rd_party/KMC/bin/kmc_tools simple ${reads[0]} ${reads[1]} union $3/kmers_tmp/kmers
        for (( i=1; i < ${#reads[@]}; i++ )); do
            ./3rd_party/KMC/bin/kmc_tools simple ${reads[$i]} $3/kmers_tmp/kmers union $3/kmers
            mv $3/kmers.kmc* $3/kmers_tmp/
        done
        mv $3/kmers_tmp/* $3/
    fi

    # count the number of lines in the given kmer list file
    MAXCOLORS=`cat $2 | sed '/^\s*$/d' | wc -l`
    echo compiling with $MAXCOLORS colors

    # build colorgram
    mkdir -p build
    cd build
    cmake ../ -DMAXCOLORS=$MAXCOLORS
    make
    cd ../

    echo build/colorgram-build $K 0 $2 $3/kmers $3/begin $3/end $4
    time build/colorgram-build $K 0 $2 $3/kmers $3/begin $3/end $4
else
    # example usage: ./colorgra_tool.sh -k=32 file_list.txt kmc-files-dir out-filename -m
    mkdir -p $3/
    rm $3/*kmc_pre
    rm $3/*kmc_suf
    # cut the beginnings, ends of the lists
    python cut_sides.py "$@"


    # #############CALCULATE THE SIDES#############

    # calculate k - 1
    K_1="$((${1/-k=/}-1))"

    # create temp tidrectory for kmc
    mkdir -p kmc_temp
    ls -d -1 $3/* | xargs -l -i ./3rd_party/KMC/bin/kmc -ci0 -fm -k$K_1 -b -m1 -cs1 {} {} kmc_temp
    # renaming the begin and end kmer databases...
    `mv $3/*begin.kmc_pre $3/begin.kmc_pre`
    `mv $3/*begin.kmc_suf $3/begin.kmc_suf`
    `mv $3/*end.kmc_pre $3/end.kmc_pre`
    `mv $3/*end.kmc_suf $3/end.kmc_suf`

    # #############CALCULATE THE EDGES#############

    # get value k
    K="$((${1/-k=/}))"
    ./3rd_party/KMC/bin/kmc -ci0 -fm -k$K -b -m1 -cs1 $2 $2 kmc_temp
    
    `mv $2.kmc_pre $3/kmers.kmc_pre`
    `mv $2.kmc_suf $3/kmers.kmc_suf`

    # count the number of lines in the given kmer list file
    MAXCOLORS=`cat $2 | sed '/^\s*$/d' | wc -l`
    MAXCOLORS=$((MAXCOLORS / 2))
    echo compiling with $MAXCOLORS colors

    # build colorgram
    mkdir -p build
    cd build
    cmake ../ -DMAXCOLORS=$MAXCOLORS
    make
    cd ../

    echo build/colorgram-build $K 0 $2 $3/kmers $3/begin $3/end $4 -m
    time build/colorgram-build $K 0 $2 $3/kmers $3/begin $3/end $4 -m
fi

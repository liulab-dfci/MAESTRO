Get source code

    git clone https://github.com/samtools/htslib.git
    cd htslib
    autoheader
    autoconf
    ./configure --disable-bz2 --disable-lzma --enable-libcurl
    make -j 20
    export HTSLIB_ROOT=`pwd`
    cd ..

    git clone https://github.com/ryanlayer/giggle.git
    cd giggle
    make
    export GIGGLE_ROOT=`pwd`
    cd ..

    wget -O gargs https://github.com/brentp/gargs/releases/download/v0.3.6/gargs_linux
    chmod +x gargs


Get database data and index

    mkdir data
    cd data
    wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz
    mkdir orig
    tar zxvf all.mnemonics.bedFiles.tgz -C orig/
    mkdir split

    pip install toolshed --user

    python $GIGGLE_ROOT/examples/rme/rename.py \
        $GIGGLE_ROOT/examples/rme/states.txt \
        $GIGGLE_ROOT/examples/rme/EDACC_NAME.txt \
        "orig/*gz" \
        "split/"

    cd split
    ls *.bed | ../gargs -p 30 "bgzip {}"
    cd ..

    mkdir split_sort

    $GIGGLE_ROOT/scripts/sort_bed "split/*gz" split_sort/ 30

    time $GIGGLE_ROOT/bin/giggle index \
        -i "split_sort/*gz" \
        -o split_sort_b \
        -f -s
    Indexed 55605005 intervals.

    real    1m19.609s
    user    1m16.027s
    sys     0m3.127s

Run experiment

    wget https://www.nature.com/nature/journal/v518/n7539/extref/nature13835-s1.xls
    # copy first five columns into a text file named gwas.txt

    mkdir gwas_results

    for D in `tail -n+2 gwas.txt | awk '{print $1;}'  | sort -u`; do
        grep $D gwas.txt \
        | awk '{OFS="\t"; print $4,$5,$5+1;}' \
        | bgzip -c > $D.bed.gz

        $GIGGLE_ROOT/bin/giggle search \
            -i rme_data/split_sort_b \
            -q $D.bed.gz -s \
        > gwas_results/$D.bed.gz.result

    done

    $GIGGLE_ROOT/scripts/giggle_heat_map.py \
        -s $GIGGLE_ROOT/examples/rme/states.txt \
        --state_names $GIGGLE_ROOT/examples/rme/short_states.txt \
        -c $GIGGLE_ROOT/examples/rme/EDACC_NAME.txt \
        -i gwas_results/Crohns_disease.bed.gz.result \
        -o Crohns_disease.bed.gz.result.pdf \
        -n $GIGGLE_ROOT/examples/rme/new_groups.txt \
        --x_size 3 \
        --y_size 11 \
        --stat combo \
        --ablines 15,26,31,43,52,60,72,82,87,89,93,101,103,116,120,122,127 \
        --group_names $GIGGLE_ROOT/examples/rme/new_groups_names.txt

    $GIGGLE_ROOT/scripts/giggle_gwas_heatmap.py \
        -i "gwas_results/*" \
        --states Enhancers,Strong_transcription \
        -s $GIGGLE_ROOT/examples/rme/states.txt \
        -c $GIGGLE_ROOT/examples/rme/EDACC_NAME.txt \
        -n $GIGGLE_ROOT/examples/rme/new_groups.txt \
        --ablines 15,26,31,43,52,60,72,82,87,89,93,101,103,116,120,122,127 \
        --group_names $GIGGLE_ROOT/examples/rme/new_groups_names.txt \
        -o gwas.pdf \
        --x_size 6 \
        --y_size 11

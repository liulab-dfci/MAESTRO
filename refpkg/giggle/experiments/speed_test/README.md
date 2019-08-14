Get source code 

    git clone https://github.com/arq5x/bedtools2.git
    cd bedtools2 
    make -j 20
    export BEDTOOLS_ROOT=`pwd`
    cd ..

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
    make -j 20
    export GIGGLE_ROOT=`pwd`
    cd ..

    wget -O gsort https://github.com/brentp/gsort/releases/download/v0.0.4/gsort_linux_amd64
    chmod +x gsort

    wget -O gargs https://github.com/brentp/gargs/releases/download/v0.3.6/gargs_linux
    chmod +x gargs

Get database data and index

    mkdir ucsc_data
    cd ucsc_data

    rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database ./

    mkdir parsed_tracks

    ls database/*sql \
    | gargs -p 30 \
        "python $GIGGLE_ROOT/examples/ucsc/parse_sql.py {} parsed_tracks/"

    mkdir parsed_tracks_sorted

    $GIGGLE_ROOT/scripts/sort_bed "parsed_tracks/*gz" parsed_tracks_sorted 30

    time $GIGGLE_ROOT/bin/giggle \
        index \
        -i "parsed_tracks_sorted/*gz" \
        -o parsed_tracks_sorted_b \
        -s \
        -f
    Indexed 6980993757 intervals.

    real    268m46.033s
    user    245m44.262s
    sys     11m51.567s

    cd ..

    mkdir rme_data
    cd rme_data
    wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz
    mkdir orig
    tar zxvf all.mnemonics.bedFiles.tgz -C orig/
    mkdir split

    pip install --user toolshed

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

    cd split_sort
    time ls *.bed.gz | ../gargs "tabix {}"

    real    2m49.350s
    user    2m13.874s
    sys     0m34.953s

    cd ..

    time $GIGGLE_ROOT/bin/giggle index \
        -i "split_sort/*gz" \
        -o split_sort_b \
        -f -s
    Indexed 55605005 intervals.

    real    1m22.409s
    user    1m18.907s
    sys     0m3.379s

    cd ..

Random query sets

    cd rme_data

    zcat split_sort/Adipose_Nuclei_Repressed_PolyComb.bed.gz \
    | cut -f1 \
    | uniq > rme_chrm_order.txt

    for s in `cat rme_chrm_order.txt`; do
        grep -w $s $GIGGLE_ROOT/data/human.hg19.genome
    done \
    > rme.human.hg19.genome

    Q_SIZES="10 100 1000 10000 100000 1000000"
    for Q_SIZE in $Q_SIZES; do
        bedtools random -n $Q_SIZE -g rme.human.hg19.genome  \
        | ../gsort /dev/stdin rme.human.hg19.genome \
        | bgzip -c \
        > rme_r$Q_SIZE.bed.gz
    done

    cd ..

    cd ucsc_data

    zcat parsed_tracks_sorted/snp147.bed.gz \
    | cut -f1 \
    | uniq \
    > ucsc_chrm_order.txt

    for s in `cat ucsc_chrm_order.txt`; do
        grep -w $s $GIGGLE_ROOT/data/human.hg19.genome
    done \
    > ucsc.human.hg19.genome

    Q_SIZES="10 100 1000 10000 100000 1000000"
    for Q_SIZE in $Q_SIZES; do
        bedtools random -n $Q_SIZE -g ucsc.human.hg19.genome  \
        | ../gsort /dev/stdin ucsc.human.hg19.genome \
        | bgzip -c \
        > ucsc_r$Q_SIZE.bed.gz
    done
    cd ..

Speed tests

    cd rme_data

    export RESULTS=rme.speed_test

    Q_SIZES="10 100 1000 10000 100000 1000000"
    for Q_SIZE in $Q_SIZES; do
        $GIGGLE_ROOT/experiments/speed_test/speed_test.sh \
            rme_r$Q_SIZE.bed.gz \
            split_sort \
            rme.human.hg19.genome
    done \
    > $RESULTS

    (cat $RESULTS  |  awk '$1=="giggle"' |  awk '{print $4;}' |  paste -sd " " -;
     cat $RESULTS  |  awk '$1=="bedtools"' | awk '{print $4;}' |  paste -sd " " -;
     cat $RESULTS  |  awk '$1=="tabix"' | awk '{print $4;}' |  paste -sd " " -) \
    | $GIGGLE_ROOT/scripts/lines.py \
        -o $RESULTS.pdf \
        --legend_loc 4 \
        --plot_width 4 \
        --plot_height 3 \
        --xticks ",10,100,1000,1e4,1e5,1e6" \
        --xlabel "Number query intervals" \
        --ylabel "Run time (s)" \
        --ylog \
        -c "#00405B,#4A7C96,#0091D6" \
        --x_min -0.1  --x_max 5.1 \
        --y_min 2e-3  --y_max 1000000 \
        --legend "GIGGLE,BEDTOOLS,TABIX" \
        --legend_loc 4

    cd ..

    cd ucsc_data

    export RESULTS=ucsc.speed_test

    Q_SIZES="10 100 1000 10000 100000 1000000"
    for Q_SIZE in $Q_SIZES; do
        $GIGGLE_ROOT/experiments/speed_test/speed_test.sh \
            ucsc_r$Q_SIZE.bed.gz \
            parsed_tracks_sorted \
            ucsc.human.hg19.genome
    done \
    > $RESULTS
    
    (cat $RESULTS  |  awk '$1=="giggle"' |  awk '{print $4;}' |  paste -sd " " -;
     cat $RESULTS  |  awk '$1=="bedtools"' | awk '{print $4;}' |  paste -sd " " -;
     cat $RESULTS  |  awk '$1=="tabix"' | awk '{print $4;}' |  paste -sd " " -) \
    | $GIGGLE_ROOT/scripts/lines.py \
        -o $RESULTS.pdf \
        --legend_loc 4 \
        --plot_width 4 \
        --plot_height 3 \
        --xticks ",10,100,1000,1e4,1e5,1e6" \
        --xlabel "Number query intervals" \
        --ylabel "Run time (s)" \
        --ylog \
        -c "#00405B,#4A7C96,#0091D6" \
        --x_min -0.1  --x_max 5.1 \
        --y_min 2e-3  --y_max 1000000 \
        --legend "GIGGLE,BEDTOOLS,TABIX" \
        --legend_loc 4

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

    mkdir rme_data
    cd rme_data
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

    mkdir split_sort_raw
    cd split_sort_raw
    cp ../split_sort/* .
    ls *.bed.gz | gargs -p 30 "bgzip -d {}"
    cd ..
    
    time $GIGGLE_ROOT/bin/giggle index \
        -i "split_sort/*gz" \
        -o split_sort_b \
        -f -s
    Indexed 55605005 intervals.

    real    1m20.940s
    user    1m17.933s
    sys     0m2.969s
    
    cd ..

Run experiment

    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1218nnn/GSM1218850/suppl/GSM1218850_MB135DMMD.peak.txt.gz

    zcat GSM1218850_MB135DMMD.peak.txt.gz \
    | awk '$8 >= 100'
    | $HTSLIB_ROOT/bgzip -c > GSM1218850_MB135DMMD.peak.q100.bed.gz

    time $GIGGLE_ROOT/bin/giggle search \
    -i rme_data/split_sort_b \
    -q GSM1218850_MB135DMMD.peak.q100.bed.gz \
    -s \
    > GSM1218850_MB135DMMD.peak.q100.bed.gz.giggle.result

    real    0m0.968s
    user    0m0.218s
    sys     0m0.123s

    $GIGGLE_ROOT/scripts/giggle_heat_map.py \
        -s $GIGGLE_ROOT/examples/rme/states.txt \
        -c $GIGGLE_ROOT/examples/rme/EDACC_NAME.txt \
        -i GSM1218850_MB135DMMD.peak.q100.bed.gz.result \
        -o GSM1218850_MB135DMMD.peak.q100.bed.gz.result.3x11.pdf \
        -n $GIGGLE_ROOT/examples/rme/new_groups.txt \
        --x_size 3 \
        --y_size 11 \
        --stat combo \
        --ablines 15,26,31,43,52,60,72,82,87,89,93,101,103,116,120,122,127 \
        --state_names $GIGGLE_ROOT/examples/rme/short_states.txt \
        --group_names $GIGGLE_ROOT/examples/rme/new_groups_names.txt
        #--no_ylabels \

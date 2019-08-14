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

    git clone https://github.com/arq5x/bits.git
    cd bits
    make
    export BITS_ROOT=`pwd`
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
    ls *.bed.gz | ../gargs -p 30 "bgzip -d {}"
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

    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1218nnn/GSM1218850/suppl/GSM1218850_MB135DMMD.peak.txt.gz

    zcat GSM1218850_MB135DMMD.peak.txt.gz \
    | awk '$8>=100' \
    | $HTSLIB_ROOT/bgzip -c > GSM1218850_MB135DMMD.peak.q100.bed.gz

    time $GIGGLE_ROOT/bin/giggle search \
    -i rme_data/split_sort_b \
    -q GSM1218850_MB135DMMD.peak.q100.bed.gz \
    -s \
    > GSM1218850_MB135DMMD.peak.q100.bed.gz.giggle.result

    real    0m0.945s
    user    0m0.151s
    sys     0m0.163s

    zcat GSM1218850_MB135DMMD.peak.q100.bed.gz \
    > GSM1218850_MB135DMMD.peak.q100.bed

    time \
    ls rme_data/split_sort_raw \
    | ./gargs -p 30 \
    "echo -ne \"{}\\t\"; $BITS_ROOT/bin/bits_test -g $BITS_ROOT/genomes/human.hg19.genome -a GSM1218850_MB135DMMD.peak.q100.bed -b rme_data/split_sort_raw/{} -n 1000" \
    > GSM1218850_MB135DMMD.peak.q100.bed.bits.result

    real    5m14.458s
    user    153m57.696s
    sys     0m15.429s

    paste \
        <(tail -n+2 GSM1218850_MB135DMMD.peak.q100.bed.gz.giggle.result | sort) \
        <(cat GSM1218850_MB135DMMD.peak.q100.bed.bits.result | sort) \
    | cut -f7,14 \
    | sed -e "s/p://g" \ 
    | $GIGGLE_ROOT/scripts/scatter.py \
        -a 0.25 \
        -o mc_obs_fisher_pval.pdf \
        --y_label "MC p-value" \
        --x_label "Fisher's exact p-value (GIGGLE)" \
        --fig_x 3 \
        --fig_y 3 \
        --x_max 1.01 --x_min 0 \
        --y_max 1.01 --y_min 0

    paste \
        <(tail -n+2 GSM1218850_MB135DMMD.peak.q100.bed.gz.giggle.result | sort) \
        <(cat GSM1218850_MB135DMMD.peak.q100.bed.bits.result | sort) \
    | cut -f7,14 \
    | sed -e "s/p://g" \ 
    | awk '{OFS="\t"; print -1*(log($1)/log(10)),-1*(log($2)/log(10));}' \
    | $GIGGLE_ROOT/scripts/scatter.py \
        -a 0.25 \
        -o mc_obs_fisher_pval.1og10.pdf \
        --y_label "-log10(MC p-value)" \
        --x_label "-log10(Fisher's exact p-value) (GIGGLE)" \
        --fig_x 3 \
        --fig_y 3 \
        --x_max 1.01 --x_min 0 \
        --y_max 1.01 --y_min 0

    paste \
        <(tail -n+2 GSM1218850_MB135DMMD.peak.q100.bed.gz.giggle.result | sort) \
        <(cat GSM1218850_MB135DMMD.peak.q100.bed.bits.result | sort) \
    | cut -f4,11,12 \
    | sed -e "s/.://g" \ 
    | awk '$1!=0 && $2!=0 && $3!=0' \
    | awk '{print $1,$2/$3;}' \
    | $GIGGLE_ROOT/scripts/scatter.py \
        -a 0.25 \
        -o mc_obs_exp-odds.pdf \
        --y_label "MC observed/expected" \
        --x_label "Odds ratio (GIGGLE)" \
        --fig_x 3 \
        --fig_y 3

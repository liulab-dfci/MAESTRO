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

    # fill out form on http://cistrome.org/db/interface.html
    # check Human_TF Human_histone Human_chromatin_accessibility Human_other
    # get TF_humar.tar.gz

    tar -xvf TF_human.tar.gz -C data/
    cd data
    mkdir tmp
    tar -zxvf TF_human.tar.gz -C tmp/
    tar -zxvf tmp/TF_human.tar.gz
    rm -rf tmp/
    rm TF_human.tar.gz
    rm ../TF_human.tar.gz
    mkdir named
    mkdir named_sort

    cd ..

    $GIGGLE_ROOT/examples/cistrome/get_qc.sh > cistrome_id_qc.txt

    $GIGGLE_ROOT/examples/cistrome/rename.py \
        -m data/TF_human_data_information.txt \
        -i data/TF_human \
        -o data/named \
        -n cistrome_id_to_name_map.txt
        

    $GIGGLE_ROOT/scripts/sort_bed "data/named/[A-J]*" data/named_sort/ 10
    $GIGGLE_ROOT/scripts/sort_bed "data/named/[K-T]*" data/named_sort/ 10
    $GIGGLE_ROOT/scripts/sort_bed "data/named/[U-Z]*" data/named_sort/ 10

    mkdir data/named_q100_sort
    ls data/named_sort/ \
    | gargs -p 10 "$GIGGLE_ROOT/examples/cistrome/get_top.sh data/named_sort/{} data/named_q100_sort 100"

    time ~/src/giggle/bin/giggle index \
        -i "data/named_q100_sort/*gz"  \
        -o data/named_q100_sort_b \
        -s -f
    Indexed 8716024 intervals.

    real    0m17.843s
    user    0m15.791s
    sys     0m1.575s

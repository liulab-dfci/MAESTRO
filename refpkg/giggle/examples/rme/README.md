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

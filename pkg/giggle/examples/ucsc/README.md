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
    make -j 20
    export GIGGLE_ROOT=`pwd`
    cd ..

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

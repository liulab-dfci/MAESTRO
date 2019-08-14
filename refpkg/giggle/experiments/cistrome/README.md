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
    | gargs -p 10 "./get_top.sh data/named_sort/{} data/named_q100_sort 100"
    time ~/src/giggle/bin/giggle index \
        -i "data/named_q100_sort/*gz"  \
        -o data/named_q100_sort_b \
        -s -f
    Indexed 8716024 intervals.

    real    0m17.843s
    user    0m15.791s
    sys     0m1.575s

    mkdir data/named_q100_sort_result/
    ls data/named_q100_sort/ \
    | gargs -p 10 '$GIGGLE_ROOT/bin/giggle search -i data/named_q100_sort_b -q data/named_q100_sort/{} -s > data/named_q100_sort_result/{}.results'

    ls data/named_q100_sort/ \
    | gargs -o -p 10 'echo -en "{}\t"; bgzip -d -c data/named_q100_sort/{} | wc -l' \
    > named_q100_sort_num_lines.txt

    # keep 5 GSM640420,GSM588576,GSM529981,GSM365930,GSM1116655,
    POL2RA=GSM1116660,GSM1116661,GSM1116656,GSM808756,GSM1143125,GSM1006876,GSM822295,GSM1006865,GSM1091914,GSM1091915,GSM1091916,GSM1091917,GSM1091918,GSM1091919,GSM1091920,GSM1091921,GSM1276019,GSM1276020,GSM1276021,GSM1276023,GSM1276024,GSM1388123,GSM1388129,GSM1533404,GSM1533405,GSM1533406,GSM1533407,GSM1533408,GSM1533409,GSM1533410,GSM1533411,GSM1533413,GSM1636933,GSM1636934,GSM1523077,GSM1523078,GSM1523079,GSM1523080 
    # GSM614622,GSM614620,GSM614619,GSM614618,GSM614613,
    RAD21=GSM614612,GSM1010791,GSM1861937,GSM1861938,GSM1861939,GSM1861940,GSM1861941,GSM186194
    # GSM614615,GSM614614,GSM808752,GSM808753,GSM822305,
    CTCF=GSM1006875,GSM1022663,GSM1022658,GSM1006878,GSM822308,GSM822309,GSM1010734,GSM1817665,GSM1817666,GSM1817667,GSM631475,GSM631476,GSM631477,GSM631478,GSM631479
    FOXA1=GSM1534740,GSM1534741,GSM1534743,GSM798436,GSM798438
    FOXM1=GSM1000996
    # tabeled MCF-7 by cistrome, but is 
    OTHER_EX=GSM631474
    # siRNA
    OTHER_EX+=,GSM1122652,GSM1122653,GSM614621,GSM986086,GSM986085,GSM986087,GSM986088
    # Transfected
    OTHER_EX+=,GSM1861942,GSM1388124,GSM1388125,GSM1388127
    # Mutant lines
    OTHER_EX+=,GSM699989,GSM699988,GSM699987
    
    ESR1_XENOGRAPH=GSM1669134,GSM1669138,GSM1669140,GSM1669142,GSM1669143,GSM1669144,GSM1669145,GSM1019131
    ESR1_TREATMENT=GSM614610,GSM365928,GSM365927,GSM798435,GSM798424,GSM798434
    ESR1_SIRNA=GSM631465,GSM631467,GSM631468,GSM986064,GSM986063,GSM986061,GSM986059,GSM986060,GSM986062,GSM1198712,GSM1198714,GSM1295590,GSM1534746,GSM1534747,GSM1534750,GSM1534751,GSM1523081,GSM1523082,GSM1967546,GSM1967547

    
    $GIGGLE_ROOT/scripts/cross.py  \
        -i "data/named_q100_sort_result/*" \
        -c data/TF_human_data_information.txt \
        --name_map cistrome_id_to_name_map.txt \
        --qc  cistrome_id_qc.txt \
        --lc named_q100_sort_num_lines.txt \
        -d <(cat data/TF_human_data_information.txt | awk -F '\t' '$4=="MCF-7" && $7!="ESR1"' | sort -k 7 | cut -f2) \
        -q <(cat data/TF_human_data_information.txt | awk -F '\t' '$4=="MCF-7" && $7=="ESR1"' | cut -f2) \
        --x_size 30 --y_size 15 \
        --db_x $POL2RA,$RAD21,$CTCF,$FOXA1,$FOXM1,$OTHER_EX \
        --q_x $ESR1_XENOGRAPH,$ESR1_TREATMENT,$ESR1_SIRNA \
        --pretty_names $GIGGLE_ROOT/examples/cistrome/pretty_names.txt \
    -o mcf-7_esr1-all.pdf

    $GIGGLE_ROOT/scripts/cross.py  \
        -i "data/named_q100_sort_result/*" \
        -c data/TF_human_data_information.txt \
        --name_map cistrome_id_to_name_map.txt \
        --qc  cistrome_id_qc.txt \
        --lc named_q100_sort_num_lines.txt \
        -d <(cat data/TF_human_data_information.txt | awk -F '\t' '$4=="MCF-7"' | sort -k 7 | cut -f2) \
        -q <(cat data/TF_human_data_information.txt | awk -F '\t' '$4=="MCF-7"' | sort -k 7 | cut -f2) \
        --x_size 100 --y_size 100 \
    -o mcf_7_x_mcf_y.pdf

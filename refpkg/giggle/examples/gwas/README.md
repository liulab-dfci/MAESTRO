Get data

    wget https://www.nature.com/nature/journal/v518/n7539/extref/nature13835-s1.xls
    # copy first five columns into a text file named gwas.txt

    for D in `tail -n+2 gwas.txt | awk '{print $1;}'  | sort -u`; do
        grep $D gwas.txt \
        | awk '{OFS="\t"; print $4,$5,$5+1;}' \
        | bgzip -c > $D.bed.gz
    done

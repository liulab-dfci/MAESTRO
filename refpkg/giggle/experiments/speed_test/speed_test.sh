#!/bin/bash

if [ "$#" -ne "3" ]; then
    echo "usage: $0 <query file> <dataset> <genome>"
    exit 1
fi

export QUERY_FILE=$1
export DATA_DIR=$2
export GENOME=$3

echo -en "giggle\t$QUERY_FILE\t$DATA_DIR\t"
/usr/bin/time -f "%e real\t%U user\t%S sys" \
    $GIGGLE_ROOT/bin/giggle search \
    -i ${DATA_DIR}_b/ \
    -q $QUERY_FILE \
>/dev/null
#| cut -d":" -f 3 | awk '{s+=$1;}END{print s;}'

echo -en "bedtools\t$QUERY_FILE\t$DATA_DIR\t"
/usr/bin/time -f "%e real\t%U user\t%S sys" \
    $BEDTOOLS_ROOT/bin/bedtools intersect \
    -sorted \
    -g $GENOME \
    -a $QUERY_FILE \
    -b $DATA_DIR/*gz \
> /dev/null
#| wc -l

echo -en "tabix\t$QUERY_FILE\t$DATA_DIR\t"
/usr/bin/time -f "%e real\t%U user\t%S sys" \
    bash -c ' ls $DATA_DIR/*gz | xargs -I {} bash -c "$HTSLIB_ROOT/tabix -R $QUERY_FILE {}" > /dev/null'
#bash -c ' ls $DATA_DIR/*gz | xargs -I {} bash -c "$HTSLIB_ROOT/tabix -R $QUERY_FILE {}" | wc -l'

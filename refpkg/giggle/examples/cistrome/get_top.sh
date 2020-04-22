#!/bin/bash

if [ "$#" -ne "3" ]; then
    echo -e "ussage:\t$0 <in file> <out dir> <min q>\n"
    exit
fi

IN_FILE=$1
OUT_DIR=$2
N=$3

B=`basename $IN_FILE`
OUT_FILE="$OUT_DIR/$B"

echo $OUT_FILE

bgzip -c -d  $IN_FILE \
| awk -v N=$N '$9 >= N' \
|  LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n | bgzip -c \
> $OUT_FILE

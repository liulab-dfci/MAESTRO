#!/bin/sh
set -euo pipefail

if [ "$#" -lt "2" ]; then
    echo "usage: `basename $0` <input path> <output dir> <threads>"
    exit 0
fi

export IN_PATH=$1
export OUT_PATH=$2

export THREADS=1

if [ "$#" -eq "3" ]; then
    export THREADS=$3
fi


SORT()
{
    set -euo pipefail
    IN=$1

    if [ ${IN: -4} == ".bed" ]; then
        BASE=`basename $IN`
        echo $BASE
        cat $IN | LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n -k3,3n | bgzip -c > $OUT_PATH/$BASE.gz
    elif [ ${IN: -7} == ".bed.gz" ]; then
        BASE=`basename $IN`
        echo $BASE
        gzip -d -c $IN | LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n -k3,3n | bgzip -c > $OUT_PATH/$BASE
    else
        echo "SKIPPING $IN"
    fi
}

export -f SORT

ls $IN_PATH | xargs -I{} -P $THREADS bash -c "SORT {}"

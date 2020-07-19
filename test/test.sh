#!/bin/bash

TEST_PREFIX=test
STD_OUTPUT_DIR=standard_output

#---------------------- block of test case 1 --------------------
echo "1. Test support of gzipped file of scATAC_10x_PeakCount.py"
# testing support of gzipped file of scATAC_10x_PeakCount.py

TEST_PEAK_GZ=test.peak.narrowPeak.gz
TEST_BC_GZ=test.barcode.txt.gz
TEST_FRAG_GZ=test.fragments.tsv.gz
TEST1_OUTPUT_DIR=${TEST_PREFIX}_peakcount_dir
TEST1_OUTPUT_FILE=${TEST1_OUTPUT_DIR}/${TEST_PREFIX}_peak_count.h5
TEST1_STD_OUTPUT=${STD_OUTPUT_DIR}/${TEST_PREFIX}_peak_count.h5

MAESTRO scatac-peakcount --peak $TEST_PEAK_GZ --fragment $TEST_FRAG_GZ --barcode $TEST_BC_GZ --species GRCh38 --cores 1 --outprefix $TEST_PREFIX -d ${TEST1_OUTPUT_DIR}

echo "1. Check if the output is the same as the standard output"
# Check if the output file is the same as standard output
d=`h5diff ${TEST1_OUTPUT_FILE} ${TEST1_STD_OUTPUT}`
if [ -z "$d" ]; then
    echo " ... success!"
else
    echo " ... failed! The first 10 lines of difference:"
    h5diff ${TEST1_OUTPUT_FILE} ${TEST1_STD_OUTPUT} | head -10
    exit 1
fi
#---------------------- end of block -----------------------------


echo "Done!"

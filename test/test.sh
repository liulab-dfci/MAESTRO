#!/bin/bash

TEST_PREFIX=test

echo "1. test support of gzipped file of scATAC_10x_PeakCount.py"

# testing support of gzipped file of scATAC_10x_PeakCount.py

TEST_PEAK_GZ=test.peak.narrowPeak.gz
TEST_BC_GZ=test.barcode.txt.gz
TEST_FRAG_GZ=test.fragments.tsv.gz

MAESTRO scatac-peakcount --peak $TEST_PEAK_GZ --fragment $TEST_FRAG_GZ --barcode $TEST_BC_GZ --species GRCh38 --cores 1 --outprefix $TEST_PREFIX
ls -l $TEST_PREFIX

#

echo "Done!"

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


#---------------------- block of test case 2 --------------------
# echo "2. Test MAESTRO scATAC-seq pipeline for 10x-genomics data"
# # Test MAESTRO scATAC-seq pipeline for 10x-genomics data

# TEST_PREFIX=${TEST_PREFIX}_scatac

# TEST_FASTQ_DIR=$HOME/Data/atac_pbmc_500_v1_fastqs_sampling
# TEST_FASTQ_PREFIX=atac_pbmc_500_v1_downsampling_0.4
# TEST_OUTPUT_DIRDIR=${TEST_PREFIX}_pipeline_dir
# TEST_GIGGLE_DIR=$HOME/Data/giggle.all
# TEST_REFERENCE_PATH=$HOME/Data/Refdata_scATAC_MAESTRO_GRCh38_1.1.0/GRCh38_genome.fa
# TEST_WHITELIST_PATH=$HOME/Data/737K-cratac-v1.txt
# TEST_PEAK_GZ=test.peak.narrowPeak.gz
# TEST_BC_GZ=test.barcode.txt.gz
# TEST_FRAG_GZ=test.fragments.tsv.gz
# TEST1_OUTPUT_DIR=${TEST_PREFIX}_peakcount_dir
# TEST1_OUTPUT_FILE=${TEST1_OUTPUT_DIR}/${TEST_PREFIX}_peak_count.h5
# TEST1_STD_OUTPUT=${STD_OUTPUT_DIR}/${TEST_PREFIX}_peak_count.h5

# MAESTRO scatac-init --format fastq --platform 10x-genomics --species GRCh38 \
# --fastq-dir $TEST_FASTQ_DIR --fastq-prefix $TEST_FASTQ_PREFIX \
# --cores 8 --directory $TEST_OUTPUT_DIRDIR --outprefix TEST_PREFIX \
# --peak-cutoff 10 --count-cutoff 100 --frip-cutoff 0.1 --cell-cutoff 10 \
# --giggleannotation $TEST_GIGGLE_DIR \
# --fasta $TEST_REFERENCE_PATH \
# --whitelist TEST_WHITELIST_PATH \
# --annotation --method both --signature human.immune.CIBERSORT --clusterpeak --rpmodel Adjusted

# cd $TEST_OUTPUT_DIRDIR;
# snakemake --cores 8
# cd ../

# echo "2. Check if the report html exists"
# # Check if the report html exists
# TEST_HTML_PATH=$TEST_OUTPUT_DIRDIR/Result/${TEST_PREFIX}_scATAC_report.html
# if [ -e $TEST_HTML_PATH ]; then
#     echo " ... succeed!"
# else
#     echo " ... failed!"
#     exit 1
# fi
#---------------------- end of block -----------------------------


echo "Done!"

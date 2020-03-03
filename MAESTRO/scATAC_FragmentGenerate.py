# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-01-17 01:43:08
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-01-17 03:00:26

import argparse
import pysam
import time, os
from collections import defaultdict

chr_list = ['chr'+str(i) for i in list(range(1,23))] 
chr_list = chr_list + ['chrX', 'chrY']

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("-B","--bam", dest = "bamfile",default = "",help = "The merged bamfile.")
    parser.add_argument("-b","--barcodefq",dest = "barcode_fastq",default = "",help = "The barcode fastq file (for sciATAC or 10x genomics, R2)")
    parser.add_argument("-O", "--outdir", dest = "outdir",default = "",help = "The output directory.")

    return parser.parse_args()


parser = CommandLineParser()
bamfile = parser.bamfile
barcode_fq = parser.barcode_fastq
outdir = parser.outdir
fragment_file = os.path.join(outdir, "fragments.tsv")


start_time = time.time()
print("Start to read R2 barcode file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
barcodefile_in = pysam.FastxFile(barcode_fq)
read_barcode_dict = {}
for barcode in barcodefile_in:
    sequence = barcode.sequence
    read_barcode_dict[barcode.name] = sequence
end_time = time.time()
print("End", end_time-start_time)


start_time = time.time()
print("Start to read mapping bam file and generate fragment file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
bamfile_in = pysam.AlignmentFile(bamfile,"rb")
fragment_out = open(fragment_file, "w")
for read in bamfile_in:
    barcode = read_barcode_dict[read.query_name]
    flag = read.flag
    if read.reference_name in chr_list and (flag & 0x2 !=0) and (flag & 0xc == 0) and (flag & 0x900 == 0) and read.mapping_quality > 30 and read.template_length < 1000 and read.template_length > 10:
        frag_list = [read.reference_name, str(read.reference_start + 4), str(read.reference_start + read.template_length - 5), barcode]
        fragment_out.write("\t".join(frag_list) + "\n")
fragment_out.close()
end_time = time.time()
print("End:", end_time-start_time)


# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-01-16 19:44:48
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-02-28 01:02:30


import argparse
import pysam
import time, os
from collections import defaultdict

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("-b","--barcodefq",dest = "barcode_fastq",default = "",help = "The barcode fastq file (for 10x genomics, R2)")
    parser.add_argument("-B","--barcodelib",dest = "barcode_library",default = "",help = "The barcode library, in which the first column is the valid barcode combinations included in experiment.")
    parser.add_argument("-O", "--outdir", dest = "outdir",default = "",help = "The output directory.")

    return parser.parse_args()

def GenerateMimatch(seq):
    base_list = ["A", "T", "C", "G"]
    seq_mut_list = []
    for i in range(len(seq)):
        seq_mut = [seq[:i] + base + seq[i+1:] for base in base_list]
        seq_mut = list(set(seq_mut))
        seq_mut_list = seq_mut_list + seq_mut
    return(seq_mut_list)

def GenerateMimatchDict(whitelist):
    barcode_dict = defaultdict(set)
    barcode_list = open(whitelist, "r").readlines()
    barcode_list = [bc.strip() for bc in barcode_list]
    with open(whitelist, "r") as filein:
        for line in filein:
            barcode = line.strip().split("\t")[0]
            barcode_mut_list = GenerateMimatch(barcode)
            for seq in barcode_mut_list:
                barcode_dict[seq].add(barcode)
    for bc in barcode_list:
        barcode_dict[bc] = {bc}
    return(barcode_dict, barcode_list)

def main():

    parser = CommandLineParser()
    barcode_fq = parser.barcode_fastq
    barcode_lib = parser.barcode_library
    outdir = parser.outdir
    barcode_correct_file = os.path.join(outdir, "barcode_correct.txt")

    if barcode_lib == "":
        # Correct barcode
        start_time = time.time()
        print("Start to correct barcode",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
        barcode_fq_in = pysam.FastxFile(barcode_fq)
        barcode_correct_file_out = open(barcode_correct_file, "w")
        for reads in barcode_fq_in:
            barcode = reads.sequence
            name = reads.name
            outstr = name + "\t" + barcode + "\t" + barcode + "\n"
            barcode_correct_file_out.write(outstr)
        barcode_correct_file_out.close()
        end_time = time.time()
        print("End", end_time-start_time)
        
    else:
        # Expand the barcode library
        start_time = time.time()
        print("Start to read barcode library file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
        barcode_lib_dict, barcode_lib_list = GenerateMimatchDict(barcode_lib)
        end_time = time.time()
        print("End", end_time-start_time)

        # Correct barcode
        start_time = time.time()
        print("Start to correct barcode",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
        barcode_fq_in = pysam.FastxFile(barcode_fq)
        barcode_correct_file_out = open(barcode_correct_file, "w")
        for reads in barcode_fq_in:
            barcode = reads.sequence
            name = reads.name
            if barcode in barcode_lib_dict:
                for bc_value in list(barcode_lib_dict[barcode]):
                    outstr = name + "\t" + barcode + "\t" + bc_value + "\n"
                    barcode_correct_file_out.write(outstr)
            else:
                continue
        barcode_correct_file_out.close()
        end_time = time.time()
        print("End", end_time-start_time)


if __name__ == "__main__":
    main()



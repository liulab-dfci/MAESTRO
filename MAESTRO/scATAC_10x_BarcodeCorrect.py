# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-01-16 19:44:48
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2021-03-30 11:45:20


import argparse
import pysam
import time, os
from collections import defaultdict

from MAESTRO.scATAC_utility import *

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("-b","--barcodefq",dest = "barcode_fastq",default = "",help = "The barcode fastq file (for 10x genomics, R2). Support gzipped fastq file.")
    parser.add_argument("-B","--barcodelib-atac",dest = "barcode_library_atac",default = "",help = "The ATAC barcode library, in which the first column is the valid barcode combinations included in experiment. Support gzipped file.")
    parser.add_argument("--barcodelib-rna",dest = "barcode_library_rna",default = "",help = "The RNA barcode library, in which the first column is the valid barcode combinations included in experiment. Support gzipped file.")
    parser.add_argument("-O", "--outdir", dest = "outdir",default = "",help = "The output directory.")

    return parser.parse_args()

def GenerateMismatch(seq):
    base_list = ["A", "T", "C", "G"]
    seq_mut_list = []
    for i in range(len(seq)):
        seq_mut = [seq[:i] + base + seq[i+1:] for base in base_list]
        seq_mut = list(set(seq_mut))
        seq_mut_list = seq_mut_list + seq_mut
    return(seq_mut_list)

def GenerateMismatchDict(whitelist_atac, whitelist_rna = ""):
    barcode_dict = defaultdict(set)

    filein = universal_open(whitelist_atac, "rt" )
    barcode_list_atac = filein.readlines()
    filein.close()
    barcode_list_atac = [bc.strip() for bc in barcode_list_atac]

    if whitelist_rna:
        filein = universal_open(whitelist_rna, "rt" )
        barcode_list_rna = filein.readlines()
        filein.close()
        barcode_list_rna = [bc.strip() for bc in barcode_list_rna]

        for i in range(len(barcode_list_atac)):
            line_atac = barcode_list_atac[i]
            barcode_atac = line_atac.strip().split("\t")[0]
            line_rna = barcode_list_rna[i]
            barcode_rna = line_rna.strip().split("\t")[0]
            barcode_mut_list = GenerateMismatch(barcode_atac)
            for seq in barcode_mut_list:
                if len(barcode_dict[seq]) == 0:
                    barcode_dict[seq].add(barcode_rna)
                else:
                    continue
        for i in range(len(barcode_list_atac)):
            line_atac = barcode_list_atac[i]
            barcode_atac = line_atac.strip().split("\t")[0]
            line_rna = barcode_list_rna[i]
            barcode_rna = line_rna.strip().split("\t")[0]
            barcode_dict[barcode_atac] = {barcode_rna}
    else:
        for i in range(len(barcode_list_atac)):
            line_atac = barcode_list_atac[i]
            barcode_atac = line_atac.strip().split("\t")[0]
            barcode_mut_list = GenerateMismatch(barcode_atac)
            for seq in barcode_mut_list:
                if len(barcode_dict[seq]) == 0:
                    barcode_dict[seq].add(barcode_atac)
                else:
                    continue
        for bc in barcode_list_atac:
            barcode_dict[bc] = {bc}

    return barcode_dict

def main():

    parser = CommandLineParser()
    barcode_fq = parser.barcode_fastq
    barcode_lib_atac = parser.barcode_library_atac
    barcode_lib_rna = parser.barcode_library_rna
    outdir = parser.outdir
    barcode_correct_file = os.path.join(outdir, "barcode_correct.txt")

    if barcode_lib_atac == "":
        # Correct barcode
        start_time = time.time()
        print("Start to correct barcode",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
        barcode_fq_in = pysam.FastxFile(barcode_fq)
        barcode_correct_file_out = open(barcode_correct_file, "w")
        for reads in barcode_fq_in:
            barcode = reads.sequence
            outstr = barcode + "\tCB\t" + barcode + "\n"
            barcode_correct_file_out.write(outstr)
        barcode_correct_file_out.close()
        end_time = time.time()
        print("End", end_time-start_time)
        
    else:
        # Expand the barcode library
        start_time = time.time()
        print("Start to read barcode library file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
        barcode_lib_dict = GenerateMismatchDict(barcode_lib_atac, barcode_lib_rna)
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
                    outstr = barcode + "\t" + "CB" + "\t" + bc_value + "\n"
                    barcode_correct_file_out.write(outstr)
            else:
                continue
        barcode_correct_file_out.close()
        end_time = time.time()
        print("End", end_time-start_time)


if __name__ == "__main__":
    main()



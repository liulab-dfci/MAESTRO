# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-24 22:26:54
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-02-28 14:24:39

import os,sys
import time
import shutil
import multiprocessing as mp
import argparse as ap
from collections import defaultdict
from functools import partial

from MAESTRO.scATAC_utility import *
from MAESTRO.scATAC_H5Process import *

tmp = randomString()

def peakcount_parser(subparsers):
    """
    Add main function peakcount argument parsers.
    """

    workflow = subparsers.add_parser("scatac-peakcount", 
        help = "Generate peak-cell binary count matrix. ")
    group_input = workflow.add_argument_group("Input arguments")
    group_input.add_argument("--peak", dest = "peak", default = "", type = str,
        help = "Location of peak file. "
        "The peak file is BED formatted with tab seperated. "
        "The first column is chromsome, the second is chromStart, and the third is chromEnd.")
    group_input.add_argument("--fragment", dest = "fragment", default = "", type = str, 
        help = "Location of fragment.tsv file. "
        "The fragments.tsv contains one line per unique fragment, with tab-separated fields as described below. "
        "Each row has 5 columns, representing chrom, chromStart, chromEnd, barcode and duplicateCount, respectively.")   
    group_input.add_argument("--barcode", dest = "barcode", default = "", type = str, 
        help = "Location of valid cell barcode file (optional). "
        "Each line of the file represents a valid barcode. "
        "If not set, the barcodes with enough read count (> --count-cutoff) "
        "in the fragment file will be used to generate peak-cell binary count matrix.")
    group_input.add_argument("--count-cutoff", dest = "count_cutoff", default = 1000, type = int, 
        help = "Cutoff for the number of count in each cell. DEFAULT: 1000.")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")
    

    group_output = workflow.add_argument_group("Output and running arguments")
    group_output.add_argument("--cores", dest = "cores", default = 8, 
        type = int, help = "Number of cores to use. DEFAULT: 8.")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")


def filter_fragment_file(barcode_file, frag_file, count_cutoff = 1000):
    """Filter fragment file and only keep the valid barcode ones."""

    barcode_list = []
    barcode_out = defaultdict(list)
    barcode_count = defaultdict(lambda: 0)

    if barcode_file == "":
        for line in open(frag_file, 'r'):
            line = line.strip().split('\t')
            barcode_out[line[3]].append(line[0]+'\t'+line[1]+'\t'+line[2])
            barcode_count[line[3]] += int(line[4])

        for k in barcode_count:
            if barcode_count[k] >= count_cutoff:
                barcode_list.append(k)

    else:
        for line in open(barcode_file,'r'):
            barcode_list.append(line.strip())
            barcode_out[line.strip()] = []

        for line in open(frag_file, 'r'):
            line = line.strip().split('\t')
            barcode_out[line[3]].append(line[0]+'\t'+line[1]+'\t'+line[2])
    
    for k in barcode_list:
        outf = open(tmp+"/"+k,'w')
        for line in sorted(barcode_out[k]):
            print(line, file=outf)
        outf.close()

    return(barcode_list)

def bedtools_intersect(barcode, peak_bed):
    """Intersect frag file with peak file to genearate binary count output."""
    os.system("bedtools intersect -wa -a " + peak_bed + " -b " + tmp + "/" + barcode + " -u > " + tmp + "/" + barcode + ".bed")
    return(tmp + "/" + barcode + ".bed")

def merge_binary_file(peak_file, count_list, count_file, genome = 'GRCh38'):
    """Merge the intersectBed result into binary count table."""

    binary_count = {}
    for line in open(peak_file, 'r'):
        line = line.strip().split('\t')
        binary_count[line[0]+'_'+line[1]+'_'+line[2]] = [0]*len(count_list)
    
    barcodes = []
    for i in range(0,len(count_list)):
        barcodes.append(count_list[i].split("/")[-1][:-4])
        for line in open(count_list[i], 'r'):
            line = line.strip().split('\t')
            if line[0]+'_'+line[1]+'_'+line[2] in binary_count:
               binary_count[line[0]+'_'+line[1]+'_'+line[2]][i] = 1

    features = list(sorted(binary_count.keys()))
    matrix = []
    for k in features:
        matrix.append(binary_count[k])
    write_10X_h5(count_file, matrix, features, barcodes, genome = genome, datatype = 'Peaks')


def peakcount(peak, fragment, barcode, cores, count_cutoff, directory, outprefix, species):

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    os.makedirs(tmp)

    count_file = os.path.join(directory, outprefix + "_peak_count.h5")

    barcode_list = filter_fragment_file(barcode, fragment, count_cutoff)

    pool = mp.Pool(processes = int(cores))
    partial_bedtools_intersect = partial(bedtools_intersect, peak_bed = peak) 
    result = pool.map_async(partial_bedtools_intersect, barcode_list)
    pool.close()
    pool.join()

    count_list = result.get()
    merge_binary_file(peak, count_list, count_file, species)
    shutil.rmtree(tmp)


# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-28 13:57:12
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-02-29 03:02:25


import os,sys
import time
import shutil
import multiprocessing as mp
import argparse as ap
from collections import defaultdict
from functools import partial

from MAESTRO.scATAC_utility import *
from MAESTRO.scATAC_H5Process import *
from MAESTRO.scATAC_10x_PeakCount import merge_binary_file

tmp = randomString()

def CommandLineParser():
    """
    Add main function peakcount argument parsers.
    """

    parser = ap.ArgumentParser(description = "Merge microfluidic QC log files into singlecell.txt. ")

    group_input = parser.add_argument_group("Input arguments")
    group_input.add_argument("--peak", dest = "peak", default = "", type = str,
        help = "Location of peak file. "
        "The peak file is BED formatted with tab seperated. "
        "The first column is chromsome, the second is chromStart, and the third is chromEnd.")
    group_input.add_argument("--bam-dir", dest = "bam_dir", default = "", type = str, 
        help = "Directory where bam files are stored.")   
    group_input.add_argument("--barcode", dest = "barcode", default = "", type = str, 
        help = "Location of valid cell barcode file (optional). "
        "Each line of the file represents a valid barcode. "
        "If not set, the barcodes with enough read count (> --count-cutoff) "
        "in the fragment file will be used to generate peak-cell binary count matrix.")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")

    group_output = parser.add_argument_group("Output and running arguments")
    group_output.add_argument("--cores", dest = "cores", default = 8, 
        type = int, help = "Number of cores to use. DEFAULT: 8.")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")

    return parser.parse_args()


def bedtools_intersect(barcode, bam_dir, peak_bed):
    """Intersect bam file with peak file to genearate binary count output."""
    if not os.path.isfile(bam_dir + "/" + barcode +".sortedByPos.rmdp.unique.bed"):
        error(bam_dir+"/"+barcode+".sortedByPos.rmdp.unique.bed not exist!")
    else:
        os.system("bedtools intersect -wa -a " + peak_bed + " -b " + bam_dir + "/" + barcode + ".sortedByPos.rmdp.unique.bed -u > " + tmp + '/' + barcode + ".bed")
    return(tmp + "/" + barcode + ".bed")


def main():

    myparser = CommandLineParser()

    peak_file = myparser.peak
    barcode_file = myparser.barcode
    directory = myparser.directory
    outprefix = myparser.outprefix
    bam_dir = myparser.bam_dir
    cores = int(myparser.cores)
    genome = myparser.species

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    os.makedirs(tmp)

    count_file = os.path.join(directory, outprefix + "_peak_count.h5")

    barcode_list = []
    for line in open(barcode_file,'r'):
        barcode_list.append(line.strip())

    pool = mp.Pool(processes = cores)
    partial_bedtools_intersect = partial(bedtools_intersect, bam_dir = bam_dir, peak_bed = peak_file) 
    result = pool.map_async(partial_bedtools_intersect, barcode_list)
    pool.close()
    pool.join()

    count_list = result.get()
    merge_binary_file(peak_file, count_list, count_file, genome)
    shutil.rmtree(tmp)

if __name__ == "__main__":
    main()

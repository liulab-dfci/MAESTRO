# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-24 22:26:54
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-12-13 20:46:06

import os,sys
import time
import shutil
import itertools
import multiprocessing as mp
import argparse as ap
import numpy as np
import scipy.sparse as sp_sparse
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
        help = "Generate peak-cell count matrix. ")
    group_input = workflow.add_argument_group("Input arguments")
    group_input.add_argument("--peak", dest = "peak", type = str, required = True,
        help = "Location of peak file. Support gzipped file."
        "The peak file is BED formatted with tab seperated. "
        "The first column is chromsome, the second is chromStart, and the third is chromEnd.")
    group_input.add_argument("--fragment", dest = "fragment", type = str, required = True,
        help = "Location of fragments.tsv file. Support gzipped file."
        "The fragments.tsv contains one line per unique fragment, with tab-separated fields as described below. "
        "Each row has 5 columns, representing chrom, chromStart, chromEnd, barcode and duplicateCount, respectively. "
        "In MAESTOR output, the file should be fragments_corrected_count.tsv. In Cell Ranger ATAC output, the file should be fragments.tsv. ")   
    group_input.add_argument("--barcode", dest = "barcode", default = "", type = str, required = False,
        help = "Location of valid cell barcode file (optional). Support gzipped file."
        "Each line of the file represents a valid barcode. "
        "If not set, the barcodes with enough read count (> --count-cutoff) "
        "in the fragment file will be used to generate peak-cell count matrix.")
    group_input.add_argument("--binary", dest = "binary", action = "store_true",
        help = "Whether or not to generate binary peak count matrix. If set, "
        "MAESTRO will binarize the peak-cell count matrix. "
        "If not (by default), MAESTRO will generate the original peak-count matrix. ")
    group_input.add_argument("--count-cutoff", dest = "count_cutoff", default = 1000, type = int, required = False,
        help = "Cutoff for the number of count in each cell. DEFAULT: 1000.")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", required = False,
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Genome assembly of either 'GRCh38' for human or 'GRCm38' for mouse. DEFAULT: GRCh38.")
    

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
        # no barcode file
        # read fragment file to get barcode
        fhd = universal_open( frag_file, "rt" )

        for line in fhd:
            line = line.strip().split('\t')
            barcode_out[line[3]].append(line[0]+'\t'+line[1]+'\t'+line[2])
            barcode_count[line[3]] += int(line[4])

        fhd.close()

        for k in barcode_count:
            if barcode_count[k] >= count_cutoff:
                barcode_list.append(k)

    else:
        # read barcode file
        fhd = universal_open( barcode_file, "rt" )

        for line in fhd:
            barcode_list.append(line.strip())
            barcode_out[line.strip()] = []
        fhd.close()

        fhd = universal_open( frag_file, "rt" )

        for line in fhd:
            line = line.strip().split('\t')
            barcode_out[line[3]].append(line[0]+'\t'+line[1]+'\t'+line[2])
        fhd.close()

    for k in barcode_list:
        outf = open(tmp+"/"+k,'w')
        for line in sorted(barcode_out[k]):
            print(line, file=outf)
        outf.close()

    return(barcode_list)

def bedtools_intersect(barcode, peak_bed):
    """Intersect frag file with peak file to genearate the count output."""
    os.system("bedtools intersect -wa -a " + peak_bed + " -b " + tmp + "/" + barcode + " -c | awk '{if ($NF>0) print $0}' > " + tmp + "/" + barcode + ".bed")
    return(tmp + "/" + barcode + ".bed")


def generate_count_matrix(count_list, peak_list, binary = False):

    peak_count = {}
    for peak in peak_list:
        peak_count[peak] = sp_sparse.dok_matrix((1, len(count_list)), dtype=np.int8)
    
    barcodes = []
    for i in range(0,len(count_list)):
        barcodes.append(count_list[i].split("/")[-1][:-4])
        for line in open(count_list[i], 'r'):
            line = line.strip().split('\t')
            if line[0]+'_'+line[1]+'_'+line[2] in peak_count:
                if binary:
                    peak_count[line[0]+'_'+line[1]+'_'+line[2]][0, i] = 1
                else:
                    peak_count[line[0]+'_'+line[1]+'_'+line[2]][0, i] = int(line[-1])

    matrix = []
    for k in peak_list:
        matrix.append(peak_count[k])
    matrix = sp_sparse.vstack(matrix)
    
    return((matrix, barcodes))


def merge_count_file(peak_file, count_list, count_file, cores, binary, genome = 'GRCh38'):
    """Merge the intersectBed result into the count table."""
    
    peak_list = []

    fhd = universal_open(peak_file, "rt")
    for line in fhd:
        line = line.strip().split('\t')
        peak_list.append(line[0]+'_'+line[1]+'_'+line[2])
    fhd.close()
    peak_list = sorted(peak_list)

    count_list_split = []
    i = 0
    while i < len(count_list):
        if i + 1000 < len(count_list):
            count_list_split.append(count_list[i : (i + 1000)])
            i += 1000
        else:
            count_list_split.append(count_list[i : len(count_list)])
            i = len(count_list)

    pool = mp.Pool(processes = int(cores))
    partial_generate_count_matrix = partial(generate_count_matrix, peak_list = peak_list, binary = binary) 
    result = pool.map_async(partial_generate_count_matrix, count_list_split)
    pool.close()
    pool.join()

    result_list = result.get()
    matrix_list = [i[0] for i in result_list]
    barcode_list = [i[1] for i in result_list]

    matrix = sp_sparse.hstack(matrix_list)
    barcodes = list(itertools.chain.from_iterable(barcode_list))

    write_10X_h5(count_file, matrix, peak_list, barcodes, genome = genome, datatype = 'Peaks')


def peakcount(peak, fragment, barcode, binary, cores, count_cutoff, directory, outprefix, species):

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
    merge_count_file(peak, count_list, count_file, cores, binary, species)
    shutil.rmtree(tmp)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 1 19:53:16 2019

@author: Chenfei Wang
"""

import os,sys
import time
import multiprocessing as mp
from MAESTRO.scATAC_utility import *
from MAESTRO.scATAC_H5Process import *

tmp = randomString()

def bedtools_intersect(barcode):
    """Intersect bam file with peak file to genearate binary count output."""
    if not os.path.isfile(sys.argv[3]+"/"+barcode+".rmdp.bed"):
        error(sys.argv[3]+"/"+barcode+".rmdp.bed not exist!")
    else:
        os.system("bedtools intersect -wa -a " + sys.argv[1] + " -b " + sys.argv[3] + '/' + barcode + ".rmdp.bed -u > " + tmp + '/' + barcode + ".bed")
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
    write_10X_h5(count_file, matrix, features, barcodes, genome = genome, type = 'Peaks')

def main():

    peak_file = sys.argv[1]
    barcode_file = sys.argv[2]
    bam_file = sys.argv[3]
    count_file = sys.argv[4]
    cores = int(sys.argv[5])
    genome = sys.argv[6]
    
    start = time.time()
    os.system("mkdir " + tmp)
    barcode_list = []
    for line in open(barcode_file,'r'):
        barcode_list.append(line.strip())
    pool = mp.Pool(processes = cores)
    result = pool.map_async(bedtools_intersect, barcode_list)
    pool.close()
    pool.join()
    count_list = result.get()
    merge_binary_file(peak_file, count_list, count_file, genome)
    os.system("rm -rf " + tmp)
    end = time.time()
    print("Peakcount Time:", end-start)    

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 1 19:53:16 2019

@author: Chenfei Wang
"""
    
import os,sys
import time
import multiprocessing as mp
from scATAC_utility import *

tmp = randomString()

def filter_fragment_file(barcode_file, frag_file):
    """Filter fragment file and only keep the valid barcode ones."""
    barcode_list = []
    barcode_out = {}
    for line in open(barcode_file,'r'):
        barcode_list.append(line.strip())
        barcode_out[line.strip()] = []

    for line in open(frag_file, 'r'):
        line = line.strip().split('\t')
        if line[3] in barcode_out:
           barcode_out[line[3]].append(line[0]+'\t'+line[1]+'\t'+line[2])
    
    for k in barcode_out.keys():
        outf = open(tmp+"/"+k,'w')
        for line in sorted(barcode_out[k]):
            print(line, file=outf)
        outf.close()

    return(barcode_list)

def bedtools_intersect(barcode):
    """Intersect frag file with peak file to genearate binary count output."""
    os.system("bedtools intersect -wa -a " + sys.argv[1] + " -b " + tmp + "/" + barcode + " -u > " + tmp + "/" + barcode + ".bed")
    return(tmp + "/" + barcode + ".bed")

def merge_binary_file(peak_file, count_list, count_file):
    """Merge the intersectBed result into binary count table."""

    binary_count = {}
    for line in open(peak_file, 'r'):
        line = line.strip().split('\t')
        binary_count[line[0]+'_'+line[1]+'_'+line[2]] = ['0']*len(count_list)
    
    header = []
    for i in range(0,len(count_list)):
        header.append(count_list[i].split("/")[-1][:-4])
        for line in open(count_list[i], 'r'):
            line = line.strip().split('\t')
            if line[0]+'_'+line[1]+'_'+line[2] in binary_count:
               binary_count[line[0]+'_'+line[1]+'_'+line[2]][i] = '1'

    outf = open(count_file, 'w')
    print("\t".join(header), file=outf)
    for k in sorted(binary_count.keys()):
        print(k+"\t"+"\t".join(binary_count[k]), file=outf)
    outf.close()

def main():

    peak_file = sys.argv[1]
    barcode_file = sys.argv[2]
    frag_file = sys.argv[3]
    count_file = sys.argv[4]
    cores = int(sys.argv[5])
    
    start = time.time()
    os.system("mkdir " + tmp)
    barcode_list = filter_fragment_file(barcode_file, frag_file)
    pool = mp.Pool(processes = cores)
    result = pool.map_async(bedtools_intersect, barcode_list)
    pool.close()
    pool.join()
    count_list = result.get()
    merge_binary_file(peak_file, count_list, count_file)
    os.system("rm -rf " + tmp)
    end = time.time()
    print("Peakcount Time:", end-start)    

if __name__ == "__main__":
    main()

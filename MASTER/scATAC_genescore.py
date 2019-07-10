#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 1 10:24:36 2019

@author: Chenfei Wang
"""

import os, sys
import math, numpy
import time
import multiprocessing as mp
from MASTER.scATAC_utility import *

Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])
chroms = ['chrY', 'chrX', 'chrM', 'chr13', 'chr12', 'chr11', 'chr10', 'chr17',\
          'chr16', 'chr15', 'chr14', 'chr19', 'chr18', 'chr22', 'chr20',\
          'chr21', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2', 'chr1',\
          'chr9', 'chr8']

def read_peak_file(peak_file):
    """Read peak count file, store it into dict, first level key is cell name, second 
       level key is chromomsome for single cell peaks, third level value is peak."""

    cell_peak = {}
    cell_list = []
    for line in open(peak_file, 'r').readlines():
        line = line.strip().split('\t')
        peak = line[0].split('_')
        if not peak[0].startswith("chr"):
            cell_list = line
            for cell in cell_list:
                cell_peak[cell] = {}
        else:
            for i in range(0,len(cell_list)):
                if int(line[i+1]) > 0:
                   if not peak[0] in cell_peak[cell_list[i]]:
                       cell_peak[cell_list[i]][peak[0]] = [[peak[0], int(peak[1]), int(peak[2])]]
                   else:
                       cell_peak[cell_list[i]][peak[0]].append([peak[0], int(peak[1]), int(peak[2])])
    
    for cell in cell_list:
        for chr in cell_peak[cell]:
            cell_peak[cell][chr].sort()
    
    return(cell_peak, cell_list)

def RP(args):
    """Multiple processing function to calculate regulation potential."""
    
    gene_info, cell_peak, gene_distance, cell = args
    score = []
    for gene in gene_info:
         gene = gene.split('@')
         if gene[3] == "+":
             gTSS = int(gene[1])
             gTTS = int(gene[2])
         elif gene[3] == "-":
             gTSS = int(gene[2])
             gTTS = int(gene[1])
         
         try:
             peaks = cell_peak[cell][gene[0]]
         except KeyError:
             peaks = []
         peaksInDistance = [abs((t[1]+t[2])/2-gTSS)*1.0/gene_distance for t in peaks if abs((t[1]+t[2])/2-gTSS) < gene_distance]
         peaksInDistance.sort()
         if len(peaksInDistance) > 10000:  # extract no more than 10k peaks
             peaksInDistance = peaksInDistance[:10000]
         score.append(str(Sg(peaksInDistance)))         
    
    return(score)

def calculate_RP_score(cell_peak, cell_list, score_file, gene_distance, gene_bed, cores):
    """Calculate regulatery potential for each gene based on the single-cell peaks."""

    gene_info = []
    for line in open(gene_bed, 'r'):
        line = line.strip().split('\t')
        if not line[0].startswith('#'):
            if line[1] in chroms:
               gene_info.append(line[1]+'@'+line[3]+'@'+line[4]+'@'+line[2]+'@'+line[5])
               # gene_info 'chrom@start@end@strand@symbol'
    gene_info = list(set(gene_info))

    args = ((gene_info, cell_peak, gene_distance, cell) for cell in cell_list)
    pool = mp.Pool(processes = cores)
    result = pool.map_async(RP, args)
    pool.close()
    pool.join()
    score_list = result.get()    
    
    score_dup = {}                              # recode gene score with duplicate names
    for i in range(0, len(gene_info)):
        score_dup[gene_info[i]] = []
        for j in range(0,len(cell_list)):
            score_dup[gene_info[i]].append(score_list[j][i])
    
    score_nondup = {}                           # remove the potential duplciate names
    for k in score_dup.keys():
        gene = k.split('@')[4]
        if not gene in score_nondup:
           score_nondup[gene] = score_dup[k]
        else:
           score_current = numpy.mean([float(t) for t in score_dup[k]])
           score_max = numpy.mean([float(t) for t in score_nondup[gene]])
           if score_current > score_max:
              score_nondup[gene] = score_dup[k]
    
    outf = open(score_file, 'w')
    print("\t".join(cell_list), file=outf)
    for k in score_nondup.keys():
        print(k+"\t"+"\t".join(score_nondup[k]), file=outf)
    outf.close()

def main():

    peak_file = sys.argv[1]
    score_file = sys.argv[2]
    gene_distance = int(sys.argv[3])
    gene_bed = sys.argv[4]
    cores = int(sys.argv[5])

    start = time.time()
    cell_peak, cell_list = read_peak_file(peak_file)
    calculate_RP_score(cell_peak, cell_list, score_file, gene_distance, gene_bed, cores)   
    end = time.time()
    print("GeneScore Time:", end-start)

if __name__ == "__main__":
    main()
    
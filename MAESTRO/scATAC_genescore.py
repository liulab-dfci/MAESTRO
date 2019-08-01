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
from MAESTRO.scATAC_utility import *

Sg_CistromeGO = lambda ldx: sum([2**(-t) for t in ldx])
chroms = ['chrY', 'chrX', 'chrM', 'chr13', 'chr12', 'chr11', 'chr10', 'chr17',\
          'chr16', 'chr15', 'chr14', 'chr19', 'chr18', 'chr22', 'chr20',\
          'chr21', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2', 'chr1',\
          'chr9', 'chr8']

def read_peak_file(peak_file):
    """Read peak count file, store it into dict, first level key is cell name, second 
       level key is chromomsome for single cell peaks, third level value is peak."""

    cell_peak = {}
    cell_list = []
    for line in open(peak_file, 'r'):
        line = line.strip().split('\t')
        peak = line[0].split('_')
        if not peak[0] in chroms:
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
    
    gene_info, cell_peak, decay_distance, cell = args
    score = []
    score_memory = {}
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
 
         peaksInDistance = [round(abs((t[1]+t[2])/2-gTSS)*1.0/decay_distance,3) for t in peaks if abs((t[1]+t[2])/2-gTSS) < 15*decay_distance]
         peaksInDistance.sort()
         ID = '_'.join([str(t) for t in peaksInDistance])
         if ID in score_memory:
             score.append(score_memory[ID])      # save results in memory to accelerate
         else:
             Sg = Sg_CistromeGO(peaksInDistance)
             score.append(Sg)         
             score_memory[ID] = Sg

    return(score)

def calculate_RP_score(cell_peak, cell_list, score_file, decay_distance, gene_bed, rp_json, cores):
    """Calculate regulatery potential for each gene based on the single-cell peaks."""

    gene_info = []
    for line in open(gene_bed, 'r'):
        line = line.strip().split('\t')
        if not line[0].startswith('#'):
            if line[1] in chroms:
               gene_info.append(line[1]+'@'+line[3]+'@'+line[4]+'@'+line[2]+'@'+line[5])
               # gene_info 'chrom@start@end@strand@symbol'
    gene_info = list(set(gene_info))

    args = ((gene_info, cell_peak, decay_distance, cell) for cell in cell_list)
    pool = mp.Pool(processes = cores)
    result = pool.map_async(RP, args)
    pool.close()
    pool.join()
    score_list = result.get()    
    
    score_dup = {}                              # Gene score with duplicate names
    for i in range(0, len(gene_info)):
        score_dup[gene_info[i]] = []
        for j in range(0,len(cell_list)):
            score_dup[gene_info[i]].append(score_list[j][i])
    
    score_nondup = {}                           # Remove the duplicate gene names
    for k in score_dup.keys():
        gene = k.split('@')[4]
        if not gene in score_nondup:
           score_nondup[gene] = score_dup[k]
        else:
           score_current = numpy.mean([t for t in score_dup[k]])
           score_max = numpy.mean([t for t in score_nondup[gene]])
           if score_current > score_max:
              score_nondup[gene] = score_dup[k]
    
    rp_bgfile = open(rp_json, 'r')              # Normalize the RP score
    rp_bg = json.load(rp_bgfile)
    rp_bgfile.close()
    for k in score_nondup.keys():
        if k in rp_bg:
            score_nondup[k] = [str(round(max(0, t - rp_bg[k]),3)) for t in score_nondup[k]]
        else:
            score_nondup[k] =  [str(round(t,3)) for t in score_nondup[k]]
    
    outf = open(score_file, 'w')
    print("\t".join(cell_list), file=outf)
    for k in sorted(score_nondup.keys()):
        print(k+"\t"+"\t".join(score_nondup[k]), file=outf)
    outf.close()

def main():

    peak_file = sys.argv[1]
    score_file = sys.argv[2]
    decay_distance = int(sys.argv[3])
    gene_bed = sys.argv[4]
    rp_json = sys.argv[5]
    cores = int(sys.argv[6])

    start = time.time()
    cell_peak, cell_list = read_peak_file(peak_file)
    calculate_RP_score(cell_peak, cell_list, score_file, decay_distance, gene_bed, rp_json, cores)   
    end = time.time()
    print("GeneScore Time:", end-start)

if __name__ == "__main__":
    main()

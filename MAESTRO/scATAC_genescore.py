#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 1 10:24:36 2019

@author: Chenfei Wang, Changxin Wan
"""

import os, sys
import math, numpy
import time
import multiprocessing as mp
from MAESTRO.scATAC_utility import *


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
                cell_peak[cell] = []
        else:
            for i in range(0, len(cell_list)):
                if int(line[i + 1]) > 0:
                    cell_peak[cell_list[i]].append((peak[0], (int(peak[1]) + int(peak[2])) / 2.0, 0, None))

    return (cell_peak, cell_list)


def RP(args):
    """Multiple processing function to calculate regulation potential."""

    Sg = lambda ldx: sum([math.exp(-0.5 - 4 * t) for t in ldx])
    gene_info, cell_peak, gene_distance = args
    score_dict = {}

    w = gene_info + cell_peak

    D = {}
    A = {}

    w.sort()
    for elem in w:
        if elem[2] == 1:
            A[elem[-1]] = [elem[0], elem[1]]
            D[elem[-1]] = []
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                if (g[0] != elem[0]) or ((elem[1] - g[1]) > gene_distance):
                    dlist.append(gene_name)
                else:
                    A[gene_name].append(
                        (elem[1] - g[1]) / gene_distance)  # peak in distance will calculate the distance
            for gene_name in dlist:
                D[gene_name] += A.pop(gene_name)[2:]

    w.reverse()
    for elem in w:
        if elem[2] == 1:
            A[elem[-1]] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                if (g[0] != elem[0]) or (-(elem[1] - g[1]) > gene_distance):
                    dlist.append(gene_name)
                else:
                    A[gene_name].append(-(elem[1] - g[1]) / gene_distance)
            for gene_name in dlist:
                D[gene_name] += A.pop(gene_name)[2:]

    for gene_name in list(A.keys()):
        D[gene_name] += A.pop(gene_name)[2:]

    for gene in D.keys():
        score_dict[gene] = Sg(D[gene])

    return (score_dict)


def calculate_RP_score(cell_peak, cell_list, score_file, gene_distance, gene_bed, cores):
    """Calculate regulatery potential for each gene based on the single-cell peaks."""

    gene_info = []
    for line in open(gene_bed, 'r'):
        line = line.strip().split('\t')
        if not line[0].startswith('#'):
            if line[2] == "+":
                gene_info.append((line[1], int(line[3]), 1, "%s@%s@%s" % (line[5], line[1], line[3])))
            else:
                gene_info.append((line[1], int(line[4]), 1, "%s@%s@%s" % (line[5], line[1], line[3])))
                # gene_info [chrom, tss, 1, gene_unique]

    gene_info = list(set(gene_info))
    genes = list(set([i[3].split("@")[0] for i in gene_info]))

    args = ((gene_info, cell_peak[cell], gene_distance) for cell in cell_list)
    pool = mp.Pool(processes=cores)
    result = pool.map_async(RP, args)
    pool.close()
    pool.join()
    score_list = result.get()

    score_cells_dict = {}
    for gene in gene_info:
        score_cells_dict[gene[3]] = []
        for icell in range(0, len(cell_list)):
            score_cells_dict[gene[3]].append(score_list[icell][gene[3]])

    score_cells_dict_dedup = {}  # recode gene score with duplicate names
    for gene in genes:
        score_cells_dict_dedup[gene] = [0]
    for gene in score_cells_dict.keys():
        if sum(score_cells_dict[gene]) >= sum(score_cells_dict_dedup[gene.split("@")[0]]):
            score_cells_dict_dedup[gene.split("@")[0]] = score_cells_dict[gene]
        else:
            continue

    outf = open(score_file, 'w')
    print("\t".join(cell_list), file=outf)
    for k in score_cells_dict_dedup.keys():
        print(k + "\t" + "\t".join(map(str, score_cells_dict_dedup[k])), file=outf)
    outf.close()

def main():

    peak_file = sys.argv[1]
    score_file = sys.argv[2]
    gene_distance = float(sys.argv[3])
    gene_bed = sys.argv[4]
    cores = int(sys.argv[5])

    start = time.time()
    cell_peak, cell_list = read_peak_file(peak_file)
    calculate_RP_score(cell_peak, cell_list, score_file, gene_distance, gene_bed, cores)
    end = time.time()
    print("GeneScore Time:", end - start)

if __name__ == "__main__":
    main()
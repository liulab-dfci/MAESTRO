#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 1 10:24:36 2019

@author: Chenfei Wang, Changxin Wan
"""

import os, sys
import math
import time
import numpy as np
import pandas as pd
import scipy.sparse as sparse
from MAESTRO.scATAC_utility import *


def RP(peaks_info, genes_info, gene_distance):
    """Multiple processing function to calculate regulation potential."""

    Sg = lambda x: math.exp(-0.5 - 4 * x)
    genes_peaks_score_array = sparse.dok_matrix((len(genes_info), len(peaks_info)), dtype=np.float64)

    w = genes_info + peaks_info

    A = {}

    w.sort()
    for elem in w:
        if elem[2] == 1:
            A[elem[-1]] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = elem[1] - g[1]
                if (g[0] != elem[0]) or (tmp_distance > gene_distance):
                    dlist.append(gene_name)
                else:
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / gene_distance)
            for gene_name in dlist:
                del A[gene_name]

    w.reverse()
    for elem in w:
        if elem[2] == 1:
            A[elem[-1]] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = g[1] - elem[1]
                if (g[0] != elem[0]) or tmp_distance > gene_distance:
                    dlist.append(gene_name)
                else:
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / gene_distance)
            for gene_name in dlist:
                del A[gene_name]

    return(genes_peaks_score_array)


def calculate_RP_score(peak_file, gene_bed, gene_distance, score_file):
    """Calculate regulatery potential for each gene based on the single-cell peaks."""

    genes_info = []
    genes_list = []
    for line in open(gene_bed, 'r'):
        line = line.strip().split('\t')
        if not line[0].startswith('#'):
            if line[2] == "+":
                genes_info.append((line[1], int(line[3]), 1, "%s@%s@%s" % (line[5], line[1], line[3])))
            else:
                genes_info.append((line[1], int(line[4]), 1, "%s@%s@%s" % (line[5], line[1], line[4])))
                # gene_info [chrom, tss, 1, gene_unique]
    genes_info = list(set(genes_info))
    for igene in range(len(genes_info)):
        tmp_gene = list(genes_info[igene])
        genes_list.append(tmp_gene[3])
        tmp_gene[3] = igene
        genes_info[igene] = tmp_gene
    genes = list(set([i.split("@")[0] for i in genes_list]))

    peaks_info = []
    cell_peaks = pd.read_csv(peak_file, sep="\t", header=0, index_col=0)
    cells_list = list(cell_peaks.columns)
    peaks_list = list(cell_peaks.index)
    cell_peaks = sparse.csc_matrix(cell_peaks.values)
    for ipeak in range(len(peaks_list)):
        peaks_tmp = peaks_list[ipeak].split("_")
        peaks_info.append([peaks_tmp[0], (int(peaks_tmp[1]) + int(peaks_tmp[2])) / 2.0, 0, ipeak])

    genes_peaks_score_dok = RP(peaks_info, genes_info, gene_distance)
    genes_peaks_score_csc = genes_peaks_score_dok.tocsc()
    genes_cells_score_csc = genes_peaks_score_csc.dot(cell_peaks)

    score_cells_dict = {}
    score_cells_sum_dict = {}
    score_cells_sum = np.asarray(genes_cells_score_csc.sum(axis=1)).ravel().tolist()
    genes_cells_score_csr = genes_cells_score_csc.tocsr()
    for igene in range(len(genes_list)):
        score_cells_dict[genes_list[igene]] = genes_cells_score_csr[igene, :].toarray().ravel().tolist()
        score_cells_sum_dict[genes_list[igene]] = score_cells_sum[igene]

    score_cells_dict_dedup = {}
    score_cells_dict_max = {}
    for gene in genes:
        score_cells_dict_max[gene] = 0

    for gene in genes_list:
        symbol = gene.split("@")[0]
        if score_cells_sum_dict[gene] >= score_cells_dict_max[symbol]:
            score_cells_dict_dedup[symbol] = score_cells_dict[gene]
            score_cells_dict_max[symbol] = score_cells_sum_dict[gene]

    outf = open(score_file, 'w')
    outf.write("\t".join(cells_list) + "\n")
    for symbol in score_cells_dict_dedup.keys():
        outf.write(symbol + "\t" + "\t".join(map(str, score_cells_dict_dedup[symbol])) + "\n")
    outf.close()

def main():

    peak_file = sys.argv[1]
    score_file = sys.argv[2]
    gene_distance = float(sys.argv[3])
    gene_bed = sys.argv[4]
    cores = int(sys.argv[5])

    start = time.time()
    calculate_RP_score(peak_file, gene_bed, gene_distance, score_file)
    end = time.time()
    print("GeneScore Time:", end - start)

if __name__ == "__main__":
    main()
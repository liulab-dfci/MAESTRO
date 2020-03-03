#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 1 10:24:36 2019

@author: Chenfei Wang, Changxin Wan
"""

import os, sys
import time
import numpy as np
import scipy.sparse as sp_sparse
import tables
import h5py
import collections

def RP(peaks_info, genes_info, decay):
    """Multiple processing function to calculate regulation potential."""

    Sg = lambda x: 2**(-x)
    gene_distance = 15 * decay
    genes_peaks_score_array = sp_sparse.dok_matrix((len(genes_info), len(peaks_info)), dtype=np.float64)

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
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
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
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
            for gene_name in dlist:
                del A[gene_name]

    return(genes_peaks_score_array)


def calculate_RP_score(cell_peaks, peaks_list, genes_info, genes_list, decay):
    """Calculate regulatery potential for each gene based on the single-cell peaks."""

    genes_info = genes_info.values.tolist()
    genes = list(set([i.split("@")[0] for i in genes_list]))

    peaks_info = []

    for ipeak, peak in enumerate(peaks_list):
        peaks_tmp = peak.rsplit("_",maxsplit=2)
        peaks_info.append([peaks_tmp[0], (int(peaks_tmp[1]) + int(peaks_tmp[2])) / 2.0, 0, ipeak])

    genes_peaks_score_dok = RP(peaks_info, genes_info, decay)
    genes_peaks_score_csr = genes_peaks_score_dok.tocsr()
    genes_cells_score_csr = genes_peaks_score_csr.dot(cell_peaks.tocsr())
 
    # genes_peaks_score_csc = genes_peaks_score_dok.tocsc()
    # genes_cells_score_csr = genes_peaks_score_csc.dot(cell_peaks).tocsr()

    score_cells_dict = {}
    score_cells_sum_dict = {}
    for igene, gene in enumerate(genes_list):
        # score_cells_dict[gene] = genes_cells_score_csr[igene, :].toarray().ravel().tolist()
        # score_cells_sum_dict[gene] = sum(score_cells_dict[gene])
        score_cells_dict[gene] = igene
        score_cells_sum_dict[gene] = genes_cells_score_csr[igene, :].sum()


    score_cells_dict_dedup = {}
    score_cells_dict_max = {}
    for gene in genes:
        score_cells_dict_max[gene] = float("-inf")

    for gene in genes_list:
        symbol = gene.split("@")[0]
        if score_cells_sum_dict[gene] > score_cells_dict_max[symbol]:
            score_cells_dict_dedup[symbol] = score_cells_dict[gene]
            score_cells_dict_max[symbol] = score_cells_sum_dict[gene]
    gene_symbol = sorted(score_cells_dict_dedup.keys())
    matrix_row = []
    for gene in gene_symbol:
        matrix_row.append(score_cells_dict_dedup[gene])

    score_cells_matrix = genes_cells_score_csr[matrix_row, :]
    score_cells_matrix = score_cells_matrix.tocsc()

    # score_cells_matrix = []
    # for gene in gene_symbol:
    #     score_cells_matrix.append(score_cells_dict_dedup[gene])
    # score_cells_matrix = np.array(score_cells_matrix)
    # score_cells_matrix = sp_sparse.csc_matrix(score_cells_matrix)

    return((score_cells_matrix, gene_symbol))

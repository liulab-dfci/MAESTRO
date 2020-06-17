# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-06-12 04:13:17
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-06-13 03:41:17

import os, sys
import time
import tables
import h5py
import collections
import numpy as np
import scipy.sparse as sp_sparse
import argparse as ap
import pandas as pd


def ExtractGeneInfo(bed):
    """Extract gene information from gene bed data frame."""

    bed['transcript'] = list(map(lambda x: x.strip().split(".")[0], bed['name'].tolist()))
    bed["chrom"] = list(map(lambda x:"chr%s"%x, bed["chrom"]))
    bed['tss'] = bed.apply(lambda x: x['txStart'] if x['strand']=='+' else x['txEnd'], axis=1)

    ### adjacent P+GB
    bed["start"] = bed.apply(lambda x: x['txStart']-2000 if x['strand']=='+' else x['txStart'], axis=1)
    bed["end"] = bed.apply(lambda x: x['txEnd']+2000 if x['strand']=='-' else x['txEnd'], axis=1)
    
    bed['promoter'] = bed.apply(lambda x: tuple([x['tss']-2000, x['tss']+2000]), axis=1)
    bed['exons'] = bed.apply(lambda x: tuple([(int(i), int(j)) for i, j in zip(x['exonStarts'].strip(',').split(','), x['exonEnds'].strip(',').split(','))]), axis=1)

    ### exon length
    bed['length'] = bed.apply(lambda x: sum(list(map(lambda i: (i[1]-i[0])/1000.0, x['exons']))), axis=1)
    bed['uid'] = bed.apply(lambda x: "%s@%s@%s"%(x['name2'], x['start'], x['end']), axis=1)
    bed = bed.drop_duplicates(subset='uid', keep="first")
    gene_info = []
    for irow, x in bed.iterrows():
        gene_info.append([x['chrom'], x['start'], x['end'], x['tss'], x['promoter'], x['exons'], x['length'], 1, x['uid']])
    ### [chrom_0, start_1, end_2, tss_3, promoter_4, exons_5, length_6, 1_7, uid_8]
    return(gene_info)


def RP_Simple(peaks_info, genes_info, decay):
    """Multiple processing function to calculate regulation potential."""

    Sg = lambda x: 2**(-x)
    gene_distance = 15 * decay
    genes_peaks_score_array = sp_sparse.dok_matrix((len(genes_info), len(peaks_info)), dtype=np.float64)

    w = genes_info + peaks_info

    A = {}

    w.sort()
    for elem in w:
        if int(elem[2]) == 1:
            A[int(elem[-1])] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = elem[1] - g[1]
                if (g[0] != elem[0]) or (tmp_distance > gene_distance):
                    dlist.append(gene_name)
                else:
                    genes_peaks_score_array[gene_name, int(elem[-1])] = Sg(tmp_distance / decay)
            for gene_name in dlist:
                del A[gene_name]

    w.reverse()
    for elem in w:
        if int(elem[2]) == 1:
            A[int(elem[-1])] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = g[1] - elem[1]
                if (g[0] != elem[0]) or tmp_distance > gene_distance:
                    dlist.append(gene_name)
                else:
                    genes_peaks_score_array[gene_name, int(elem[-1])] = Sg(tmp_distance / decay)
            for gene_name in dlist:
                del A[gene_name]

    return(genes_peaks_score_array)


def RP_AddExon(peaks_info, genes_info_full, genes_info_tss, decay):
    """Multiple processing function to calculate regulation potential."""

    Sg = lambda x: 2**(-x)
    checkInclude = lambda x, y: all([x>=y[0], x<=y[1]])
    gene_distance = 15 * decay
    genes_peaks_score_array = sp_sparse.dok_matrix((len(genes_info_full), len(peaks_info)), dtype=np.float64)
        
    w = genes_info_tss + peaks_info
    A = {}
    
    w.sort()
    for elem in w:
        if elem[-3] == 1:
            A[elem[-1]] = elem
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = elem[1] - g[1]
                if all([g[0]==elem[0], any(list(map(checkInclude, [elem[1]]*len(g[5]), list(g[5]))))]):
                    genes_peaks_score_array[gene_name, elem[-1]] = 1.0 / g[-4]
                elif all([g[0]==elem[0], tmp_distance <= gene_distance]):
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
                else:
                    dlist.append(gene_name)
            for gene_name in dlist:
                del A[gene_name]

    w.reverse()
    for elem in w:
        if elem[-3] == 1:
            A[elem[-1]] = elem
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = g[1] - elem[1]
                if all([g[0]==elem[0], any(list(map(checkInclude, [elem[1]]*len(g[5]), list(g[5]))))]):
                    genes_peaks_score_array[gene_name, elem[-1]] = 1.0 / g[-4]
                if all([g[0]==elem[0], tmp_distance <= gene_distance]):
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
                else:
                    dlist.append(gene_name)
            for gene_name in dlist:
                del A[gene_name]
    
    return(genes_peaks_score_array)


def RP_AddExonRemovePromoter(peaks_info, genes_info_full, genes_info_tss, decay):
    """Multiple processing function to calculate regulation potential."""

    Sg = lambda x: 2**(-x)
    checkInclude = lambda x, y: all([x>=y[0], x<=y[1]])
    gene_distance = 15 * decay
    genes_peaks_score_array = sp_sparse.dok_matrix((len(genes_info_full), len(peaks_info)), dtype=np.float64)
    peaks_info_inbody = []
    peaks_info_outbody = []
    
    w = genes_info_full + peaks_info
    A = {}

    w.sort()
#     print(w[:100])
    for elem in w:
        if elem[-3] == 1:
            A[elem[-1]] = elem
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                ### NOTE: main change here
                ### if peak center in the gene area
                if all([g[0]==elem[0], elem[1]>=g[1], elem[1]<=g[2]]):
                    ### if peak center in the exons
                    if any(list(map(checkInclude, [elem[1]]*len(g[5]), list(g[5])))):
                        genes_peaks_score_array[gene_name, elem[-1]] = 1.0 / g[-4]
                        peaks_info_inbody.append(elem)
                    ### if peak cencer in the promoter
                    elif checkInclude(elem[1], g[4]):
                        tmp_distance = abs(elem[1]-g[3])
                        genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
                        peaks_info_inbody.append(elem)
                    ### intron regions
                    else:
                        continue
                else:
                    dlist.append(gene_name)
            for gene_name in dlist:
                del A[gene_name]
    
    ### remove genes in promoters and exons
    peaks_info_set = [tuple(i) for i in peaks_info]
    peaks_info_inbody_set = [tuple(i) for i in peaks_info_inbody]
    peaks_info_outbody_set = list(set(peaks_info_set)-set(peaks_info_inbody_set))
    peaks_info_outbody = [list(i) for i in peaks_info_outbody_set]
    
    # print("peaks number: ", len(peaks_info_set))
    # print("peaks number in gene promoters and exons: ", len(set(peaks_info_inbody_set)))
    # print("peaks number out gene promoters and exons:", len(peaks_info_outbody_set))
    
    w = genes_info_tss + peaks_info_outbody
    A = {}
    
    w.sort()
    for elem in w:
        if elem[-3] == 1:
            A[elem[-1]] = elem
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = elem[1] - g[1]
                if all([g[0]==elem[0], tmp_distance <= gene_distance]):
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
                else:
                    dlist.append(gene_name)
            for gene_name in dlist:
                del A[gene_name]

    w.reverse()
    for elem in w:
        if elem[-3] == 1:
            A[elem[-1]] = elem
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = g[1] - elem[1]
                if all([g[0]==elem[0], tmp_distance <= gene_distance]):
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
                else:
                    dlist.append(gene_name)
            for gene_name in dlist:
                del A[gene_name]
    
    return(genes_peaks_score_array)


def calculate_RP_score(cell_peaks, peaks_list, gene_bed_df, genes_list, decay, model):
    """Calculate regulatery potential for each gene based on the single-cell peaks."""

    # cell_peaks = pd.read_csv(peak_file, sep="\t", header=0, index_col=0)
    # cell_peaks[cell_peaks>1] = 1
    # cells_list = list(cell_peaks.columns)
    # peaks_list = [peak for peak in cell_peaks.index if peak.split("_")[1].isdigit()]
    # cell_peaks = sp_sparse.csc_matrix(cell_peaks.loc[peaks_list, :].values)
    if model == "Simple":
        peaks_info = []
        genes_info = gene_bed_df.values.tolist()
        genes = list(set([i.split("@")[0] for i in genes_list]))
        for ipeak, peak in enumerate(peaks_list):
            peaks_tmp = peak.rsplit("_",maxsplit=2)
            peaks_info.append([peaks_tmp[0], (int(peaks_tmp[1]) + int(peaks_tmp[2])) / 2.0, 0, ipeak])

        genes_peaks_score_dok = RP_Simple(peaks_info, genes_info, decay)
    else:
        genes_list = []
        genes_info = ExtractGeneInfo(gene_bed_df)
        genes_info_tss = list()
        genes_info_full = list() ### [chrom, tss, start, end, 1, unique_id]
        
        for igene in range(len(genes_info)):
            tmp_gene = genes_info[igene]
            genes_list.append(tmp_gene[-1])
            genes_info_full.append(tmp_gene + [igene])
            genes_info_tss.append([tmp_gene[0], tmp_gene[3], tmp_gene[1], tmp_gene[2]] + tmp_gene[4:] + [igene])
            ### add index at the end of gene symbol
        genes = list(set([i.split("@")[0] for i in genes_list]))
        
        # print("genes_info_full", genes_info_full[:2])
        # print("genes_info_tss", genes_info_tss[:2])
        
        peaks_info = []
        
        for ipeak, peak in enumerate(peaks_list):
            peaks_tmp = peak.rsplit("_", maxsplit=2)
            peaks_info.append([peaks_tmp[0], (int(peaks_tmp[1])+int(peaks_tmp[2]))/2.0, int(peaks_tmp[1]), int(peaks_tmp[2]), 0, peak, ipeak])
            # peaks_info [chrom, center, start, end, 0, uid, ipeak]
        # print("peaks_info", peaks_info[:2])
        if model == "Adjusted":
            genes_peaks_score_dok = RP_AddExonRemovePromoter(peaks_info, genes_info_full, genes_info_tss, decay)


    genes_peaks_score_csr = genes_peaks_score_dok.tocsr()
    genes_cells_score_csr = genes_peaks_score_csr.dot(cell_peaks.tocsr())

    # genes_peaks_score_csc = genes_peaks_score_dok.tocsc()
    # genes_cells_score_csr = genes_peaks_score_csc.dot(cell_peaks).tocsr()
    
    # genes_cells_score_lil = genes_cells_score_csc.tolil()

    score_cells_dict = {}
    score_cells_sum_dict = {}
    # for icell in range(len(cells_list)):
    #     genes_cells_score_lil[:, icell] = np.array(normMinMax(genes_cells_score_lil[:, icell].toarray().ravel().tolist())).reshape((len(genes_list), 1))
    # genes_cells_score_csr = genes_cells_score_lil.tocsr()
    for igene, gene in enumerate(genes_list):
        # score_cells_dict[gene] = list(map(lambda x: x - bgrp_dict[gene], genes_cells_score_csr[igene, :].toarray().ravel().tolist()))
        # score_cells_dict[gene] = genes_cells_score_csr[igene, :].toarray().ravel().tolist()
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
    
    return((score_cells_matrix, gene_symbol))


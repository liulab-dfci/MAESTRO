# -*- coding: utf-8 -*-
# @Author: Chenfei Wang, Changxin Wan
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-23 19:48:03
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-07-26 16:04:49


import os, sys
import time
import tables
import h5py
import re
import collections
import numpy as np
import scipy.sparse as sp_sparse
import argparse as ap
import pandas as pd

from pkg_resources import resource_filename

from MAESTRO.scATAC_utility import *
from MAESTRO.scATAC_H5Process import *

def genescore_parser(subparsers):
    """
    Add main function init-scatac argument parsers.
    """

    workflow = subparsers.add_parser("scatac-genescore", 
        help = "Calculate gene score according to scATAC peak count.")
    group_input = workflow.add_argument_group("Input arguments")
    group_input.add_argument("--format", dest = "format", default = "", 
        choices = ["h5", "mtx", "plain"], 
        help = "Format of the count matrix file.")
    group_input.add_argument("--peakcount", dest = "peakcount", default = "", 
        help = "Location of peak count matrix file. "
        "Peak count matrix with peaks as rows and cells as columns. "
        "If the format is 'h5' or 'plain', users need to specify the name of the count matrix file "
        "and row names should be like 'chromosome_peakstart_peakend', such as 'chr10_100020591_100020841'. "
        "If the format is 'mtx', the 'matrix' should be the name of .mtx formatted matrix file, such as 'matrix.mtx'.")
    group_input.add_argument("--feature", dest = "feature", default = "", 
        help = "Location of feature file (required for the format of 'mtx'). "
        "Features correspond to row indices of count matrix. "
        "The feature file should be the peak bed file with 3 columns. For example, peaks.bed.")
    group_input.add_argument("--barcode", dest = "barcode", default = "", 
        help = "Location of barcode file (required for the format of 'mtx'). "
        "Cell barcodes correspond to column indices of count matrix. For example, barcodes.tsv. ")
    group_input.add_argument("--genedistance", dest = "genedistance", default = 10000, type = int, 
        help = "Gene score decay distance, could be optional from 1kb (promoter-based regulation) "
        "to 10kb (enhancer-based regulation). DEFAULT: 10000.")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")
    group_input.add_argument("--model", dest = "model", default = "Enhanced", 
        choices = ["Simple", "Enhanced"], type = str, 
        help = "The RP model to use to calaculate gene score. "
        "For each gene, simple model sums over the impact of all regulatory elements within the up/dowm-stream of TSS. "
        "On the basis of simple model, enhanced model gives the regulatory elements within the exon region a higher weight, "
        "and also excludes the regulatory elements overlapped with another gene (the promoter and exon of a nearby gene). "
        "See the MAESTRO paper for more details. DEFAULT: Enhanced.")

    group_output = workflow.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")


def ExtractGeneInfo(gene_bed):
    """Extract gene information from gene bed file."""

    bed = pd.read_csv(gene_bed, sep="\t", header=0, index_col=False)
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
    
    print("peaks number: ", len(peaks_info_set))
    print("peaks number in gene promoters and exons: ", len(set(peaks_info_inbody_set)))
    print("peaks number out gene promoters and exons:", len(peaks_info_outbody_set))
    
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


def calculate_RP_score(peakmatrix, features, barcodes, gene_bed, decay, score_file, model):
    """Calculate regulatery potential for each gene based on the single-cell peaks."""

    genes_info = []
    genes_list = []

    peaks_info = []
    cell_peaks = peakmatrix
    peaks_list = features
    cells_list = barcodes
    # cell_peaks = pd.read_csv(peak_file, sep="\t", header=0, index_col=0)
    # cell_peaks[cell_peaks>1] = 1
    # cells_list = list(cell_peaks.columns)
    # peaks_list = [peak for peak in cell_peaks.index if peak.split("_")[1].isdigit()]
    # cell_peaks = sp_sparse.csc_matrix(cell_peaks.loc[peaks_list, :].values)
    if model == "Simple":
        fhd = universal_open(gene_bed, 'rt')
        for line in fhd:
            line = line.strip().split('\t')
            if not line[0].startswith('#'):
                if line[3] == "+":
                    genes_info.append(("chr" + line[2], int(line[4]), 1, "%s@%s@%s" % (line[12], "chr" + line[2], line[4])))
                else:
                    genes_info.append(("chr" + line[2], int(line[5]), 1, "%s@%s@%s" % (line[12], "chr" + line[2], line[5])))
                    # gene_info [chrom, tss, 1, gene_unique]
        fhd.close()
        genes_info = list(set(genes_info))
        for igene in range(len(genes_info)):
            tmp_gene = list(genes_info[igene])
            genes_list.append(tmp_gene[3])
            tmp_gene[3] = igene
            genes_info[igene] = tmp_gene
        genes = list(set([i.split("@")[0] for i in genes_list]))

        for ipeak, peak in enumerate(peaks_list):
            peaks_tmp = peak.decode().rsplit("_",maxsplit=2)
            peaks_info.append([peaks_tmp[0], (int(peaks_tmp[1]) + int(peaks_tmp[2])) / 2.0, 0, ipeak])

        genes_peaks_score_dok = RP_Simple(peaks_info, genes_info, decay)
    else:
        genes_info = ExtractGeneInfo(gene_bed)
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
            peaks_tmp = peak.decode().rsplit("_", maxsplit=2)
            peaks_info.append([peaks_tmp[0], (int(peaks_tmp[1])+int(peaks_tmp[2]))/2.0, int(peaks_tmp[1]), int(peaks_tmp[2]), 0, peak, ipeak])
            # peaks_info [chrom, center, start, end, 0, uid, ipeak]
        # print("peaks_info", peaks_info[:2])
        ### change here
        # if model == "Exon+":
        #     genes_peaks_score_dok = RP_AddExon(peaks_info, genes_info_full, genes_info_tss, decay)
        if model == "Enhanced":
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
    # score_cells_matrix = []
    # for gene in gene_symbol:
    #     score_cells_matrix.append(score_cells_dict_dedup[gene])

    write_10X_h5(score_file, score_cells_matrix, gene_symbol, cells_list, genome=gene_bed.split("/")[-1].split("_")[0], datatype="Gene")

    # outf = open(score_file, 'w')
    # outf.write("\t".join(cells_list) + "\n")
    # for symbol in score_cells_dict_dedup.keys():
    #     outf.write(symbol + "\t" + "\t".join(map(str, score_cells_dict_dedup[symbol])) + "\n")
    # outf.close()

def genescore(fileformat, directory, outprefix, peakcount, feature, barcode, genedistance, species, model = "Enhanced"):

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass
    
    annotation_path = resource_filename('MAESTRO', 'annotations')
    # annotation_path = os.path.join(os.path.dirname(__file__), 'annotations')
    genebed = os.path.join(annotation_path, species + "_refgenes.txt")
    decay = float(genedistance)
    score_file = os.path.join(directory, outprefix + "_gene_score.h5")

    if fileformat == "plain":
        matrix_dict = read_count(peakcount)
        peakmatrix = matrix_dict["matrix"]
        peakmatrix = sp_sparse.csc_matrix(peakmatrix, dtype=np.int8)
        features = matrix_dict["features"]
        features = [f.encode() for f in features]
        barcodes = matrix_dict["barcodes"]

    elif fileformat == "h5":
        scatac_count = read_10X_h5(peakcount)
        peakmatrix = scatac_count.matrix
        features = scatac_count.names.tolist()
        features = [re.sub("\W", "_", feature.decode()) for feature in features]
        features = [feature.encode() for feature in features]
        barcodes = scatac_count.barcodes.tolist()

    elif fileformat == "mtx":
        matrix_dict = read_10X_mtx(matrix_file = peakcount, feature_file = feature, barcode_file = barcode, datatype = "Peak")
        peakmatrix = matrix_dict["matrix"]
        features = matrix_dict["features"]
        features = [f.encode() for f in features]
        barcodes = matrix_dict["barcodes"]

    calculate_RP_score(peakmatrix, features, barcodes, genebed, decay, score_file, model)

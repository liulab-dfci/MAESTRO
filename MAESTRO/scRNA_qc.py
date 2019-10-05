#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 11:01:43 2019

@author: Dongqing Sun
"""

import sys, os
import collections
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import scipy.io as sp_io
from MAESTRO.scATAC_H5Process import *


# FeatureBCMatrix = collections.namedtuple('FeatureBCMatrix', ['feature_ids', 'feature_names', 'barcodes', 'matrix'])

# def get_matrix_from_h5(filename):
#     with h5py.File(filename) as f:
#         if u'version' in f.attrs:
#             if f.attrs['version'] > 2:
#                 raise ValueError('Matrix HDF5 file format version (%d) is an newer version that is not supported by this function.' % version)
#         else:
#             raise ValueError('Matrix HDF5 file format version (%d) is an older version that is not supported by this function.' % version)
        
#         feature_ids = [x.decode('ascii', 'ignore') for x in f['matrix']['features']['id']]
#         feature_names = [x.decode('ascii', 'ignore') for x in f['matrix']['features']['name']]        
#         barcodes = list(f['matrix']['barcodes'][:])
#         matrix = sp_sparse.csc_matrix((f['matrix']['data'], f['matrix']['indices'], f['matrix']['indptr']), shape=f['matrix']['shape'])
#         return FeatureBCMatrix(feature_ids, feature_names, barcodes, matrix)

def FilterCell(rawmatrix, count_cutoff, gene_cutoff, outprefix, platform, genome):
    if platform == "10x-genomics":
        expmatrix = rawmatrix.matrix
        raw_feature_ids = rawmatrix.feature_ids
        raw_feature_names = rawmatrix.feature_names
        raw_barcodes = rawmatrix.barcodes

        count_per_cell = np.asarray(expmatrix.sum(axis=0))
        genes_per_cell = np.asarray((expmatrix > 0).sum(axis=0))
        count_gene = np.concatenate((count_per_cell,genes_per_cell), axis=0)
        np.savetxt(outprefix + "_count_gene_stat.txt", count_gene.T, delimiter=" ", fmt = "%d")
        
        passed_cell = np.logical_and(count_per_cell > count_cutoff, genes_per_cell > gene_cutoff)
        gene = [True]*expmatrix.shape[0]
        passed_cell_matrix = expmatrix[np.ix_(gene, passed_cell.tolist()[0])]

        passed_barcodes = np.array(raw_barcodes)[passed_cell.tolist()[0]].tolist()
        passed_barcodes = [bc.decode('utf-8') for bc in passed_barcodes]

        write_10X_h5(outprefix + "_filtered_gene_count_matrix.h5", matrix = passed_cell_matrix, features = raw_feature_names, barcodes = passed_barcodes, genome = genome, type = 'Gene Expression')
        # sp_io.mmwrite(outprefix + "_filtered_feature_bc_matrix/matrix.mtx", passed_cell_matrix, field = "integer")

        # passed_barcodes = np.array(raw_barcodes)[passed_cell.tolist()[0]].tolist()
        # passed_barcodes = [bc.decode('utf-8') for bc in passed_barcodes]
        # outbarcode = open(outprefix + "_filtered_feature_bc_matrix/barcodes.tsv", "w")
        # outbarcode.write("\n".join(passed_barcodes))
        # outbarcode.close()

        # feature_df = pd.DataFrame({"id":raw_feature_ids, "name":raw_feature_names, "type":"Gene Expression"})
        # feature_df.to_csv(outprefix + "_filtered_feature_bc_matrix/features.tsv", index = False, header = False, sep = "\t")
        # os.system("gzip " + outprefix + "_filtered_feature_bc_matrix/*")

    else:
        expmatrix = np.array(rawmatrix)
        count_per_cell = np.asarray(expmatrix.sum(axis=0))
        genes_per_cell = np.asarray((expmatrix > 0).sum(axis=0))
        count_gene = np.array([count_per_cell,genes_per_cell])
        np.savetxt(outprefix + "_count_gene_stat.txt", count_gene.T, delimiter=" ", fmt = "%d")

        passed_cell = np.logical_and(count_per_cell > count_cutoff, genes_per_cell > gene_cutoff)
        gene = [True]*expmatrix.shape[0]
        passed_cell_matrix = expmatrix[np.ix_(gene, passed_cell.tolist())]
        passed_cell_name = rawmatrix.columns[passed_cell]

        write_10X_h5(outprefix + "_filtered_gene_count_matrix.h5", matrix = passed_cell_matrix, features = rawmatrix.index, barcodes = passed_cell_name, genome = genome, type = 'Gene Expression')

        # passed_cell_matrix_df = pd.DataFrame(passed_cell_matrix)
        # passed_cell_matrix_df.columns = passed_cell_name
        # passed_cell_matrix_df.index = rawmatrix.index

        # passed_cell_matrix_df.to_csv(outprefix + "_filtered_gene_count_matrix.txt", index = True, header = True, sep = "\t")


    # filtered_cell = np.logical_not(passed_cell)
    
    # plt.figure(figsize=(4, 4))
    # p1 = plt.scatter(log_count_per_cell[passed_cell], genes_per_cell[passed_cell], linewidths=0, s=2, c="r")
    # p2 = plt.scatter(log_count_per_cell[filtered_cell], genes_per_cell[filtered_cell], linewidths=0, s=2, c="b")
    # plt.legend(handles=[p1,p2], labels=["cell", "non-cells"],loc = "upper left", edgecolor="none")
    # plt.xlabel("Count (log10)")
    # plt.ylabel("Genes covered")
    # plt.subplots_adjust(left = 0.2, right = 0.95, bottom = 0.15, top = 0.89)
    # plt.savefig(outprefix + "_count_gene.png", dpi=300)


def main():
    
    raw_matrix_file = sys.argv[1]
    platform = sys.argv[2]
    outpre = sys.argv[3]
    outdir = sys.argv[4]
    genome = sys.argv[5]
    
    if platform == "10x-genomics":
        raw_matrix = get_matrix_from_h5(raw_matrix_file)
        FilterCell(raw_matrix, 1000, 500, os.path.join(outdir, outpre), platform)
    elif platform == "Smartseq2":
        raw_matrix_df = pd.read_csv(raw_matrix_file, sep = "\t", header = 0, index_col = 0)
        df_rownames = list(raw_matrix_df.index)
        df_rownames = [i.split('_')[1] for i in df_rownames]
        raw_matrix_df.index = df_rownames
        df_colnames = list(raw_matrix_df.columns)
        df_colnames = [i.split('/')[-1].split('.genes.results')[0] for i in df_colnames]
        raw_matrix_df.columns = df_colnames
        FilterCell(raw_matrix_df, 1000, 500, os.path.join(outdir, outpre), platform, genome)
    else:
        raw_matrix_df = pd.read_csv(raw_matrix_file, sep = "\t", header = 0, index_col = 0)
        FilterCell(raw_matrix_df, 1000, 500, os.path.join(outdir, outpre), platform, genome)


if __name__ == "__main__":
    main()

    
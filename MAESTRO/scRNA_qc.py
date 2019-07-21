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
import matplotlib.pyplot as plt


FeatureBCMatrix = collections.namedtuple('FeatureBCMatrix', ['feature_ids', 'feature_names', 'barcodes', 'matrix'])

def get_matrix_from_h5(filename):
    with h5py.File(filename) as f:
        if u'version' in f.attrs:
            if f.attrs['version'] > 2:
                raise ValueError('Matrix HDF5 file format version (%d) is an newer version that is not supported by this function.' % version)
        else:
            raise ValueError('Matrix HDF5 file format version (%d) is an older version that is not supported by this function.' % version)
        
        feature_ids = [x.decode('ascii', 'ignore') for x in f['matrix']['features']['id']]
        feature_names = [x.decode('ascii', 'ignore') for x in f['matrix']['features']['name']]        
        barcodes = list(f['matrix']['barcodes'][:])
        matrix = sp_sparse.csc_matrix((f['matrix']['data'], f['matrix']['indices'], f['matrix']['indptr']), shape=f['matrix']['shape'])
        return FeatureBCMatrix(feature_ids, feature_names, barcodes, matrix)

def CountGenePlot(expmatrix, count_cutoff, gene_cutoff, outprefix):
    count_per_cell = np.asarray(expmatrix.sum(axis=0))
    genes_per_cell = np.asarray((expmatrix > 0).sum(axis=0))
    log_count_per_cell = np.log10(count_per_cell+1)
    count_gene = np.concatenate((count_per_cell,genes_per_cell), axis=0)
    np.savetxt(outprefix + "_count_gene_stat.txt", count_gene.T, delimiter=" ", fmt = "%d")
    
    # passed_cell = np.logical_and(count_per_cell > count_cutoff, genes_per_cell > gene_cutoff)
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
    
    if platform == "10xGenomics":
        raw_matrix = get_matrix_from_h5(raw_matrix_file).matrix
    else:
        raw_matrix_df = pd.read_csv(raw_matrix_file, sep = "\t", header = 0, index_col = 0)
        raw_matrix = np.array(raw_matrix_df)
    CountGenePlot(raw_matrix, 500, 200, os.path.join(outdir, outpre))


if __name__ == "__main__":
    main()

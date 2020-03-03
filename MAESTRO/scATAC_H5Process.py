# -*- coding: utf-8 -*-
# @Author: Chenfei Wang
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-23 19:44:05
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-02-28 01:15:51


import os
import collections
import tables
import h5py
import numpy
import scipy.io
import csv
import gzip
import scipy.sparse as sp_sparse
import argparse as ap
import pandas as pd


def mtxtoh5_parser(subparsers):
    """
    Add main function init-scatac argument parsers.
    """

    workflow = subparsers.add_parser("mtxtoh5", 
        help = "Convert 10X mtx format matrix to HDF5 format.")
    group_input = workflow.add_argument_group("Input files arguments")
    group_input.add_argument("--type", dest = "datatype", default = "Peak", 
        choices = ["Peak", "Gene"], 
        help = "Type of the count matrix (Peak for scATAC-seq and Gene for scRNA-seq). DEFAULT: Peak.")
    group_input.add_argument("--matrix", dest = "matrix", default = "matrix.mtx", 
        help = "Location of .mtx formatted matrix file. DEFAULT: matrix.mtx.")
    group_input.add_argument("--feature", dest = "feature", default = "features.tsv", 
        help = "Location of feature file (required for the format of 'mtx'). Features correspond to row indices of count matrix. "
        "If the type is Peak, please provide the peak bed file with 3 columns. "
        "If the type is Gene, each row should include gene name. DEFAULT: features.tsv.")
    group_input.add_argument("--gene-column", dest = "gene_column", default = 2, type = int,
        help = "If the type is 'Peak', please specify which column of the feature file to use for gene names. DEFAULT: 2.")
    group_input.add_argument("--barcode", dest = "barcode", default = "barcodes.tsv", 
        help = "Location of barcode file. Cell barcodes correspond to column indices of count matrix. DEFAULT: barcodes.tsv. ")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")
    

    group_output = workflow.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")


def counttoh5_parser(subparsers):
    """
    Add main function init-scatac argument parsers.
    """

    workflow = subparsers.add_parser("counttoh5", 
        help = "Convert plain text count table to HDF5 format.")
    group_input = workflow.add_argument_group("Input files arguments")
    group_input.add_argument("--type", dest = "datatype", default = "Peak", 
        choices = ["Peak", "Gene"], 
        help = "Type of the count matrix (Peak for scATAC-seq and Gene for scRNA-seq). DEFAULT: Peak.")
    group_input.add_argument("--count", dest = "count", default = "", 
        help = "Location of plain text count table file. "
        "The count file need to be tab-separated with column and row names. "
        "If the type is Peak, the row name should be like 'chromosome_peakstart_peakend'. "
        "For example, 'chr10_100020591_100020841'.")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")
    

    group_output = workflow.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")


def mergeh5_parser(subparsers):
    """
    Add main function init-scatac argument parsers.
    """

    workflow = subparsers.add_parser("mergeh5", 
        help = "Merge 10X HDF5 files.")
    group_input = workflow.add_argument_group("Input files arguments")
    group_input.add_argument("--type", dest = "datatype", default = "Peak", 
        choices = ["Peak", "Gene"], 
        help = "Type of the count matrix (Peak for scATAC-seq and Gene for scRNA-seq). DEFAULT: Peak.")
    group_input.add_argument("--h5", dest = "h5_list", default = [], nargs = "+",
        help = "Location of HDF5 count files. Multiple files should be separated by space. "
        "For example, --h5 A.h5 B.h5 C.h5 ")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")
    

    group_output = workflow.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")


FeatureBCMatrix = collections.namedtuple('FeatureBCMatrix', ['ids', 'names', 'barcodes', 'matrix'])


def read_10X_h5(filename):
    """Read 10X HDF5 files, support both gene expression and peaks."""
    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, 'matrix')
        except tables.NoSuchNodeError:
            print("Matrix group does not exist in this file.")
            return None
        feature_group = getattr(group, 'features')
        ids = getattr(feature_group, 'id').read()
        names = getattr(feature_group, 'name').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return FeatureBCMatrix(ids, names, barcodes, matrix)


def write_10X_h5(filename, matrix, features, barcodes, genome = 'GRCh38', datatype = 'Peak'):
    """Write 10X HDF5 files, support both gene expression and peaks."""
    f = h5py.File(filename, 'w')
    if datatype == 'Peak':
       M = sp_sparse.csc_matrix(matrix, dtype=numpy.int8)
    else:
       M = sp_sparse.csc_matrix(matrix, dtype=numpy.float32)
    B = numpy.array(barcodes, dtype='|S200')
    P = numpy.array(features, dtype='|S100')
    GM = numpy.array([genome]*len(features), dtype='|S10')
    FT = numpy.array([datatype]*len(features), dtype='|S100')
    AT = numpy.array(['genome'], dtype='|S10')
    mat = f.create_group('matrix')
    mat.create_dataset('barcodes', data=B)
    mat.create_dataset('data', data=M.data)
    mat.create_dataset('indices', data=M.indices)
    mat.create_dataset('indptr', data=M.indptr)
    mat.create_dataset('shape', data=M.shape)
    fet = mat.create_group('features')
    fet.create_dataset('_all_tag_keys', data=AT)
    fet.create_dataset('feature_type', data=FT)
    fet.create_dataset('genome', data=GM)
    fet.create_dataset('id', data=P)
    fet.create_dataset('name', data=P)
    f.close()

def merge_10X_h5(directory, outprefix, h5list, genome = 'GRCh38', datatype = 'Gene'):
    """Merge 10X HDF5 files, h5 filenames should be provided as list."""

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    if datatype == "Peak":
        filename = os.path.join(directory, outprefix + "_peak_count.h5")
    else:
        filename = os.path.join(directory, outprefix + "_gene_count.h5")
    
    mlist = []
    for file in h5list:
        mlist.append(read_10X_h5(file))
    
    dflist = []
    for i in range(0,len(mlist)):
        barcode_i = numpy.array([h5list[i].strip(".h5")+"@"+t.decode('UTF-8') for t in mlist[i].barcodes.tolist()], dtype='|S200')
        dflist.append(pd.DataFrame(mlist[i].matrix.toarray(), index = mlist[i].ids, columns = barcode_i))

    dfmerge = pd.concat(dflist, axis = 1, join = "outer")
    dfmerge_numpy = dfmerge.fillna(0).to_numpy()
    features = dfmerge.index.tolist()
    barcodes = dfmerge.columns.tolist()
        
    write_10X_h5(filename, dfmerge_numpy, features, barcodes, genome, datatype)


def read_10X_mtx(matrix_file, feature_file, barcode_file, datatype, gene_column = 2):
    """Convert 10x mtx as matrix."""

    matrix = scipy.io.mmread(matrix_file)
    matrix = sp_sparse.csc_matrix(matrix, dtype=numpy.int32)

    if feature_file.split('.')[-1] == 'gz' or feature_file.split('.')[-1] == 'gzip':
        feature_in = gzip.open(feature_file, "r")
    else:
        feature_in = open(feature_file, "r")
    features = feature_in.readlines() 
    if datatype == "Peak":
        features = ["_".join(feature.strip().split("\t")[0:3]) for feature in features]
    else:
        if type(features[0]) == str:
            features = [feature.strip().split("\t")[gene_column-1] for feature in features]
        if type(features[0]) == bytes:
            features = [feature.decode().strip().split("\t")[gene_column-1] for feature in features]

    if barcode_file.split('.')[-1] == 'gz' or barcode_file.split('.')[-1] == 'gzip':
        barcode_in = gzip.open(barcode_file, "r")
    else:
        barcode_in = open(barcode_file, "r")
    barcodes = barcode_in.readlines() 
    if type(barcodes[0]) == str:
        barcodes = [barcode.strip().split("\t")[0] for barcode in barcodes]
    if type(barcodes[0]) == bytes:
        barcodes = [barcode.decode().strip().split("\t")[0] for barcode in barcodes]

    return {"matrix": matrix, "features": features, "barcodes": barcodes}


def mtx_2_h5(directory, outprefix, matrix_file, feature_file, barcode_file, gene_column = 2, genome = 'GRCh38', datatype = 'Peak'):
    """Convert 10x mtx format matrix to HDF5."""

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    if datatype == "Peak":
        filename = os.path.join(directory, outprefix + "_peak_count.h5")
    else:
        filename = os.path.join(directory, outprefix + "_gene_count.h5")

    matrix_dict = read_10X_mtx(matrix_file = matrix_file, feature_file = feature_file, barcode_file = barcode_file, datatype = datatype, gene_column = gene_column)

    write_10X_h5(filename = filename, matrix = matrix_dict["matrix"], 
        features = matrix_dict["features"], barcodes = matrix_dict["barcodes"], genome = genome, datatype = datatype)


def read_count(count_file):
    """Read count table as matrix."""

    infile = open(count_file, 'r').readlines()
    barcodes = infile[0].strip().split('\t')
    features = []
    matrix = []
    for line in infile[1:]:
        line = line.strip().split('\t')
        features.append(line[0])
        matrix.append([float(t) for t in line[1:]])

    return {"matrix": matrix, "features": features, "barcodes": barcodes}


def count_2_h5(directory, outprefix, count_file, genome = 'GRCh38', datatype = 'Peak'):
    """Convert plain text count table to HDF5."""

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    if datatype == "Peak":
        filename = os.path.join(directory, outprefix + "_peak_count.h5")
    else:
        filename = os.path.join(directory, outprefix + "_gene_count.h5")

    matrix_dict = read_count(count_file)


    write_10X_h5(filename = filename, matrix = matrix_dict["matrix"], 
        features = matrix_dict["features"], barcodes = matrix_dict["barcodes"], genome = genome, datatype = datatype)


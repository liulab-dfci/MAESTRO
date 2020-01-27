import collections
import scipy.sparse as sp_sparse
import tables
import h5py
import numpy
import scipy.io
import csv

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
 
def write_10X_h5(filename, matrix, features, barcodes, genome = 'GRCh38', type = 'Peaks'):
    """Write 10X HDF5 files, support both gene expression and peaks."""
    f = h5py.File(filename, 'w')
    if type == 'Peaks':
       M = sp_sparse.csc_matrix(matrix, dtype=numpy.int8)
    else:
       M = sp_sparse.csc_matrix(matrix, dtype=numpy.float32)
    B = numpy.array(barcodes, dtype='|S100')
    P = numpy.array(features, dtype='|S100')
    GM = numpy.array([genome]*len(features), dtype='|S10')
    FT = numpy.array([type]*len(features), dtype='|S100')
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

def merge_10X_h5(filename, h5list, genome = 'GRCh38', type = 'Peaks'):
    """Merge 10X HDF5 files, h5 filenames should be provided as list."""
    mlist = []
    for file in h5list:
        mlist.append(read_10X_h5(file))
    
    features = mlist[0].ids
    barcodes = []
    matrix = []
    for i in range(0,len(mlist)):
        new_barcodes = []
        for b in mlist[i].barcodes:
            new_barcodes.append(h5list[i].split('/')[-1][:-3]+"@"+b.decode('UTF-8'))
        barcodes.append(numpy.array(new_barcodes, dtype='|S100'))
        matrix.append(mlist[i].matrix)
        
    write_10X_h5(filename, sp_sparse.hstack(matrix), features, numpy.concatenate(barcodes), genome, type)

def mtx_2_h5(filename, matrix_file, feature_file, barcode_file, genome = 'GRCh38', type = 'Peaks'):
    """Convert 10x mtx format matrix to HDF5."""
    matrix = scipy.io.mmread(matrix_file)
    matrix = sp_sparse.csc_matrix(matrix, dtype=numpy.int32)
    features = [row[0] for row in csv.reader(open(feature_file), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(open(barcode_file), delimiter="\t")]
    
    write_10X_h5(filename, matrix, features, barcodes, genome, type)

def count_2_h5(filename, count_file, genome = 'GRCh38', type = 'Peaks'):
    """Convert plain text count table to HDF5."""
    infile = open(count_file, 'r').readlines()
    barcodes = infile[0].strip().split('\t')
    features = []
    matrix = []
    for line in infile[1:]:
        line = line.strip().split('\t')
        features.append(line[0])
        matrix.append([float(t) for t in line[1:]])

    write_10X_h5(filename, matrix, features, barcodes, genome, type)



        

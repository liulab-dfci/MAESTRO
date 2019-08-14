import toolshed as ts
import os.path as op
import gzip
import sys
from itertools import imap

if len(sys.argv) != 5:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <states> <EDACC names> <input dir> <output dir>')
    sys.exit(1)

states_file_name = sys.argv[1]
edacc_names_file_name = sys.argv[2]
input_dir = sys.argv[3]
output_dir = sys.argv[4]

name_lookup = dict(x.strip().split() for x in ts.nopen(states_file_name))
file_lookup = dict(x.strip().split() for x in ts.nopen(edacc_names_file_name))

def process_file(fname):
    out_files = {}
    f = op.basename(fname).split("_")[0]
    fl = file_lookup[f]

    tmpl = output_dir + fl + "_%s.bed"

    for l in ts.nopen(fname):
        toks = l.rstrip().split("\t", 4)
        name = name_lookup[toks[3]]

        fname = tmpl % name
        if not fname in out_files:
            out_files[fname] = open(fname, 'w')
            print fname
        out_files[fname].write(l)


import glob
files = glob.glob(input_dir)

import multiprocessing
#multiprocessing.Pool(12).map(process_file, files)
for i in imap(process_file, files):
    pass

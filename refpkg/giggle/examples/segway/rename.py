import toolshed as ts
import os.path as op
import gzip
import sys
from itertools import imap

if len(sys.argv) != 3:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <input dir> <output dir>')
input_dir = sys.argv[1]
output_dir = sys.argv[2]

def process_file(fname):
    out_files = {}
    f = op.basename(fname).split(".")[0]
    tmpl = output_dir + '/' + f + "_%s.bed"

    for l in ts.nopen(fname):
        name = l.rstrip().split("\t", 4)[3].split('_')[1]

        fname = tmpl % name
        if not fname in out_files:
            out_files[fname] = open(fname, 'w')
            print fname
        out_files[fname].write(l)


import glob
files = glob.glob(input_dir)

import multiprocessing
multiprocessing.Pool(12).map(process_file, files)
#for i in imap(process_file, files):
    #pass

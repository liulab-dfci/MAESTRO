#!/usr/bin/python
import glob
import sys
import math
from optparse import OptionParser
import os

parser = OptionParser()

parser.add_option("--pretty_names",
                  dest="pretty_name_file",
                  help="Map of GEO id to pretty names for axis labels")

parser.add_option( "--q_x",
                  dest="q_x",
                  help="GEO ids to exclude from the query")


parser.add_option( "--db_x",
                  dest="db_x",
                  help="GEO ids to exclude from the DB")

parser.add_option("-d",
                  "--db_file",
                  dest="db_file",
                  help="File with database GEO sample ids to consider")

parser.add_option("-q",
                  "--query_file",
                  dest="query_file",
                  help="File with query GEO sample ids to consider")

parser.add_option("-c",
                  "--cistrome_data",
                  dest="cistrome_data_file",
                  help="Cistrome data")

parser.add_option("--qc",
                  dest="qc_data_file",
                  help="File with GEO id to file name map")

parser.add_option("--lc",
                  dest="line_count_file",
                  help="Number of lines in db/query file")

parser.add_option("--name_map",
                  dest="name_map_file",
                  help="File with cistrome name map")

parser.add_option("-i",
                  "--input_dir",
                  dest="input_dir",
                  help="Input directory")

parser.add_option("-o",
                  "--output",
                  dest="output_file",
                  help="Output file name")

parser.add_option("--x_size",
                  dest="x_size",
                  type="int",
                  default=10,
                  help="Figure x size (Default 10)")

parser.add_option("--y_size",
                  dest="y_size",
                  type="int",
                  default=30,
                  help="Figure x size (Default 30)")

(options, args) = parser.parse_args()

#if not options.output_file:
#    parser.error('Output file not given')
#if not options.input_dir:
#    parser.error('Iput directory not given')
#if not options.query_file:
#    parser.error('Querys id file not given')
#if not options.db_file:
#    parser.error('Database id file not given')

pretty_names = {}
if options.pretty_name_file:
    for l in open(options.pretty_name_file, 'r'):
        A = l.rstrip().split('\t')
        pretty_names[A[0]] = A[1:]

cistrome_to_pretty_name_map = {}
for l in open(options.name_map_file, 'r'):
    A = l.rstrip().split()
    cistrome_to_pretty_name_map[A[0]] = A[1]

pretty_name_to_geoid_map={}
cistrome_data_by_geoid={}
cistrome_id_to_geoid={}
for l in open(options.cistrome_data_file, 'r'):
    A = l.rstrip().split('\t')
    cistrome_data_by_geoid[A[1]] = {}
    cistrome_data_by_geoid[A[1]]['cis_id'] = A[0]
    cistrome_data_by_geoid[A[1]]['geo_id'] = A[1]
    cistrome_data_by_geoid[A[1]]['species'] = A[2]
    cistrome_data_by_geoid[A[1]]['cell_line'] = A[3]
    cistrome_data_by_geoid[A[1]]['tissue'] = A[4]
    cistrome_data_by_geoid[A[1]]['cell_type'] = A[5]
    cistrome_data_by_geoid[A[1]]['protein'] = A[6]
    cistrome_data_by_geoid[A[1]]['orig_file_name'] = A[7]

    pretty_name_to_geoid_map[cistrome_to_pretty_name_map[A[7]] + '.gz'] = A[1] 

    cistrome_id_to_geoid[A[0]] = A[1]


qc_by_geoid=None
if options.qc_data_file:
    qc_by_geoid={}
    for l in open(options.qc_data_file, 'r'):
        A = l.rstrip().split('\t')
        if A[0] == 'id':
            continue

        qc_by_geoid[cistrome_id_to_geoid[A[0]]] = { \
                'map' : A[1], \
                'peaks' : A[2], \
                'fastqc' : A[3], \
                'frip' : A[4], \
                'pbc' : A[5], \
                'motif_judge' : A[6], \
                'dhs' : A[7]}

line_count_by_geoid={}
for l in open(options.line_count_file, 'r'):
    A = l.rstrip().split('\t')
    line_count_by_geoid[pretty_name_to_geoid_map[A[0]]] = int(A[1])

DB_ids=[]
if options.db_file:
    for l in open(options.db_file, 'r'):
        DB_ids.append(l.rstrip().split()[0])

Q_ids=[]
if options.query_file:
    for l in open(options.query_file, 'r'):
        Q_ids.append(l.rstrip().split()[0])

D = []
X_names = []
Y_names = []

results = {}

Q_ids = [Q_id for Q_id in Q_ids if line_count_by_geoid[Q_id] > 20]
if qc_by_geoid:
    Q_ids = [Q_id for Q_id in Q_ids if qc_by_geoid[Q_id]['peaks'] == 'true']
    Q_ids = [Q_id for Q_id in Q_ids if qc_by_geoid[Q_id]['frip'] == 'true']

if options.q_x:
    Q_x = options.q_x.split(',')
    Q_ids = [Q_id for Q_id in Q_ids if Q_id not in Q_x]

DB_ids = [DB_id for DB_id in DB_ids if line_count_by_geoid[DB_id] > 20]
if qc_by_geoid:
    DB_ids = [DB_id for DB_id in DB_ids if qc_by_geoid[DB_id]['peaks'] == 'true']
    DB_ids = [DB_id for DB_id in DB_ids if qc_by_geoid[DB_id]['frip'] == 'true']


if options.db_x:
    DB_x = options.db_x.split(',')
    DB_ids = [DB_id for DB_id in DB_ids if DB_id not in DB_x]


for file_name in glob.glob(options.input_dir):
    query_file_name = '.'.join(file_name.split('/')[-1].split('.')[:-1])
    query_id = pretty_name_to_geoid_map[query_file_name]

    if query_id not in Q_ids:
        continue

    results[query_id] = {}

    for l in open(file_name, 'r'):
        A = l.rstrip().split('\t')
        if A[-1] == 'combo_score':
            continue

        db_id = pretty_name_to_geoid_map[A[0].split('/')[-1]]

        if  db_id not in DB_ids:
            continue

        results[query_id][db_id] = float(A[-1])


D = []

for Q_id in Q_ids:
    d = []
    for DB_id in DB_ids:
        if DB_id not in results[Q_id]:
            print Q_id, len(results[Q_id].keys()), 'GSM1024799' in results[Q_id]
            sys.exit(1)
        d.append(results[Q_id][DB_id])
    D.append(d)

Y_names = []
Y_colors = []
for Q_id in Q_ids:
    if Q_id in pretty_names:
        Y_names.append(pretty_names[Q_id][0])
        if len(pretty_names[Q_id]) > 1:
            Y_colors.append(pretty_names[Q_id][1])
        else:
            Y_colors.append('black')
    else:
        #Y_names.append(cistrome_data_by_geoid[Q_id]['protein'] + ' ' + \
                        #cistrome_data_by_geoid[Q_id]['geo_id'] )
        Y_names.append(cistrome_data_by_geoid[Q_id]['protein'])
        Y_colors.append('black')


X_names = []
X_colors = []
for DB_id in DB_ids:
    if DB_id in pretty_names:
        X_names.append(pretty_names[DB_id][0])
        if len(pretty_names[DB_id]) > 1:
            X_colors.append(pretty_names[DB_id][1])
        else:
            X_colors.append('black')
    else:
        #X_names.append(cistrome_data_by_geoid[DB_id]['protein'] + ' ' + \
                        #cistrome_data_by_geoid[DB_id]['geo_id'] )
        X_names.append(cistrome_data_by_geoid[DB_id]['protein'])
        X_colors.append('black')

#Y_names = Q_ids
#X_names = DB_ids

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

NP_D=np.zeros([len(Y_names),len(X_names)])

c = 0
for d in D:
    s = 0
    for v in d:
        NP_D[c,s] = v
        s+=1
    c+=1

print c,s

fig, ax = plt.subplots(figsize=(options.x_size,options.y_size), \
                               dpi=300)
from matplotlib import colors as mcolors
from matplotlib.colors import Normalize

_seismic_data = ( (0.0, 0.0, 0.3),
                  (0.0, 0.0, 1.0),

                  (1.0, 1.0, 1.0),

                  (1.0, 0.0, 0.0),
                  (0.5, 0.0, 0.0))

hm = mcolors.LinearSegmentedColormap.from_list( \
        name='red_white_blue', \
        colors=_seismic_data, N=256)

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

from pylab import get_cmap

data = NP_D
print data.min(), data.max()
plt.pcolor(data, 
           cmap=get_cmap("Reds"))
           #cmap=get_cmap("Rlues"))
           #cmap=hm)
           #norm = MidpointNormalize(midpoint=100),

ax.set_frame_on(False)
ax.xaxis.tick_top()
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.set_xticks(np.arange(len(X_names)) + 0.5, minor=False)
ax.set_yticks(np.arange(len(Y_names)) + 0.5, minor=False)
ax.set_xticklabels(X_names,rotation=90)


i = 0
for xtick in ax.get_xticklabels():
    xtick.set_color(X_colors[i])
    i+=1

i = 0
for ytick in ax.get_yticklabels():
    # colors from bottom to top
    ytick.set_color(Y_colors[len(Y_colors) - 1 -i])
    i+=1


if options.pretty_name_file:
    ax.set_yticklabels(Y_names[::-1])
else:
    ax.set_yticklabels(Y_names)
plt.ylim((0,len(Y_names)))
plt.xlim((0,len(X_names)))
cbar = plt.colorbar(fraction=0.046, pad=0.01, aspect=50)
cbar.ax.tick_params(labelsize=10)

plt.savefig(options.output_file,bbox_inches='tight')

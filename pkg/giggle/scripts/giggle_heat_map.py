#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-b",
                  action="store_true",
                  default=False,
                  dest="black",
                  help="black background")

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

parser.add_option("--no_labels",
                  dest="no_labels",
                  action="store_true",
                  default=False,
                  help="Do not label x and y axis");

parser.add_option("--no_ylabels",
                  dest="no_ylabels",
                  action="store_true",
                  default=False,
                  help="Do not label y axis");


parser.add_option("--stat",
                  dest="stat",
                  default="combo",
                  help="Stat to plot (odds, sig, combo) (Default: combo)")

parser.add_option("-i",
                  "--input",
                  dest="input_file",
                  help="GIGGLE result file name with '-s'")

parser.add_option("--ablines",
                  dest="ablines")

parser.add_option("-o",
                  "--output",
                  dest="output_file",
                  help="Output file name")

parser.add_option("-s",
                  "--states",
                  dest="states_file",
                  help="States file name")

parser.add_option("-n",
                  "--cell_names",
                  dest="cells_names_file",
                  help="Cells names file name")

parser.add_option("--sort",
                  dest="cells_sort_file",
                  help="Cells sort order file name")

parser.add_option("--state_names",
                  dest="state_names",
                  help="States names")

parser.add_option("--group_names",
                  dest="group_names",
                  help="Group names")

parser.add_option("-c",
                  "--cells",
                  dest="cells_file",
                  help="Cells file name")


(options, args) = parser.parse_args()

if not options.output_file:
    parser.error('Output file not given')
if not options.input_file:
    parser.error('Input file not given')
if not options.states_file:
    parser.error('States file not given')
if not options.cells_file:
    parser.error('Cells file not given')

if options.stat not in ['odds', 'sig', 'combo']:
    parser.error('Stat "' + options.stat + '" not supported')

states = []
for l in open(options.states_file, 'r'):
    states.append(l.rstrip().split('\t')[1])

cells = []
for l in open(options.cells_file, 'r'):
    cells.append(l.rstrip().split('\t')[1])

names = []
if options.cells_names_file:
    for l in open(options.cells_names_file, 'r'):
        names.append(l.rstrip().split('\t')[1])

sorts = []
if options.cells_sort_file:
    for l in open(options.cells_sort_file, 'r'):
        sorts.append(l.rstrip().split('\t')[1])

M = {}

for c in cells:
    M[c] = {}
    for s in states:
        M[c][s] = 0

for l in open(options.input_file,'r'):
    if l[0] == '#':
        continue
    A = l.rstrip().split('\t')
    odds = float(A[3])
    sig = float(A[4])
    combo = float(A[7])

    c = ''
    s = ''
    found = False
    for i in cells:
        if '/' + i in A[0]:
            assert not found, (A[0], i)
            c = i
            found = True

    found = False
    for i in states:
        if c + '_' + i in A[0]:
            assert not found, (A[0], i)
            s = i
            found = True

    assert M[c][s] == 0, (A[0], c, s, M[c][s])

    if options.stat == 'sig':
        M[c][s] = sig
    elif options.stat == 'combo':
        M[c][s] = combo
    elif options.stat == 'odds':
        M[c][s] = odds

D=np.zeros([len(cells),len(states)])

c = 0
for cell in cells:
    s = 0
    for state in states:
        D[c,s] = M[cell][state]
        s+=1
    c+=1

column_labels = cells
row_labels = states
data = D


fig, ax = plt.subplots(figsize=(options.x_size,options.y_size), \
                       dpi=300)
if options.black:
    fig =  matplotlib.pyplot.figure(figsize=(options.x_size,options.y_size), \
                                    dpi=300, \
                                    facecolor='black')

    ax = fig.add_subplot(1,1,1,axisbg='k')

#plt.pcolor(data, cmap=plt.cm.Blues)
#plt.pcolor(data, cmap=plt.cm.bwr)

if options.cells_sort_file:
    idxs = np.argsort(sorts)
    sorts=np.array(sorts)[idxs]
    data=data[idxs]
elif options.cells_names_file:
    idxs = np.argsort(names)
    names=np.array(names)[idxs]
    data=data[idxs]

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


#norm=matplotlib.colors.LogNorm(vmin=data.min(), vmax=data.max()), \
print data.min(), data.max()

plt.pcolor(data, \
           norm = MidpointNormalize(midpoint=0),
           cmap=hm)
           #cmap=plt.cm.seismic)
ax.xaxis.tick_top()

ax.xaxis.set_ticks_position('none') 
ax.yaxis.set_ticks_position('none') 

if options.ablines:
    for l in [int(x) for x in options.ablines.split(',')]:
        ax.axhline(l, linestyle='--', color='k') 

if options.no_labels:
    ax.set_xticklabels("")
    ax.set_yticklabels("")
elif options.no_ylabels:
    ax.set_yticklabels("")

    if options.state_names:
        state_names = []
        for l in open(options.state_names, 'r'):
            state_names.append(l.rstrip())
        if options.black:
            plt.xticks(np.arange(0.5,len(state_names)+0.5,1.0),state_names,rotation=90, fontsize=10, color='white')
        else:
            plt.xticks(np.arange(0.5,len(state_names)+0.5,1.0),state_names,rotation=90, fontsize=10)
    else:
        if options.black:
            plt.xticks(np.arange(0.5,len(states)+0.5,1.0),states,rotation=90, fontsize=10, color='white')
        else:
            plt.xticks(np.arange(0.5,len(states)+0.5,1.0),states,rotation=90, fontsize=10)
else:
    if options.group_names:

        group_names = []
        for l in open(options.group_names, 'r'):
            group_names.append(l.rstrip())
        group_names.reverse()
        if options.black:
            plt.yticks(np.arange(0.5,len(group_names)+0.5,1.0),[x.decode("utf8") for x in group_names], fontsize=10, color='white')
        else:
            plt.yticks(np.arange(0.5,len(group_names)+0.5,1.0),[x.decode("utf8") for x in group_names], fontsize=10)

    elif options.cells_names_file:
        if options.black:
            plt.yticks(np.arange(0.5,len(names)+0.5,1.0),[x.decode("utf8") for x in names], fontsize=10, color='white')
        else:
            plt.yticks(np.arange(0.5,len(names)+0.5,1.0),[x.decode("utf8") for x in names], fontsize=10)

    else: 
        if options.black:
            plt.yticks(np.arange(0.5,len(cells)+0.5,1.0),cells, fontsize=10, color='white')
        else:
            plt.yticks(np.arange(0.5,len(cells)+0.5,1.0),cells, fontsize=10)

    if options.state_names:
        state_names = []
        for l in open(options.state_names, 'r'):
            state_names.append(l.rstrip())
        if options.black:
            plt.xticks(np.arange(0.5,len(state_names)+0.5,1.0),state_names,rotation=90, fontsize=10, color='white')
        else:
            plt.xticks(np.arange(0.5,len(state_names)+0.5,1.0),state_names,rotation=90, fontsize=10)
    else:
        if options.black:
            plt.xticks(np.arange(0.5,len(state_names)+0.5,1.0),state_names,rotation=90, fontsize=10, color='white')
        else:
            plt.xticks(np.arange(0.5,len(states)+0.5,1.0),states,rotation=90, fontsize=10)

plt.ylim((0,len(cells)))
plt.xlim((0,len(states)))


#cbar_min = -1 * max([abs(data.min()), data.max()])
#cbar_max = max([abs(data.min()), data.max()])

cbar = plt.colorbar(fraction=0.046, pad=0.04, ticks=[data.min(), 0, data.max()])
cbar.ax.tick_params(labelsize=10) 

if options.black:
    cbar.ax.tick_params(color='white', labelcolor='white')
    plt.savefig(options.output_file,bbox_inches='tight',\
                facecolor=fig.get_facecolor(),\
                transparent=True)
else:
    plt.savefig(options.output_file,bbox_inches='tight')

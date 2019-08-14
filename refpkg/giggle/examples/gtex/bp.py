#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab
import random
from optparse import OptionParser
import seaborn

delim = '\t'
parser = OptionParser()

parser.add_option("-l",
                  "--log_y",
                  action="store_true", dest="logy", default=False,
                  help="Use log scale for y-axis")

parser.add_option("-o",
                  "--output_file",
                  dest="output_file",
                  help="Data file")

parser.add_option("--y_min",
                  dest="min_y",
                  help="Min y value")

parser.add_option("--y_max",
                  dest="max_y",
                  help="Max y value")

parser.add_option("--line_style",
                  dest="line_style",
                  default=".",
                  help="Line style")



(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

Y=[]
row_names = []
col_names = ''
for l in sys.stdin:
    if col_names == '':
        col_names = l.rstrip().split()
        continue
    Y.append([])
    A = l.rstrip().split()
    row_names.append(A[0])
    for y in A[1:]:
        Y[-1].append([float(x) for x in y.split(',')])

matplotlib.rcParams.update({'font.size': 12})
fig = matplotlib.pyplot.figure(figsize=(4,14),dpi=300)
fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)
N = len(Y)
i = 1
for y in Y:
    ax = fig.add_subplot(N,1,i)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    #ax.plot(range(len(y)),y,options.line_style,color='black', linewidth=1)
    ax.boxplot(y)
    ticks = ax.set_yticklabels([])
    ax.text(0.5,max([max(x) for x in y]), row_names[i-1], fontsize=8, va='top', ha='left')

    if i != 1:
        ax.set_xticklabels([])
    else:
        ax.xaxis.tick_top()
        ax.set_xticklabels(col_names)
        matplotlib.pyplot.xticks(rotation=90, fontsize=8)

    i+=1

matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight')

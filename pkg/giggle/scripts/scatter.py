#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab
import random
from optparse import OptionParser

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'

delim = '\t'
parser = OptionParser()

parser.add_option("-a",
                  "--alpha",
                  type="float",
                  default=1.0,
                  help="Alpha value for points")


parser.add_option("-l",
                  "--log_y",
                  action="store_true", dest="logy", default=False,
                  help="Use log scale for y-axis")

parser.add_option("-o",
                  "--output_file",
                  dest="output_file",
                  help="Data file")

parser.add_option("--title",
                  dest="title",
                  help="Plot title")


parser.add_option("--x_label",
                  dest="x_label",
                  help="X axis label")

parser.add_option("--y_label",
                  dest="y_label",
                  help="Y axis label")



parser.add_option("--y_min",
                  dest="min_y",
                  help="Min y value")

parser.add_option("--y_max",
                  dest="max_y",
                  help="Max y value")

parser.add_option("--x_min",
                  dest="min_x",
                  help="Min x value")

parser.add_option("--x_max",
                  dest="max_x",
                  help="Max x value")


parser.add_option("--line_style",
                  dest="line_style",
                  default=".",
                  help="Line style")

parser.add_option("--fig_x",
                  dest="fig_x",
                  type="int",
                  default=5,
                  help="Figure width")

parser.add_option("--fig_y",
                  dest="fig_y",
                  type="int",
                  default=5,
                  help="Figure height")

parser.add_option("-b",
                  action="store_true", 
                  default=False,
                  dest="black",
                  help="black background")

parser.add_option("-c",
                  "--color",
                  dest="color",
                  default="black",
                  help="Color")

parser.add_option("--trend",
                  action="store_true", 
                  default=False,
                  dest="trend",
                  help="Trend line")

parser.add_option("--point_size",
                  dest="point_size",
                  type="int",
                  default=1,
                  help="Scatter point size (defualt 1)")



(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

X=[]
Y=[]
E=[]
for l in sys.stdin:
    #a = l.rstrip().split(delim)
    a = l.rstrip().split()
    if len(a) == 1:
        Y.append(float(a[0]))
    if len(a) == 2:
        X.append(float(a[0]))
        Y.append(float(a[1]))
    if len(a) == 3:
        X.append(float(a[0]))
        Y.append(float(a[1]))
        E.append(float(a[2]))


matplotlib.rcParams.update({'font.size': 12})
#fig = matplotlib.pyplot.figure(figsize=(options.fig_x,options.fig_y),dpi=300)
fig = 1
if options.black:
    fig = matplotlib.pyplot.figure(\
            figsize=(options.fig_x,options.fig_y),\
            dpi=300,\
            facecolor='black')
else:
    fig = matplotlib.pyplot.figure(\
            figsize=(options.fig_x,options.fig_y),\
            dpi=300)

fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)

#ax = fig.add_subplot(1,1,1)
ax = 1
if options.black:
    ax = fig.add_subplot(1,1,1,axisbg='k')
else:
    ax = fig.add_subplot(1,1,1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

if options.black:
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.title.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.xaxis.label.set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

if len(X) == 0:
    ax.plot(range(len(Y)),\
            Y,\
            options.line_style,color=options.color, \
            s=1,\
            linewidth=1,
            alpha=options.alpha)
else:
    ax.plot(X,Y,options.line_style,color=options.color, linewidth=1, \
            alpha=options.alpha)
#ax.scatter(X,Y,s=options.point_size,color=options.color)

#print X
#print Y
#ax.errorbar(X, Y, yerr=E, fmt='-o', color=options.color)
if len(E)!=0:
    ax.errorbar(X,Y,yerr=E, linestyle="None", color='gray')


if options.logy:
    ax.set_yscale('log')

if ((options.max_y) and (options.min_y)):
    ax.set_ylim(float(options.min_y),float(options.max_y))

if ((options.max_x) and (options.min_x)):
    ax.set_xlim(float(options.min_x),float(options.max_x))

#if len(X) != 0:
#    ax.set_xticks([float(x) for x in X])
#    #ax.set_xticklabels

if options.x_label:
    ax.set_xlabel(options.x_label)

if options.y_label:
    ax.set_ylabel(options.y_label)

if options.title:
    ax.set_title(options.title)

if options.trend:
    z = np.polyfit(X, Y, 1)
    p = np.poly1d(z)
    ax.plot(X,p(Y),'r--',color='red')





#matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight')

if options.black:
    matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight',\
            facecolor=fig.get_facecolor(),\
              transparent=True)
else:
    matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight')

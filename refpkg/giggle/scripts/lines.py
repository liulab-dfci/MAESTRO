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
rcParams['legend.numpoints'] = 1


delim = '\t'
parser = OptionParser()

parser.add_option("--xticks",
                  dest="xticks",
                  help="CSV of x tick lables")


parser.add_option("--numyticks",
                  dest="numyticks",
                  help="Number of Y ticks")


parser.add_option("-c",
                  "--colors",
                  dest="colors",
                  default="blue,green,red,magenta,black",
                  help="Color CSV (blue,green,red,magenta,black)")


parser.add_option("-o",
                  "--output_file",
                  dest="output_file",
                  help="Data file")

parser.add_option("-l",
                  "--line_stype",
                  default='-o',
                  dest="line_style",
                  help="Line style (Default '-o')")

parser.add_option("--plot_width",
                  dest="plot_width",
                  type='int',
                  default=2,
                  help="Line width (Default 1)")

parser.add_option("--plot_height",
                  dest="plot_height",
                  type='int',
                  default=3,
                  help="Line width (Default 1)")


parser.add_option("-w",
                  "--line_width",
                  dest="line_width",
                  type='int',
                  default=1,
                  help="Line width (Default 1)")

parser.add_option("-a",
                  "--ablines",
                  dest="ablines",
                  help="Ablines (e.g. h:10:-o:white)")

parser.add_option("--xbins",
                  dest="xbins",
                  help="Number of ticks on x-axis")



parser.add_option("--legend",
                  dest="legend",
                  help="Comma sperated legend")

parser.add_option("--legend_loc",
                  dest="legend_loc",
                  type="int",
                  default=4,
                  help="Legend location")


parser.add_option("--title",
                  dest="title",
                  help="Title")

parser.add_option("--xlabel",
                  dest="xlabel",
                  help="X axis label")

parser.add_option("--ylabel",
                  dest="ylabel",
                  help="Y axis label")

parser.add_option( "--xlog",
                  action="store_true", 
                  default=False,
                  dest="xlog",
                  help="X axis log")


parser.add_option( "--ylog",
                  action="store_true", 
                  default=False,
                  dest="ylog",
                  help="Y axis log")

parser.add_option("-X",
                  action="store_true", 
                  default=False,
                  dest="X",
                  help="X values includeded (Y line i, X line i+1)")

parser.add_option("-b",
                  action="store_true", 
                  default=False,
                  dest="black",
                  help="black background")

parser.add_option("--x_max",
                  dest="max_x",
                  type="float",
                  help="Max x value")

parser.add_option("--x_min",
                  dest="min_x",
                  type="float",
                  help="Min x value")

parser.add_option("--y_max",
                  dest="max_y",
                  type="float",
                  help="Max y value")

parser.add_option("--y_min",
                  dest="min_y",
                  type="float",
                  help="Min y value")

(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

matplotlib.rcParams.update({'font.size': 12})
fig = 1
if options.black:
    fig = matplotlib.pyplot.figure(\
            figsize=(options.plot_width,options.plot_height),\
            dpi=300,\
            facecolor='black')
else:
    fig = matplotlib.pyplot.figure(\
            figsize=(options.plot_width,options.plot_height),\
            dpi=300)
if options.xbins:
    matplotlib.pyplot.locator_params(axis = 'x', nbins = options.xbins)
#fig.subplots_adjust(top=1.0)

ax = 1
if options.black:
    ax = fig.add_subplot(1,1,1,axisbg='k')
else:
    ax = fig.add_subplot(1,1,1)

ax.set_axisbelow(True)
ax.tick_params(labelsize=10)


#colors = [ 'blue', \
#    'green', \
#    'red', \
#    'magenta', \
#    'black']
colors=options.colors.split(',')

color_i = 0
plts=[]
lines = sys.stdin.readlines()

if (options.X):
    for i in range(len(lines))[::2]:
        Y = [float(x) for x in lines[i].rstrip().split()]
        X = [float(x) for x in lines[i+1].rstrip().split()]
        p, = ax.plot(X,\
                     Y,\
                     options.line_style,\
                     color=colors[color_i],\
                     linewidth=options.line_width,
                     markeredgewidth=0)
        plts.append(p)
        color_i = (color_i + 1) % len(colors)
else:
    for i in range(len(lines))[::1]:
        Y = [float(x) for x in lines[i].rstrip().split()]
        p, = ax.plot(range(len(Y)), \
                     Y,\
                     options.line_style,\
                     color=colors[color_i], \
                     linewidth=options.line_width,
                     markeredgewidth=0)
        plts.append(p)
        color_i = (color_i + 1) % len(colors)

if options.ylog:
    ax.set_yscale('log')

if options.xlog:
    ax.set_xscale('log')

if options.legend:
    l1=ax.legend(plts, options.legend.split(","), frameon=False, fontsize=10,loc=options.legend_loc)
    if options.black:
        for text in l1.get_texts():
            matplotlib.pyplot.setp(text,color='white')


if options.title:
    #matplotlib.pyplot.suptitle(options.title)
    ax.set_title(options.title)

if options.xlabel:
    ax.set_xlabel(options.xlabel)

if options.ylabel:
    ax.set_ylabel(options.ylabel)

if options.max_x:
    ax.set_xlim(xmax=options.max_x)
if options.min_x:
    ax.set_xlim(xmin=options.min_x)

if options.max_y:
    ax.set_ylim(ymax=options.max_y)
if options.min_y:
    ax.set_ylim(ymin=options.min_y)

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
#ax.legend()


if options.ablines:
    #h:10:-o:white
    for ab in options.ablines.split(","):
        ab_d,ab_pos,ab_ls,ab_col = ab.split(':')
        ab_pos=float(ab_pos)
        if ab_d == 'v':
            ab_min, ab_max = ax.get_ylim()
            ax.plot([ab_pos,ab_pos],[ab_min,ab_max],ab_ls,c=ab_col)
        else:
            ab_min, ab_max = ax.get_xlim()
            ax.plot([ab_min,ab_max],[ab_pos,ab_pos],ab_ls,c=ab_col)


if options.numyticks:
    matplotlib.pyplot.locator_params(axis = 'y', nbins = int(options.numyticks))

if options.xticks:
    ax.set_xticklabels(options.xticks.split(','))

if options.black:
    matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight',\
            facecolor=fig.get_facecolor(),\
              transparent=True)
else:
    matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight')

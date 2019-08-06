#!/usr/bin/env python
import sys
import math
from operator import itemgetter
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab
import random
from optparse import OptionParser
from matplotlib import colors as mcolors

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'

delim = '\t'


def log2fc(v):
    if v == 0.0 :
        return 0.0

    if v < 1:
        v = 1.0/v
        return -math.log(v,2)

    return math.log(v,2)

def neglog10p(v):
    if v == 0.0 :
        return 6.0
    return -math.log(v,10)


parser = OptionParser()

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

D = []
for l in sys.stdin:
    if l[0] == '#':
        continue
    A = l.rstrip().split()
    ratio = float(A[3])
    pval = float(A[4])
    D.append([abs(log2fc(ratio)* neglog10p(pval)), log2fc(ratio), neglog10p(pval)])

D.sort(key=itemgetter(0), reverse=False)
D = [x + [i] for i, x in enumerate(D)]
X=[]
Y=[]
R=[]
E=[]

for d in D:
    #R.append(d[3]**2)
    R.append(1+d[3])
    X.append(d[1])
    Y.append(d[2])


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


_seismic_data = ( (0.0, 0.0, 0.3),
                  (0.0, 0.0, 1.0),
                  (1.0, 1.0, 1.0),
                  (1.0, 0.0, 0.0),
                  (0.5, 0.0, 0.0))

def mapc(rgb_l):
    return [x/256.0 for x in rgb_l] 

aliceblue = (240,248,255)
lavender = (230,230,250)
powderblue = (176,224,230)
lightblue = (173,216,230)
lightskyblue = (135,206,250)
skyblue = (135,206,235)
deepskyblue = (0,191,255)
lightsteelblue = (176,196,222)
dodgerblue = (30,144,255)
cornflowerblue = (100,149,237)
steelblue = (70,130,180)
cadetblue = (95,158,160)
mediumslateblue = (123,104,238)
slateblue = (106,90,205)
darkslateblue = (72,61,139)
royalblue = (65,105,225)
blue = (0,0,255)
mediumblue = (0,0,205)
darkblue = (0,0,139)
navy = (0,0,128)
midnightblue = (25,25,112)

_seismic_data = (
                 mapc((201,201,243)),
                 mapc((201,201,243)),
                 mapc((201,201,243)),
                 mapc((201,201,243)),
                 mapc((201,201,243)),
                 mapc((201,201,243)),
                 mapc((139,139,229)),
                 mapc((139,139,229)),
                 mapc((139,139,229)),
                 mapc((139,139,229)),
                 mapc((139,139,229)),
                 mapc((11,11,49)))


hm = mcolors.LinearSegmentedColormap.from_list( \
    name='red_white_blue', \
    colors=_seismic_data, N=256)

sc = None
if len(X) == 0:
    sc = ax.plot(range(len(Y)),Y,options.line_style,color=options.color, s=1,linewidth=1)
else:
    #sc = ax.scatter(X,Y,c=R,lw = 0,cmap=matplotlib.pyplot.cm.Blues)
    #sc = ax.scatter(X,Y,c=R,lw = 0,cmap=matplotlib.pyplot.cm.seismic)
    sc = ax.scatter(X,Y,c=R,lw = 0,cmap=hm)

if len(E)!=0:
    ax.errorbar(X,Y,yerr=E, linestyle="None", color='gray')


if options.logy:
    ax.set_yscale('log')

if ((options.max_y) and (options.min_y)):
    ax.set_ylim(float(options.min_y),float(options.max_y))

if ((options.max_x) and (options.min_x)):
    ax.set_xlim(float(options.min_x),float(options.max_x))

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

print len(R)
#cbar = matplotlib.pyplot.colorbar(fraction=0.046, pad=0.04)
#labels=[1,500,1000,1500,1905]
labels=[1,1000,1905]
#cbar = matplotlib.pyplot.colorbar(sc, ticks=labels,fraction=0.015, pad=0.04)
cbar = matplotlib.pyplot.colorbar(sc, ticks=labels, shrink=.25, aspect=7, pad=0.0)
cbar.ax.set_yticklabels([1905,1000,1], size=10)  # horizontal colorbar
#cbar.ax.set_title('Rank',y=1.08)

#cbar.ax.set_xticklabels(list(labels).reverse())
#cbar.ax.set_xticklabels(labels)

if options.black:
    matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight',\
            facecolor=fig.get_facecolor(),\
              transparent=True)
else:
    matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight')

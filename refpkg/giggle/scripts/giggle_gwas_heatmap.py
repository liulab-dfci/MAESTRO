#!/usr/bin/python
import glob
import sys
import math
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np
import os

parser = OptionParser()

parser.add_option("-o",
                  "--output",
                  dest="output_file",
                  help="Output file name")

parser.add_option("-i",
                  "--input_dir",
                  dest="input_dir",
                  help="Input directory")

parser.add_option("-s",
                  "--states_file",
                  dest="states_file",
                  help="Genomic states")

parser.add_option("--states",
                  dest="states_list",
                  help="CSV or states to consider")

parser.add_option("-n",
                  "--cell_names",
                  dest="cells_names_file",
                  help="Cells names file name")

parser.add_option("-c",
                  "--cells",
                  dest="cells_file",
                  help="Cells file name")

parser.add_option("--group_names",
                  dest="group_names",
                  help="Group names")

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

parser.add_option("-b",
                  action="store_true",
                  default=False,
                  dest="black",
                  help="black background")

parser.add_option("--ablines",
                  dest="ablines")

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

(options, args) = parser.parse_args()
if not options.input_dir:
    parser.error('Input directory not given')
if not options.states_file:
    parser.error('States file not given')
if not options.states_list:
    parser.error('States file not given')
if not options.cells_file:
    parser.error('Cells file not given')
if not options.cells_names_file:
    parser.error('Cells names file not given')

viz_states = options.states_list.split(',')

states = []
for l in open(options.states_file, 'r'):
    states.append(l.rstrip().split('\t')[1])

sorted_states = states[::]
sorted_states.sort(lambda x,y: cmp(len(y), len(x)))

cells = []
for l in open(options.cells_file, 'r'):
    cells.append(l.rstrip().split('\t')[1])

names = []
for l in open(options.cells_names_file, 'r'):
    names.append(l.rstrip().split('\t')[1])


#
#group_names = {}
#for l in open(options.group_names_file, 'r'):
#    A = l.rstrip().split('\t')
#    group_names[A[0]] = A[1]
#
#cell_names = {}
#cells_to_groups = {}
#groups_to_cells = {}
#names = []
#for l in open(options.cells_file, 'r'):
#    #cells.append(l.rstrip().split('\t')[1])
#    A = l.rstrip().split('\t')
#    names.append(A[1])
#    cell_names[A[0]] = A[1]
#    if A[0] in group_names:
#        cells_to_groups[A[1]] = group_names[A[0]]
#        groups_to_cells[group_names[A[0]]] = A[1]
#
#all_states = []
#for l in open(options.states_file_name, 'r'):
#    A = l.rstrip().split('\t')
#    all_states.append(A[1])
#
#all_states.sort(lambda x,y: cmp(len(y), len(x)))
#
#for state in states:
#    if state not in all_states:
#        print 'Unknown state: ' + state
#        sys.exit(1)
#
results = {}

for file_name in glob.glob(options.input_dir):
    trait = file_name.split('/')[-1].split('.')[0]

    results[trait] = {}

    for l in open(file_name, 'r'):
        A = l.rstrip().split('\t')
        if A[-1] == 'combo_score':
            continue

        curr_state = ''
        for state in sorted_states:
            if state in A[0]:
                curr_state = state
                break

        tissue = A[0].split('/')[-1][:-7][:-1*(len(curr_state) + 1)]

        if tissue not in results[trait]:
            results[trait][tissue] = {}

        if curr_state in viz_states:
            results[trait][tissue][curr_state] = float(A[7])

traits = results.keys()
D=np.zeros([len(cells),2*len(traits)])

print len(traits)

print '\n'.join(traits)

sorted_traits = ['Celiac_disease',
                 'Autoimmune_thyroiditis',
                 'Primary_biliary_cirrhosis',
                 'Rheumatoid_arthritis',
                 'Asthma',
                 'Allergy',
                 'Kawasaki_disease',
                 'Behcets_disease',
                 'Vitiligo',
                 'Alopecia_areata',
                 'Systemic_lupus_erythematosus',
                 'Systemic_sclerosis',
                 'Type_1_diabetes',
                 'Crohns_disease',
                 'Ulcerative_colitis',
                 'Ankylosing_spondylitis',
                 'Atopic_dermatitis',
                 'Primary_sclerosing_cholangitis',
                 'Juvenile_idiopathic_arthritis',
                 'Psoriasis',
                 'Multiple_sclerosis',
                 'Renal_function_related_traits_BUN',
                 'Fasting_glucose_related_traits',
                 'Migraine',
                 'HDL_cholesterol',
                 'Platelet_counts',
                 'LDL_cholesterol',
                 'Restless_legs_syndrome',
                 'Bone_mineral_density',
                 'Red_blood_cell_traits',
                 'Chronic_kidney_disease',
                 'Progressive_supranuclear_palsy',
                 'Alzheimers_combined',
                 'Liver_enzyme_levels_gamma_glutamyl_transferase',
                 'Urate_levels',
                 'Type_2_diabetes',
                 'Creatinine_levels',
                 'Triglycerides',
                 'C_reactive_protein' ]

c = 0
for cell in cells:
    t = 0
    for trait in sorted_traits:
    #for trait in traits:
        #D[c,t] = M[cell][state]
        D[c,t] = results[trait][cell][viz_states[0]]
        t+=1
    for trait in sorted_traits:
        D[c,t] = results[trait][cell][viz_states[1]]
        t+=1
    c+=1

data = D
idxs = np.argsort(names)
names=np.array(names)[idxs]
data=data[idxs]


fig, ax = plt.subplots(figsize=(options.x_size,options.y_size), \
                       dpi=300)
if options.black:
    fig =  matplotlib.pyplot.figure(figsize=(options.x_size,options.y_size), \
                                    dpi=300, \
                                    facecolor='black')

    ax = fig.add_subplot(1,1,1,axisbg='k')


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

ax.axvline(21, linestyle='--', color='k') 
ax.axvline(39, linestyle='-', color='k') 
ax.axvline(39+21, linestyle='--', color='k') 


ax.set_xticklabels("")
ax.set_yticklabels("")

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

#    if options.state_names:
#        state_names = []
#        for l in open(options.state_names, 'r'):
#            state_names.append(l.rstrip())
#        if options.black:
#            plt.xticks(np.arange(0.5,len(state_names)+0.5,1.0),state_names,rotation=90, fontsize=10, color='white')
#        else:
#            plt.xticks(np.arange(0.5,len(state_names)+0.5,1.0),state_names,rotation=90, fontsize=10)
#    else:
#        if options.black:
#            plt.xticks(np.arange(0.5,len(state_names)+0.5,1.0),state_names,rotation=90, fontsize=10, color='white')
#        else:
#            plt.xticks(np.arange(0.5,len(states)+0.5,1.0),states,rotation=90, fontsize=10)

plt.ylim((0,len(cells)))
plt.xlim((0,2*len(traits)))
cbar = plt.colorbar(fraction=0.026, pad=0.02, ticks=[data.min(), 0, data.max()])
cbar.ax.tick_params(labelsize=10) 

if options.black:
    cbar.ax.tick_params(color='white', labelcolor='white')
    plt.savefig(options.output_file,bbox_inches='tight',\
                facecolor=fig.get_facecolor(),\
                transparent=True)
else:
    plt.savefig(options.output_file,bbox_inches='tight')

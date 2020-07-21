# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-07-17 18:54:18
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-07-17 20:14:07


import argparse
import pysam
import time, os
from collections import defaultdict
from contextlib import ExitStack


chr_list = ['chr'+str(i) for i in list(range(1,23))] 
chr_list = chr_list + ['chrX', 'chrY']

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "Split the fragement file according to the cluster.")
    parser.add_argument("-F","--frag", dest = "fragfile", default = "",help = "The fragment file with barcodes and counts.")
    parser.add_argument("-C","--cluster", dest = "clusterfile", default = "",help = "The barcode-cluster file generated in the step of analysis.")
    parser.add_argument("-O", "--outdir", dest = "outdir", default = "", help = "The output directory.")

    return parser.parse_args()


parser = CommandLineParser()
fragfile = parser.fragfile
clusterfile = parser.clusterfile
outdir = parser.outdir

try:
    os.makedirs(outdir)
except OSError:
    pass

start_time = time.time()
print("Start to read the cell-cluter file file.",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
cell_cluster_fict = {}
with open(clusterfile, "r") as cluster_in:
    for line in cluster_in:
        items = line.strip().split("\t")
        cell_cluster_fict[items[0]] = items[1]
end_time = time.time()
print("End:", end_time-start_time)

cluster_list = list(set(list(cell_cluster_fict.values())))
cluster_index_dict = {}
for i in range(len(cluster_list)):
    cluster_index_dict[cluster_list[i]] = i

start_time = time.time()
print("Start to split the fragement file.",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
with ExitStack() as stack:
    out_files = [stack.enter_context(open(os.path.join(outdir, cluster + ".bed"), "w")) for cluster in cluster_list]
    with open(fragfile, "r") as frag_in:
        for line in frag_in.readlines():
            items = line.strip().split("\t")
            if items[3] in cell_cluster_fict:
                cluster = cell_cluster_fict[items[3]]
                index = cluster_index_dict[cluster]
                if int(items[4]) == 1:
                    out_files[index].write("\t".join(items[0:3]) + "\n")
                else:
                    for i in range(int(items[4])):
                        out_files[index].write("\t".join(items[0:3]) + "\n")
end_time = time.time()
print("End:", end_time-start_time)


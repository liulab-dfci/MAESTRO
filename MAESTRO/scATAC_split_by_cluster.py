#! /usr/bin/env python

import argparse
import pysam
import time, os
import sys
from collections import defaultdict
from contextlib import ExitStack
from MAESTRO.scATAC_utility import *


def CommandLineParser():
    parser=argparse.ArgumentParser(description = "Split the fragement file according to the cluster and sample.")
    parser.add_argument("-F","--frag", dest = "fragfile", default = "",help = "The fragment file with barcodes and counts.")
    parser.add_argument("-C","--cluster", dest = "clusterfile", default = "",help = "The barcode-cluster file generated in the step of analysis.")
    parser.add_argument("-S","--split", dest = "split", default = "by_cluster", choices = ["by_cluster", "by_sample_cluster"], type = str, help = "Users defined split option. DEFAULT: by_cluster.")
    parser.add_argument("-O", "--outdir", dest = "outdir", default = "output/", help = "The output directory.")
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()



parser = CommandLineParser()
fragfile = parser.fragfile
clusterfile = parser.clusterfile
split = parser.split
outdir = parser.outdir

try:
    os.makedirs(outdir)
except OSError:
    pass

if split == "by_cluster":
    cell_cluster_fict = {}

    start_time = time.time()
    print("Start to read the barcode-cluster file.",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    with universal_open(clusterfile, "r") as cluster_in:
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
    print("Start to read the fragment file.",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    with ExitStack() as stack:
        out_files = [stack.enter_context(open(os.path.join(outdir, cluster + ".tsv"), "w")) for cluster in cluster_list]
        with universal_open(fragfile, "rt") as frag_in:
            for line in frag_in:
                items = line.strip().split("\t")
                if items[3] in cell_cluster_fict:
                    cluster = cell_cluster_fict[items[3]]
                    index = cluster_index_dict[cluster]
                    out_files[index].write("\t".join(items[0:4]) + "\n")
    end_time = time.time()
    print("End:", end_time-start_time)

elif split == "by_sample_cluster":
    sample_cluster_fict = {}

    start_time = time.time()
    print("Start to read the barcode-cluster file.",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    with universal_open(clusterfile, "rt") as cluster_in:
        for line in cluster_in:
            items = line.strip().split("\t")
            samples = items[0].split("@")
            sample_cluster = items[1] + "@" + samples[0]
            #cell_cluster_fict[items[0]] = items[1]
            sample_cluster_fict[items[0]] = sample_cluster
    end_time = time.time()
    print("End:", end_time-start_time)
    #cluster_list = list(set(list(cell_cluster_fict.values())))
    #cluster_index_dict = {}
    #for i in range(len(cluster_list)):
        #cluster_index_dict[cluster_list[i]] = i


    sample_cluster_list = list(set(list(sample_cluster_fict.values())))
    sample_cluster_list.sort()
    sample_cluster_index_dict = {}
    for i in range(len(sample_cluster_list)):
        sample_cluster_index_dict[sample_cluster_list[i]] = i

    start_time = time.time()
    print("Start to read the fragment file.",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    with ExitStack() as stack:
        out_files = [stack.enter_context(open(os.path.join(outdir, sample_cluster + ".tsv"), "w")) for sample_cluster in sample_cluster_list]
        with universal_open(fragfile, "rt") as frag_in:
            for line in frag_in:
                items = line.strip().split("\t")
                if items[3] in sample_cluster_fict:
                    sample_cluster = sample_cluster_fict[items[3]]
                    index = sample_cluster_index_dict[sample_cluster]
                    out_files[index].write("\t".join(items[0:4]) + "\n")

    end_time = time.time()
    print("End:", end_time-start_time)

else:
    print("split method not defined correctly.")

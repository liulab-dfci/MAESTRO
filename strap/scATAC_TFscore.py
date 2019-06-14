#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:32:00 2019

@author: Chenfei Wang
"""

import sys
import re
import json
import random
import math
import scipy.stats as stats
import multiprocessing as mp
import time

def filter_tf_peaks(peak_file, gc_file, tfindex_file, tfpercent):
    """ Filter the TF peaks based on the overlap with single-cell peaks. """
    
    gc_anno = {}                                                   # store gc annotations
    for line in open(gc_file, 'r'):
        line = line.strip().split('\t')
        gc_anno['_'.join(line[:3])] = [line[3],int(float(line[4]))]
    
    all_peak = {}                                                 # store gc percentage for allpeaks
    sc_peak = {}                                                  # store peak index for each single-cell
    cells = []
    for line in open(peak_file, 'r'):
        line = line.strip().split('\t')
        if not line[0][:3] == 'chr':
           cells = line
           for cell in cells:
               sc_peak[cell] = []
        else:
           pos = line[0].split('_')
           id = pos[0]+'_'+str(int((int(pos[1])+int(pos[2]))/2000)*1000)+'_'+str(int(((int(pos[1])+int(pos[2]))/2000+1))*1000) # only peak center matters
           if id in gc_anno:
               all_peak[int(gc_anno[id][0])] = gc_anno[id][1]
               for i in range(1, len(line)):
                   if line[i] != '0':
                       sc_peak[cells[i-1]].append(int(gc_anno[id][0]))  

    inf = open(tfindex_file, 'r')
    tfindex = json.load(inf)
    inf.close()
    tf_peak = {}
    total_peak = set(all_peak.keys())
    total_peak_l = len(set(all_peak.keys()))
    for k in tfindex.keys():
        if re.search("[a-zA-Z]+",k):
            max_op = 0
            max_id = 'id'
            for dataset in tfindex[k]:
                op = len(set(tfindex[str(dataset)])&total_peak)*1.0/total_peak_l
                if op > max_op:
                   max_op = op
                   max_id = str(dataset)
            if max_op > tfpercent:                                # only keep those TFs with at least certain percent overlap with total peaks
               tf_peak[k] = set(tfindex[max_id])

    return(all_peak, sc_peak, tf_peak)

def generate_bg_peaks(all_peak, sc_peak):
    """ Generate the background peaks based on the gc content. """

    anno_db = {}
    for k in all_peak.keys():
        if not all_peak[k] in anno_db:
           anno_db[all_peak[k]] = [k]
        else:
           anno_db[all_peak[k]].append(k)

    bg_peaks = []
    for peak in sc_peak:
        bg_peaks.append(int(random.sample(anno_db[all_peak[peak]],1)[0]))   # generate random background peaks with same gc content with sc peaks
    return set(bg_peaks)

def FP(args):
    """Multiple processing function to calculate fisher-exact test p-value."""
    
    tf_peak, tf_list, all_peak, sc_peak, cell = args
    sc_peak_set = set(sc_peak[cell])
    sc_len = len(sc_peak_set)
    bg_peak_set = generate_bg_peaks(all_peak, sc_peak_set)
    bg_len = len(bg_peak_set)
    pvalue = []
    for tf in tf_list:
        sc_tf = len(sc_peak_set & tf_peak[tf])
        bg_tf = len(bg_peak_set & tf_peak[tf])
        pvalue.append(str(-1*math.log(stats.fisher_exact([[sc_len - sc_tf, sc_tf], [bg_len - bg_tf, bg_tf]])[1]+1E-300,10)))
    
    return(pvalue)

def calculate_fisher_pvalue(score_file, all_peak, sc_peak, tf_peak, cores):
    """ Calculate the TF enrichment in single-cells based on CistomeDB. """
    
    cell_list = list(sc_peak.keys())
    tf_list = list(tf_peak.keys())
    
    args = ((tf_peak, tf_list, all_peak, sc_peak, cell) for cell in cell_list)
    pool = mp.Pool(processes = cores)
    result = pool.map_async(FP, args)
    pool.close()
    pool.join()
    pvalue_list = result.get()    
    
    outf = open(score_file, 'w')
    print("\t".join(cell_list), file=outf)
    for i in range(0,len(tf_list)):
        outline = [tf_list[i]]
        for j in range(0,len(cell_list)):
            outline.append(pvalue_list[j][i])
        print("\t".join(outline), file=outf)
    outf.close()

def main():

    peak_file = sys.argv[1]
    score_file = sys.argv[2]
    gc_file = sys.argv[3]
    tfindex_file = sys.argv[4]
    tfpercent = float(sys.argv[5])
    cores = int(sys.argv[6])

    start = time.time()
    all_peak, sc_peak, tf_peak = filter_tf_peaks(peak_file, gc_file, tfindex_file, tfpercent)
    calculate_fisher_pvalue(score_file, all_peak, sc_peak, tf_peak, cores)
    end = time.time()
    print("TFScore Time:", end-start)

if __name__ == "__main__":
    main()

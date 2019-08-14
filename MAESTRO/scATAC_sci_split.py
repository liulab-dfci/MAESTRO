#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 15:34:21 2019

@author: Dongqing Sun
"""
import argparse
import numpy as np
import multiprocessing as mp
import pysam
import time, os

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("-S","--sam", dest = "samfile",default = "",help = "The merged samfile.")
    parser.add_argument("-P","--process",dest = "process_num",default = 1,help = "Number of processes. Default: 1")
    parser.add_argument("-B","--barcode",dest = "barcode_library",default = "",help = "The barcode library, in which the first column is the valid barcode combinations included in experiment.")
    return parser.parse_args()


# calculate the hamming distance between two sequences
def hamming_distance(seq1,seq2):
    ham_dis = 0
    for ch1,ch2 in zip(seq1,seq2):
        if ch1 != ch2:
            ham_dis += 1
    return ham_dis


# determine the real barcode
def assign_barcode(barcode):
    if barcode in barcode_lib_list:
        valid_barcode = barcode
    else:
        barcode_dis_list = [hamming_distance(barcode, bc) for bc in barcode_lib_list]
        if min(barcode_dis_list) <= 3:
            barcode_similar_index = np.where(np.array(barcode_dis_list) == min(barcode_dis_list))
            barcode_similar = [barcode_lib_list[i] for i in barcode_similar_index[0].tolist()]
            valid_barcode = barcode_similar
        else:
            valid_barcode = ""
    return (barcode, valid_barcode)

# split the sam file according to barcodes
def split(index_se):
    for bc in list(barcode_dict)[index_se[0]:index_se[1]]:
        fuzzy_barcode_list = barcode_dict[bc]
        read_list = []
        for barcode in fuzzy_barcode_list:
            split_samfile_in = open("Result/BWA/tmp/" + barcode + ".sam","r").readlines()
            read_list += [read.strip() for read in split_samfile_in]
        if len(read_list) > 200:
            read_str = "\n".join(read_list)
            outfile = open("Result/BWA/Split/" + bc + ".sam","w")
            outfile.write(str(samfile_in.header))
            outfile.write(read_str)
            outfile.close()
        else:
            continue

parser = CommandLineParser()
process_num = int(parser.process_num)
samfile = parser.samfile
barcode_file = parser.barcode_library


start_time = time.time()
print("Start to split bam file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
samfile_in = pysam.AlignmentFile(samfile,"rb")
barcode_list = []
os.mkdir("Result/BWA/tmp/")
for read in samfile_in:
    barcode = read.query_name.strip().split("\t")[0].split(":")[1].split("-")[-1]
    out_bam = open("Result/BWA/tmp/" + barcode + ".sam", "a")
    out_bam.write(read.to_string() + "\n")
    out_bam.close()
    barcode_list.append(barcode)
end_time = time.time()
print("End:", end_time-start_time)


start_time = time.time()
print("Start to read barcode library",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
barcode_lib_list = open(barcode_file).readlines()
barcode_lib_list = [bc.strip().split("\t")[0] for bc in barcode_lib_list]
end_time = time.time()
print("End:", end_time-start_time)

barcode_list = list(set(barcode_list))

start_time = time.time()
print("Start to assign barcodes",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
process_num = 8
pool = mp.Pool(processes=process_num)
result = pool.map_async(assign_barcode, barcode_list)
pool.close()
pool.join()
end_time = time.time()
print("End:", end_time-start_time)

barcode_pair_list = result.get()
barcode_dict = {}
for pair in barcode_pair_list:
    if isinstance(pair[1], str):
        if pair[1]:
            if pair[1] not in barcode_dict:
                barcode_dict[pair[1]] = [pair[0]]
            else:
                barcode_dict[pair[1]].append(pair[0])
        else:
            continue
    else:
        for bc in pair[1]:
            if bc not in barcode_dict:
                barcode_dict[bc] = [pair[0]]
            else:
                barcode_dict[bc].append(pair[0])

start_time = time.time()
print("Start to split sam files according to barcode",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
barcode_num = len(barcode_dict.keys())
i = 0
index_list = []
while i < barcode_num:
    if i + int(barcode_num/process_num) + 1 < barcode_num:
        index_list.append((i, i + int(barcode_num/process_num) + 1))
    else:
        index_list.append((i, barcode_num))
    i = i + int(barcode_num/process_num) + 1

pool = mp.Pool(processes=process_num)
result = pool.map_async(split, index_list)
pool.close()
pool.join()
end_time = time.time()
print("End:", end_time-start_time)

os.system("rm -r Result/BWA/tmp/")

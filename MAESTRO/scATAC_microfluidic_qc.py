#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 13:44:50 2019

@author: ChenfeiWang, Dongqing Sun
"""

import sys,os

def main():

    log_folder = sys.argv[1]
    out_file = sys.argv[2]
    output = {}
    for logfile in os.listdir(log_folder):
        if logfile.endswith('.mapping.log'):
            sample = logfile[:-12]
            total,mapped,duplicate,uniq,mito,promoters,peaks = 0,0,0,0,0,0,0
            line_id = 0
            for line in open(log_folder + logfile).readlines():
                line = line.strip().split(' ')
                line_id += 1
                if line_id == 1:
                    total = line[0]
                if line_id == 5:
                    mapped = line[0]
                if line_id == 4:
                    duplicate = line[0]
                if line_id == 14:
                    uniq = line[0]
                if line_id == 15:
                    mito = line[0]
                if line_id == 16:
                    promoters = line[0]
                if line_id == 17:
                    peaks = line[0]
            output[sample] = [total,mapped,duplicate,uniq,mito,promoters,peaks]

    outf = open(out_file,'w')
    outf.write('sample\ttotal\tmapped\tduplicate\tuniq\tmito\tpromoters\tpeaks\n')
    for k in output:
        outf.write(k+'\t'+'\t'.join(output[k])+'\n')
    outf.close()

if __name__ == "__main__":
    main()
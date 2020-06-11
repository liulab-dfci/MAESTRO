# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-06-07 16:22:11
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-06-07 18:30:51


import argparse
import pysam
import time, os
from collections import defaultdict


def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("-B","--bam", dest = "bamfile",default = "",help = "The merged bamfile.")
    parser.add_argument("-T","--tagfile",dest = "tagfile",default = "",help = "The tag file.")
    parser.add_argument("-C","--cbtag",dest = "cbtag",default = "",help = "The tag of cell barcode.")
    parser.add_argument("-O", "--outdir", dest = "outdir",default = "",help = "The output directory.")
    parser.add_argument("-P", "--prefix", dest = "prefix",default = "",help = "The output prefix.")

    return parser.parse_args()

parser = CommandLineParser()
bamfile = parser.bamfile
tagfile = parser.tagfile
cbtag = parser.cbtag
outdir = parser.outdir
prefix = parser.prefix
bamfile_new = os.path.join(outdir, prefix + ".bam")

# Read barcode correct file
start_time = time.time()
print("Start to read tag file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
barcode_tag_dict = {}
with open(tagfile, "r") as tagfile_in:
    for line in tagfile_in:
        line_list = line.strip().split("\t")
        barcode_tag_dict[line_list[0]] = (line_list[1], line_list[2])
end_time = time.time()
print("End", end_time-start_time)


# Add barcode
start_time = time.time()
print("Start to add tag",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
bamfile_in = pysam.AlignmentFile(bamfile,"rb")
with pysam.AlignmentFile(bamfile_new, "wb", template=bamfile_in) as bamfile_out:
    for read in bamfile_in:
        barcode = read.get_tag(cbtag)
        flag = read.flag
        if barcode in barcode_tag_dict:
            read.set_tag(barcode_tag_dict[barcode][0], barcode_tag_dict[barcode][1], value_type="Z")
            bamfile_out.write(read)
    bamfile_out.close()
end_time = time.time()
print("End:", end_time-start_time)

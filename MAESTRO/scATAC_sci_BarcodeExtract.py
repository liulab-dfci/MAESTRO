# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-02 13:42:10
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-02-03 23:11:23


import argparse
import pysam
import time, os
import re
from collections import defaultdict


def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("--R1", dest = "r1_in", default = "", help = "The original r1 fastq file")
    parser.add_argument("--R2", dest = "r2_in", default = "", help = "The original r2 fastq file")
    parser.add_argument("-O", "--outdir", dest = "outdir",default = "",help = "The output directory.")
    return parser.parse_args()


parser = CommandLineParser()
r1_in_file = parser.r1_in
r2_in_file = parser.r2_in
outdir = parser.outdir
prefix = r1_in_file.split("/")[-1].split(".")[0].split("_")[0]
r1_out_file = os.path.join(outdir, prefix + "_R1.fastq")
r2_out_file = os.path.join(outdir, prefix + "_R2.fastq")
r3_out_file = os.path.join(outdir, prefix + "_R3.fastq")


start_time = time.time()
print("Start to pre-process fastq file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
with pysam.FastxFile(r1_in_file) as r1_in, pysam.FastxFile(r2_in_file) as r2_in, open(r1_out_file, "w") as r1_out, open(r2_out_file, "w") as r2_out, open(r3_out_file, "w") as r3_out:
    for entry in r1_in:
        r3_entry = next(r2_in)
        name = entry.name
        sequence = entry.sequence
        comment = entry.comment
        quality = entry.quality
        barcode = name.split(":")[1].split("-")[-1]
        barcode_entry = pysam.FastxRecord()
        # write r1 
        r1_out.write(str(entry) + "\n")
        # write r2
        barcode_entry.name = name
        barcode_entry.sequence = barcode
        barcode_entry.comment = comment
        barcode_entry.quality = quality[0:len(barcode)]
        r2_out.write(str(barcode_entry) + "\n")
        # write r3
        r3_out.write(str(r3_entry) + "\n")

end_time = time.time()
print("End", end_time-start_time)


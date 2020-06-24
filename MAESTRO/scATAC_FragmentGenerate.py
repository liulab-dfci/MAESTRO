# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-01-17 01:43:08
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-06-10 23:32:28

import argparse
import pysam
import time, os
from collections import defaultdict

chr_list = ['chr'+str(i) for i in list(range(1,23))] 
chr_list = chr_list + ['chrX', 'chrY']

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("-B","--bam", dest = "bamfile", default = "",help = "The merged bamfile.")
    # parser.add_argument("-b","--barcodefq",dest = "barcode_fastq",default = "",help = "The barcode fastq file (for sciATAC or 10x genomics, R2)")
    parser.add_argument("-O", "--outdir", dest = "outdir", default = "", help = "The output directory.")
    parser.add_argument("--addtag", dest = "addtag", default = "", help = "The tag cell barcodes will be added to. If not set, no tag will be added.")
    parser.add_argument("--CBtag", dest = "CBtag", default = "", help = "Where the cell identities are stored. For example, 'CB'. "
        "If not set, cell identities will be inferred from read name.")
    parser.add_argument("--count", dest = "count", action = "store_true", help = "Whether to count the fragments.")

    return parser.parse_args()


parser = CommandLineParser()
bamfile = parser.bamfile
# barcode_fq = parser.barcode_fastq
outdir = parser.outdir
addtag = parser.addtag
CBtag = parser.CBtag
count = parser.count
fragment_file = os.path.join(outdir, "fragments.tsv")

try:
    os.makedirs(outdir)
except OSError:
    pass

# start_time = time.time()
# print("Start to read R2 barcode file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
# barcodefile_in = pysam.FastxFile(barcode_fq)
# read_barcode_dict = {}
# for barcode in barcodefile_in:
#     sequence = barcode.sequence
#     read_barcode_dict[barcode.name] = sequence
# end_time = time.time()
# print("End", end_time-start_time)


start_time = time.time()
if addtag != "":
    bamfile_new = bamfile.split(".bam")[0] + "." + addtag + "added.bam"
    print("Start to read mapping bam file, add tag and generate fragment file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    bamfile_in = pysam.AlignmentFile(bamfile,"rb")
    fragment_out = open(fragment_file, "w")
    with pysam.AlignmentFile(bamfile_new, "wb", template=bamfile_in) as bamfile_out:
        if CBtag == "":
            for read in bamfile_in:
                barcode = read.query_name.split(":")[0]
                flag = read.flag
                read.set_tag(addtag, barcode, value_type="Z")
                bamfile_out.write(read)
                if read.reference_name in chr_list and (flag & 0x2 !=0) and (flag & 0xc == 0) and (flag & 0x900 == 0) and read.mapping_quality > 30 and read.template_length < 1000 and read.template_length > 10:
                    frag_list = [read.reference_name, str(read.reference_start + 4), str(read.reference_start + read.template_length - 5), barcode]
                    fragment_out.write("\t".join(frag_list) + "\n")
            fragment_out.close()
            bamfile_out.close()
        else:
            for read in bamfile_in:
                if read.has_tag(CBtag):
                    barcode = read.get_tag(CBtag)
                    flag = read.flag
                    read.set_tag(addtag, barcode, value_type="Z")
                    bamfile_out.write(read)
                    if read.reference_name in chr_list and (flag & 0x2 !=0) and (flag & 0xc == 0) and (flag & 0x900 == 0) and read.mapping_quality > 30 and read.template_length < 1000 and read.template_length > 10:
                        frag_list = [read.reference_name, str(read.reference_start + 4), str(read.reference_start + read.template_length - 5), barcode]
                        fragment_out.write("\t".join(frag_list) + "\n")
            fragment_out.close()
            bamfile_out.close()
else:
    print("Start to read mapping bam file, and generate fragment file",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    bamfile_in = pysam.AlignmentFile(bamfile,"rb")
    fragment_out = open(fragment_file, "w")
    if CBtag == "":
        for read in bamfile_in:
            barcode = read.query_name.split(":")[0]
            flag = read.flag
            if read.reference_name in chr_list and (flag & 0x2 !=0) and (flag & 0xc == 0) and (flag & 0x900 == 0) and read.mapping_quality > 30 and read.template_length < 1000 and read.template_length > 10:
                frag_list = [read.reference_name, str(read.reference_start + 4), str(read.reference_start + read.template_length - 5), barcode]
                fragment_out.write("\t".join(frag_list) + "\n")
        fragment_out.close()
    else:
        for read in bamfile_in:
            if read.has_tag(CBtag):
                barcode = read.get_tag(CBtag)
                flag = read.flag
                if read.reference_name in chr_list and (flag & 0x2 !=0) and (flag & 0xc == 0) and (flag & 0x900 == 0) and read.mapping_quality > 30 and read.template_length < 1000 and read.template_length > 10:
                    frag_list = [read.reference_name, str(read.reference_start + 4), str(read.reference_start + read.template_length - 5), barcode]
                    fragment_out.write("\t".join(frag_list) + "\n")
        fragment_out.close()
end_time = time.time()
print("End:", end_time-start_time)

if count:
    fragment_file_sorted = os.path.join(outdir, "fragments_corrected_sorted.tsv")
    fragment_file_count = os.path.join(outdir, "fragments_corrected_count.tsv")
    start_time = time.time()
    print("Start to count fragment",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    cmd = "sort -k1,1 -k2,2 -k3,3 -k4,4 -V %s > %s;bedtools groupby -i %s -g 1,2,3,4 -c 4 -o count > %s" %(fragment_file, fragment_file_sorted, fragment_file_sorted, fragment_file_count)
    os.system(cmd)
    end_time = time.time()
    print("End:", end_time-start_time)

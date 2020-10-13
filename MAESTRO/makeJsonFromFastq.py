#!/usr/bin/env python3


import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description = "This script is used to generate a json file from a folder containing \
	the fastq files for multi-sample processing. The fastq name format should be in a format similar to 10x. \
	e.g. If the file name is atac_pbmc_1k_v1_S1_L001_I1_001.fastq.gz, atac_pbmc_1k will be used as the sample name")
parser.add_argument("data_type", help="Required. type of data: either scatac or scrnaseq")
parser.add_argument("fastq_dir", help="Required. the FULL path to the fastq folder")
args = parser.parse_args()

assert args.data_type == "scrnaseq" or args.data_type == "scatac", "please specify scatac or scrnaseq for data_type"
assert args.fastq_dir is not None, "please provide the path to the fastq folder"


## default dictionary is quite useful!

FILES = defaultdict(lambda: defaultdict(list))

## build the dictionary with full path for each fastq.gz file
for root, dirs, files in os.walk(args.fastq_dir):
	for file in files:
		if file.endswith("fastq.gz"):
			full_path = join(root, file)
			if args.data_type == "scrnaseq":
				#R1 will be sample barcode, R2 will be reverse reads, I1 will be the index 
				m = re.search(r"([A-Z0-9a-z_]+)_S[0-9]_(L[0-9]{3})_([IR][12])_[0-9]+.fastq.gz", file)
				if m:
					sample = m.group(1)
					lane = m.group(2)
					reads = m.group(3)  
					FILES[sample][reads].append(full_path)
			elif args.data_type == "scatac":
				m = re.search(r"([A-Z0-9a-z_]+)_S[0-9]_(L[0-9]{3})_([IR][123])_[0-9]+.fastq.gz", file)
				if m:
					sample = m.group(1)
					lane = m.group(2)
					# I1 will be sample index, R1 and R3 for forward and reverse read, R2 is the cell barcode
					reads = m.group(3)  
					FILES[sample][reads].append(full_path)

			

# make sure R1 and R2 from different lanes are ordered in the same way
# e.g. L001_R1 should pair with L001_R2, L002_R1 pair with L002_R2        

FILES_sorted = defaultdict(lambda: defaultdict(list))

for sample in FILES.keys():
		for read in FILES[sample]:
			FILES_sorted[sample][read] = sorted(FILES[sample][read])



print()
print ("total {} unique samples will be processed".format(len(FILES.keys())))
print ("------------------------------------------")
for sample in FILES_sorted.keys():
	for read in sorted(FILES_sorted[sample]):
		print ("{sample} {read} has {n} fastq".format(sample = sample, read = read, n = len(FILES_sorted[sample][read])))
print ("------------------------------------------")
print("check the samples.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES_sorted, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)


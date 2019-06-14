#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 17:26:44 2019

@author: Chenfei Wang, Dongqing Sun
"""

import argparse
import sys

def rmLinker_parseargs():
	"""
	Parsing arguments for rmLinker
	"""
	parser = argparse.ArgumentParser(description='Remove linkers from barcode sequence.')
	reqgroup = parser.add_argument_group(title='Required arguments', description="Required")
	reqgroup.add_argument('-f','--fastq', help="Fastq file which store barcode reads")
	reqgroup.add_argument('-o','--output', help="Name of output fastq file")
	reqgroup.add_argument('-s','--start', help="Start position of barcode, separate by comma")
	reqgroup.add_argument('-e','--end', help="End position of barcode, separate by comma")

	args = parser.parse_args()
	args.start = args.start.split(",")
	args.end = args.end.split(",")
	return(args)

def rmLinker():
	"""
	Filter barcode sequences from sequencing reads
	"""
	args=rmLinker_parseargs()
	inp = open(args.fastq, "r")
	out = open(args.output, "w")

	for line in inp.readlines():
		if not (line.startswith("@") | line.startswith("+")):
			i = 0
			tmp = ""
			while i<len(args.start):
				tmp = tmp + line[int(args.start[i]):int(args.end[i])]
				i = i + 1
			line = tmp + "\n"
		out.write(line)

	inp.close()
	out.close()


if __name__ == '__main__':
	try:
		rmLinker()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) Bye!\n")
		sys.exit(0)

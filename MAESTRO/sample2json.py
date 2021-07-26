# -*- coding: utf-8 -*-
# @Author: Ming Tang, Gali Bai
# @E-mail: gali.bai@hotmail.com
# @Date:   2020-02-23 19:40:27
# @Last Modified by:   Gali Bai
# @Last Modified time: 2021-06-09 17:34:46

import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict
import sys

def sample_parser(subparsers):
    """
    Add main function init-sample argument parsers.
    """
    workflow = subparsers.add_parser("samples-init", help = "Initialize samples.json file in the current directory.")

    group_input = workflow.add_argument_group("Input files arguments")
    group_input.add_argument("--assay_type", dest = "assay_type", choices = ["scrna", "scatac"],
        help="Required: type of assays.")
    group_input.add_argument("--platform", dest = "platform", default = "10x-genomics",
        choices = ["10x-genomics","sci-ATAC-seq","microfluidic"],
        help="Required: type of sequencing platform")
    group_input.add_argument("--data_type", dest = "data_type", choices = ["fastq", "fragment", "bam"],
        help="Required: type of data.")
    group_input.add_argument("--data_dir", dest = "data_dir",
        help="Required. the FULL path to the fastq folder or the fragment folder")

def sample_json(args):
    """
    Generate samples.json file.
    """
    assert args.assay_type == "scrna" or args.assay_type == "scatac", "please specify scatac or scrna for assay_type"
    assert args.data_type == "fastq" or args.data_type == "fragment" or args.data_type == "bam", "please specify fastq, fragment or bam for data_type"
    assert args.data_dir is not None, "please provide the path to the fastq or fragment folder"

    if args.assay_type == "scrna":
        if not args.data_type == "fastq":
            print("for scrna assay type, only fastq data type is supported")
            sys.exit(1)
    elif args.assay_type == "scatac":
        if (not args.data_type == "fragment") and (not args.data_type == "fastq"):
            print("for scatac assay type, only fastq or fragment data type are supported")
            sys.exit(1)

    ## default dictionary is quite useful!

    ## build the dictionary with full path for each fastq.gz file
    if args.data_type == "fastq":
        FILES = defaultdict(lambda: defaultdict(list))
        for root, dirs, files in os.walk(args.data_dir):
            for file in files:
                if file.endswith("fastq.gz"):
                    full_path = join(root, file)
                    if args.assay_type == "scrna":
                        #R1 will be sample barcode, R2 will be reverse reads, I1 will be the index
                        m = re.search(r"([A-Z0-9a-z_]+)_S[0-9]+_(L[0-9]{3})_([IR][12])_[0-9]+.fastq.gz", file)
                        if m:
                            sample = m.group(1)
                            lane = m.group(2)
                            reads = m.group(3)
                            FILES[sample][reads].append(full_path)
                    elif args.assay_type == "scatac" and args.platform == "10x-genomics":
                        m = re.search(r"([A-Z0-9a-z_]+)_S[0-9]+_(L[0-9]{3})_([IR][123])_[0-9]+.fastq.gz", file)
                        if m:
                            sample = m.group(1)
                            lane = m.group(2)
                            # I1 will be sample index, R1 and R3 for forward and reverse read, R2 is the cell barcode
                            reads = m.group(3)
                            FILES[sample][reads].append(full_path)
                    elif args.assay_type == "scatac" and args.platform == "microfluidic":
                        m = re.search(r"([A-Z0-9a-z_]+)_([12])+.fastq.gz",file)
                        if m:
                            sample = m.group(1)
                            reads = m.group(2)
                            FILES[sample][reads].append(full_path)
                        else:
                            print("Please provide fastq file with format prefix_1.fastq.gz and prefix_2.fastq.gz")
                    elif args.assay_type == "scatac" and args.platform == "sci-ATAC-seq":
                        m = re.search(r"([A-Z0-9a-z_]+)_([123])+.fastq.gz",file)
                        if m:
                            sample = m.group(1)
                            reads = m.group(2)
                            FILES[sample][reads].append(full_path)
                        else:
                            print("Please provide fastq file with format prefix_R1.fastq.gz, prefix_R2.fastq.gz (barcode) and prefix_R3.fastq.gz")


                elif file.endswith("fastq"):
                    full_path = join(root, file)
                    if args.assay_type == "scrna":
                        #R1 will be sample barcode, R2 will be reverse reads, I1 will be the index
                        m = re.search(r"([A-Z0-9a-z_]+)_S[0-9]+_(L[0-9]{3})_([IR][12])_[0-9]+.fastq", file)
                        if m:
                            sample = m.group(1)
                            lane = m.group(2)
                            reads = m.group(3)
                            FILES[sample][reads].append(full_path)
                    elif args.assay_type == "scatac" and args.platform == "10x-genomics":
                        m = re.search(r"([A-Z0-9a-z_]+)_S[0-9]+_(L[0-9]{3})_([IR][123])_[0-9]+.fastq", file)
                        if m:
                            sample = m.group(1)
                            lane = m.group(2)
                            # I1 will be sample index, R1 and R3 for forward and reverse read, R2 is the cell barcode
                            reads = m.group(3)
                            FILES[sample][reads].append(full_path)
                    elif args.assay_type == "scatac" and args.platform == "microfluidic":
                        m = re.search(r"([A-Z0-9a-z_]+)_([12])+.fastq",file)
                        if m:
                            sample = m.group(1)
                            reads = m.group(2)
                            FILES[sample][reads].append(full_path)
                        else:
                            print("Please provide fastq file with format prefix_1.fastq and prefix_2.fastq")
                    elif args.assay_type == "scatac" and args.platform == "sci-ATAC-seq":
                        m = re.search(r"([A-Z0-9a-z_]+)_([123])+.fastq",file)
                        if m:
                            sample = m.group(1)
                            reads = m.group(2)
                            FILES[sample][reads].append(full_path)
                        else:
                            print("Please provide fastq file with format prefix_R1.fastq, prefix_R2.fastq (barcode) and prefix_R3.fastq")
        FILES_sorted = defaultdict(lambda: defaultdict(list))
        # make sure R1 and R2 from different lanes are ordered in the same way
        # e.g. L001_R1 should pair with L001_R2, L002_R1 pair with L002_R2
        for sample in FILES.keys():
            for read in FILES[sample]:
                FILES_sorted[sample][read] = sorted(FILES[sample][read])

    elif args.data_type == "fragment":
        FILES = defaultdict(str)
        for root, dirs, files in os.walk(args.data_dir):
            for file in files:
                if file.endswith("tsv.gz"):
                    full_path = join(root, file)
                    if args.assay_type == "scatac":
                        #R1 will be sample barcode, R2 will be reverse reads, I1 will be the index
                        m = re.search(r"([A-Z0-9a-z_]+)_fragments.tsv.gz", file)
                        if m:
                            sample = m.group(1)
                            FILES[sample] = full_path
                        else:
                            sample = file.replace('.tsv.gz')
                            FILES[sample] = full_path
                elif file.endswith("tsv"):
                    full_path = join(root, file)
                    if args.assay_type == "scatac":
                        m = re.search(r"([A-Z0-9a-z_]+)_fragments.tsv", file)
                        if m:
                            sample = m.group(1)
                            FILES[sample] = full_path
                        else:
                            sample = file.replace('.tsv')
                            FILES[sample] = full_path
                else:
                    print("please provide correct fragment file with .tsv or .tsv.gz as extensions.")
    elif args.data_type == "bam":
        FILES = defaultdict(str)
        for root, dirs, files in os.walk(args.data_dir):
            for file in files:
                if file.endswith(".bam"):
                    full_path = join(root, file)
                    m = re.search(r"([A-Z0-9a-z_]+).bam", file)
                    if m:
                        sample = m.replace('.bam')
                        FILES[sample] = full_path
                    else:
                        print("Please provide bam file with correct format.")

    if args.data_type == "fastq":
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
    elif args.data_type =="fragment":
        print()
        print ("total {} unique samples will be processed".format(len(FILES.keys())))
        print ("------------------------------------------")
        for sample in FILES.keys():
            print ("{sample}'s fragment file is {fragment}".format(sample = sample, fragment = FILES[sample]))
        print ("------------------------------------------")
        print("check the samples.json file for fragment file belongs to each sample")
        print()
        js = json.dumps(FILES, indent = 4, sort_keys=True)
        open('samples.json', 'w').writelines(js)
    elif args.data_type =="bam":
        print()
        print ("total {} unique samples will be processed".format(len(FILES.keys())))
        print ("------------------------------------------")
        for sample in FILES.keys():
            print ("{sample}'s bam file is {bam}".format(sample = sample, bam = FILES[sample]))
        print ("------------------------------------------")
        print("check the samples.json file for bam file belongs to each sample")
        print()
        js = json.dumps(FILES, indent = 4, sort_keys=True)
        open('samples.json', 'w').writelines(js)

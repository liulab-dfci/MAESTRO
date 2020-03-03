#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 22:36:07 2019

@author: Dongqing Sun, Chenfei Wang
"""

import os, re
from pkg_resources import resource_filename

SCRIPT_PATH = os.path.dirname(__file__)
RSCRIPT_PATH = resource_filename('MAESTRO', 'R')


def getfastq_10x(fastqdir, fastqprefix):
    allfiles = os.listdir(fastqdir)
    fastqfiles = []
    pattern = fastqprefix + "\S*.fastq\S*"
    for file in allfiles:
        if re.match(pattern, file):
            fastqfiles.append(file)
        else:
            pass
    r1_fastq = []
    pattern = "\S*R1\S*"
    for file in fastqfiles:
        if re.search(pattern, file):
            r1_fastq.append(file)
        else:
            pass
    r1_fastq = sorted(r1_fastq)
    r2_fastq = list(map(lambda x: re.sub("_R1_", "_R2_", x), r1_fastq))
    transcript = ""
    barcode = ""
    decompress = "-"
    if len(set(r2_fastq) & set(fastqfiles)) == len(r1_fastq):
        r2_fastq = list(map(lambda x: os.path.join(fastqdir, x), r2_fastq))
        r1_fastq = list(map(lambda x: os.path.join(fastqdir, x), r1_fastq))
        transcript = ",".join(r2_fastq)
        barcode = ",".join(r1_fastq)
    else:
        print("Invalid fastq files!")
    if r1_fastq[0].endswith("gz"):
        decompress = "zcat"
    else:
        decompress = "-"

    return({"transcript": transcript, "barcode": barcode, "decompress": decompress})

def getfastq_dropseq(fastqdir, barcode, transcript):
    barcode_files = split(barcode, ",")
    transcript_files = split(transcript, ",")
    barcode_files = sorted(barcode_files)
    transcript_files = sorted(transcript_files)
    
    transcript = ""
    barcode = ""
    decompress = "-"
    if len(set(transcript_files)) == len(set(barcode_files)):
        barcode_files = list(map(lambda x: os.path.join(fastqdir, x), barcode_files))
        transcript_files = list(map(lambda x: os.path.join(fastqdir, x), transcript_files))
        transcript = ",".join(transcript_files)
        barcode = ",".join(barcode_files)
    else:
        print("Invalid fastq files!")

    if barcode_files[0].endswith("gz"):
        decompress = "zcat"
    else:
        decompress = "-"

    return({"transcript": transcript, "barcode": barcode, "decompress": decompress})

def get_fastqfile(fastqpath):
    files = os.listdir(fastqpath)
    fastq_1 = []
    fastq_2 = []
    fastqs = []
    for f in files:
        if f.endswith('.fastq'):
            fastqs.append(f)
            if f.endswith('_1.fastq'):
                fastq_1.append(fastqpath + f)
            elif f.endswith('_2.fastq'):
                fastq_2.append(fastqpath + f)

    fastq_1 = sorted(fastq_1)
    fastq_2 = sorted(fastq_2)
    fastqs = sorted(fastqs)
    
    fastqstr = ''
    if len(fastq_1) != 0:
        if len(fastq_1) == len(fastq_2) and len(fastqs) == len(fastq_1) + len(fastq_2):
            fastqstr = ','.join(fastq_1) + ' ' + ','.join(fastq_2)
        else:
            print("Invalid fastq files!")
    else:
        if len(fastq_1) == len(fastq_2) and len(fastqs) != 0:
            fastqstr = ','.join(fastqs)
        else:
            print("Invalid fastq files!")

    return(fastqstr)

def get_fastqid(fastqpath):
    files = os.listdir(fastqpath)
    fastq_1 = []
    fastq_2 = []
    fastqs = []
    for f in files:
        if f.endswith('.fastq'):
            fastqs.append(f)
            if f.endswith('_1.fastq'):
                fastq_1.append(f)
            elif f.endswith('_2.fastq'):
                fastq_2.append(f)

    fastq_1 = sorted(fastq_1)
    fastq_2 = sorted(fastq_2)
    fastqs = sorted(fastqs)
    
    fastqidstr = ''
    if len(fastq_1) != 0:
        if len(fastq_1) == len(fastq_2) and len(fastqs) == len(fastq_1) + len(fastq_2):
            samplelist = ['ID:'+ i[0:len(i)-8] for i in fastq_1]
            fastqidstr = ' , '.join(samplelist)
        else:
            print("Invalid fastq files!")
    else:
        if len(fastq_1) == len(fastq_2) and len(fastqs) != 0:
            samplelist = ['ID:' + i[0:len(i)-6] for i in fastqs]
            fastqidstr = ' , '.join(fastqs)
        else:
            print("Invalid fastq files!")

    return(fastqidstr)

def get_fastqlist(fastqpath):
    files = os.listdir(fastqpath)
    fastq_1 = []
    fastq_2 = []
    fastqs = []
    for f in files:
        if f.endswith('.fastq'):
            fastqs.append(f)
            if f.endswith('_1.fastq'):
                fastq_1.append(f)
            elif f.endswith('_2.fastq'):
                fastq_2.append(f)

    fastq_1 = sorted(fastq_1)
    fastq_2 = sorted(fastq_2)
    fastqs = sorted(fastqs)
    samplelist = []
    if len(fastq_1) != 0:
        if len(fastq_1) == len(fastq_2) and len(fastqs) == len(fastq_1) + len(fastq_2):
            samplelist = [i[0:len(i)-8] for i in fastq_1]
        else:
            print("Invalid fastq files!")
    else:
        if len(fastq_1) == len(fastq_2) and len(fastqs) != 0:
            samplelist = [i[0:len(i)-6] for i in fastqs]
        else:
            print("Invalid fastq files!")

    return(samplelist)

def get_bamfile(bampath):
    files = os.listdir(bampath)
    bams = []
    for f in files:
        if f.endswith('Aligned.sortedByReads.out.bam'):
            bams.append(bampath + f)
    bams = sorted(bams)

    bam_str = ' '.join(bams)
    return(bam_str)

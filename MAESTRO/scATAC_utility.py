#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 22:36:07 2019

@author: Dongqing Sun, Chenfei Wang
"""

import sys, os, time
import logging
import subprocess
import random, string
import re
import gzip
from subprocess import call as subpcall
from pkg_resources import resource_filename

error   = logging.critical
warn    = logging.warning
ENV_PATH = resource_filename('MAESTRO', 'env')
SCRIPT_PATH = os.path.dirname(__file__)
RSCRIPT_PATH = resource_filename('MAESTRO', 'R')

def Info(infoStr):
    print("[%s] %s" %(time.strftime('%H:%M:%S'), infoStr))
    
def run_cmd(command):
    subpcall (command, shell = True)

def run_pip(command):
    command = subprocess.Popen(command,
                               shell = True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               )

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
    r2_fastq = list(map(lambda x: re.sub("_R1", "_R2", x), r1_fastq))
    r3_fastq = list(map(lambda x: re.sub("_R1", "_R3", x), r1_fastq))
    transcript = ""
    barcode = ""
    command = "-"
    if len(set(r2_fastq) & set(fastqfiles)) == len(r1_fastq) and len(set(r3_fastq) & set(fastqfiles)) == len(r1_fastq):
        r3_fastq = list(map(lambda x: os.path.join(fastqdir, x), r3_fastq))
        r2_fastq = list(map(lambda x: os.path.join(fastqdir, x), r2_fastq))
        r1_fastq = list(map(lambda x: os.path.join(fastqdir, x), r1_fastq))
        r3 = " ".join(r3_fastq)
        r2 = " ".join(r2_fastq)
        r1 = " ".join(r1_fastq)
    else:
        print("Invalid fastq files!")
    if is_gzip(r1_fastq[0]):
        command = "zcat"
    else:
        command = "cat"

    return({"r1": r1, "r3": r3, "barcode": r2, "command": command})


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

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def is_gzip ( filename ):
    """Check if the file is gzipped.
    """
    with gzip.open( filename, 'r' ) as f:
        try:
            f.read(1)
        except OSError:
            return False
    return True

def universal_open ( filename, mode ):
    if is_gzip( filename ):
        return gzip.open( filename, mode )
    else:
        return open( filename, mode )
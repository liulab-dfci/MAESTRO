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
from subprocess import call as subpcall
from pkg_resources import resource_filename

error   = logging.critical
warn    = logging.warning
ENV_PATH = resource_filename('strap', 'env')
SCRIPT_PATH = os.path.dirname(__file__)
RSCRIPT_PATH = resource_filename('strap', 'R')

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


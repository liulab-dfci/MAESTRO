#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 13:10:51 2019

@author: Chenfei Wang, Dongqing Sun
"""

import sys,os
from distutils.core import setup


# try:
#     from setuptools import setup, find_packages
# except ImportError:
#     print("Could not load setuptools. Please install the setuptools package.")

def install_drseq():
    curdir = os.getcwd()
    os.chdir("refpkg/Dr.seq.1.2.0")
    os.system("python setup.py install")
    os.chdir(curdir)
    print("Installation of Dr.seq is DONE")

def install_giggle():
    curdir = os.getcwd()
    os.chdir("refpkg/giggle")
    os.system("make")
    os.chdir(curdir)
    print("Installation of GIGGLE is DONE")

def install_rabit():
    curdir = os.getcwd()
    os.chdir("refpkg/Rabit")
    os.system("./configure --prefix=" + curdir)
    os.system("make")
    os.system("make install")
    os.chdir(curdir)
    print("Installation of Rabit is DONE")

def install_rpackage():
    os.system("Rscript MAESTRO/R/MAESTRO_install.R")
    print("Installation of required R packages is DONE")

def main():
    # install_drseq()
    install_rpackage()
    install_giggle()
    # install_rabit()
    
    setup(
        name = "MAESTRO",
        version = "1.1.0",
        package_dir = {'MAESTRO':'MAESTRO'},
        packages = ['MAESTRO'],
        package_data={'MAESTRO':['Snakemake/scRNA/*', 'Snakemake/scATAC/*', 'R/*', 'env/*', 'annotations/*', 'html/*', '']},
        data_files = [('bin', ['refpkg/giggle/bin/giggle'])],
        scripts = ['MAESTRO/MAESTRO'],
        include_package_data = True,
        
        author = "Chenfei Wang, Dongqing Sun",
        author_email = "",
        description = "MAESTRO(Model-based AnalysEs of Single-cell Transcriptome and RegulOme) is a comprehensive "
        "single-cell RNA-seq and ATAC-seq analysis suit built using snakemake.",
        license = "GPL-3.0",
        url = "https://github.com/chenfeiwang/MAESTRO",
        
        # entry_points = {"console_scripts": ["strap = strap:main"]},
        classifiers = [
            "Development Status :: 4 - Beta",
            #"Development Status :: 5 - Production/Stable",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GPL License",
            "Natural Language :: English",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ]
    )


if __name__ == "__main__":
    main()



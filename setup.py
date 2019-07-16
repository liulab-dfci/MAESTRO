#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 13:10:51 2019

@author: Chenfei Wang, Dongqing Sun
"""

import sys,os

try:
    from setuptools import setup, find_packages
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.")

def install_drseq():
    curdir = os.getcwd()
    os.chdir("pkg/Dr.seq.1.2.0")
    os.system("python setup.py install")
    os.chdir(curdir)

def install_giggle():
    curdir = os.getcwd()
    os.chdir("pkg/giggle")
    os.system("make")
    os.chdir(curdir)

def main():
    install_drseq()
    install_giggle()
    setup(
        name = "MASTER",
        version = "0.0.1",
        package_dir = {'MASTER':'MASTER'},
        packages = ['MASTER'],
        package_data={'MASTER':['Snakemake/scRNA/*', 'Snakemake/scATAC/*', 'R/*', 'env/*', 'annotations/*', '']},
        scripts = ['MASTER/MASTER'],
        include_package_data = True,
        
        author = "Chenfei Wang, Dongqing Sun",
        author_email = "",
        description = "MASTER (model-based analysis of single-cell transcriptome and regulome) is a comprehensive "
        "quality control, analysis and visualization pipeline for scATAC-seq and scRNA-seq",
        license = "GPL-3.0",
        url = "https://github.com/chenfeiwang/MASTER",
        
        zip_safe = False,
        install_requires=[],
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



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

def main():
    install_drseq()
    setup(
        name = "strap",
        version = "0.0.1",
        package_dir = {'strap':'strap'},
        packages = ['strap'],
        package_data={'strap':['Snakemake/scRNA/*', 'Snakemake/scATAC/*', 'R/*', 'env/*', 'annotations/*']},
        scripts = ['strap/strap'],
        include_package_data = True,
        
        author = "Chenfei Wang, Dongqing Sun",
        author_email = "",
        description = "STRAP (single-cell transcriptome and transcriptome analysis pipeline) is a comprehensive "
        "quality control, analysis and visualization pipeline for scATAC-seq and scRNA-seq",
        license = "MIT",
        url = "https://bitbucket.org/chenfeiwang/strap",
        
        zip_safe = False,
        install_requires=["h5py", "scipy"],
        # entry_points = {"console_scripts": ["strap = strap:main"]},
        classifiers = [
            "Development Status :: 4 - Beta",
            #"Development Status :: 5 - Production/Stable",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ]
    )


if __name__ == "__main__":
    main()



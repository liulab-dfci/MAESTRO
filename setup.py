# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-06-12 03:53:45
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-06-15 21:58:56


import sys,os
# from distutils.core import setup


try:
    from setuptools import setup, find_packages
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.")

def main():
    setup(
        name = "MAESTRO",
        version = "1.2.1",
        package_dir = {'MAESTRO':'MAESTRO'},
        packages = ['MAESTRO'],
        package_data={'MAESTRO':['Snakemake/scRNA/*', 'Snakemake/scATAC/*', 'Snakemake/integrate/*', 'R/*', 'annotations/*', 'html/*', '']},
        #data_files = [(os.path.join(sys.prefix,'bin'), ['refpkg/giggle/bin/giggle'])],
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
        ],
        #install_requires=["sinto>=0.7.1",],
        #setup_requires=["sinto>=0.7.1",],        
    )


if __name__ == "__main__":
    main()

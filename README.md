# MASTAR

**MASTAR**(**M**odel-based **A**nalysis of **S**ingle-cell **T**ranscriptome **A**nd **R**egulome) is a comprehensive scATAC-seq and scRNA-seq analysis tool built using [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) which allows for ease of use. It can apply to different platforms, such as Smart-seq2, drop-seq, SPLiT-seq, microwell-seq, Chang&Greenleaf and JayShedure protocol for ATAC-seq, and 10X gemomics. STRAP combines the use of several dozen tools, suites, and packages to create a complete pipeline that takes scATAC-seq and scRNA-seq analysis from raw sequencing data(fastq files or count 
table ) all the way through alignment, quality control, unsupervised analyses, differential expression,annotation, downstream analysis. The results are compiled in a labelled tSNE plot or heatmap.

!!![picture]

## System requirements:
Some of the tools that STRAP uses, e.g. STAR and Seurat are very memory intensive programs. Therefore we recommend the following system requirements for STRAP.

### Minimal system requirements:
We recommend that you run STRAP on a server that has at least 30GB of ram. This will allow for a single-threaded STRAP run (on human samples).

### Recommended system requirements:
We recommend that you have at least 128GB of ram and at least a 4-core CPU if you want to run STRAP n multi-threaded mode (which will speed up the workflow significantly). 

## Installation

### Installing Miniconda3:

We will be using the [Miniconda3](http://conda.pydata.org/miniconda.html) package management system to manage all of the software packages that __STRAP__ is dependent on. 

Use following commands to the install Minicoda3：
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
__NOTE__: you will only have to install Minicoda3 once.  

### Installing the STRAP conda environments**

We are now ready to use CONDA to install the software packages which STRAP is dependent on.
```
wget https://bitbucket.org/strap
tar -xf strap.tar.gz
cd strap/envs
conda env create -f environment.yml -n STRAP
```
__NOTE__: you will only have to install the STRAP conda environments once.


# MAESTRO

**MAESTRO**(**M**odel-based **A**nalys**E**s of **S**ingle-cell **T**ranscriptome and **R**egul**O**me) is a comprehensive scATAC-seq and scRNA-seq analysis tool built using [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) which allows for ease of use. It can apply to different platforms, such as Smart-seq2, drop-seq, SPLiT-seq, microwell-seq, Chang&Greenleaf and JayShedure protocol for ATAC-seq, and 10X gemomics. STRAP combines the use of several dozen tools, suites, and packages to create a complete pipeline that takes scATAC-seq and scRNA-seq analysis from raw sequencing data(fastq files or count 
table ) all the way through alignment, quality control, unsupervised analyses, differential expression,annotation, downstream analysis. The results are compiled in a labelled tSNE plot or heatmap.

## 0. System requirements
Some of the tools that STRAP uses, e.g. STAR and Seurat are very memory intensive programs. Therefore we recommend the following system requirements for STRAP.

**Minimal system requirements:**
We recommend that you run STRAP on a server that has at least **x**GB of ram. This will allow for a single-threaded STRAP run (on human samples).

**Recommended system requirements:**
We recommend that you have at least **x**GB of ram and at least a **x**-core CPU if you want to run STRAP n multi-threaded mode (which will speed up the workflow significantly). 


## 1.Installation

### 1.1 Installing Cell Ranger

__STRAP__ depent on the Cell Ranger for the mapping of the data genertaed by 10X genomic. So you need to install the Cell Ranger for the first step.Use this [link](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) to the install Cell Ranger.

### 1.2 Installing Miniconda3

We will be using the [Miniconda3](http://conda.pydata.org/miniconda.html) package management system to manage all of the software packages that __STRAP__ is dependent on. 

Use following commands to the install Minicoda3：

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

*NOTE*: you will only have to install Minicoda3 once.  

### 1.3 Installing the STRAP by conda

We are now ready to use CONDA to install the STRAP.

```
conda install -c dongqingsun strap
```

*NOTE*: you will only have to install the STRAP conda environments once.

## 2. Galleries & Tutorials (click on the image for details)

[![](image/ATAC.png)](./example/STRAP_ATAC_infrastructure.md)
[![](image/RNA.png)](./example/STRAP_RNA_infrastructure.md)
[![](image/INTERGRATE.png)](./example/STRAP_INTERGRATE_infrastructure.md)


## 3. FAQs
[Q1:]()             
[Q2:]()          
[Q3:]()   


## 4. Citation


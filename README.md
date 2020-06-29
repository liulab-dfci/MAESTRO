# MAESTRO

![GitHub](https://img.shields.io/github/license/liulab-dfci/MAESTRO)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/liulab-dfci/MAESTRO)
[![Conda](https://img.shields.io/conda/dn/liulab-dfci/maestro?label=Conda%20downloads)](https://anaconda.org/liulab-dfci/maestro)
[![Docker Pulls](https://img.shields.io/docker/pulls/winterdongqing/maestro)](https://hub.docker.com/repository/docker/winterdongqing/maestro)

**MAESTRO**(**M**odel-based **A**nalys**E**s of **S**ingle-cell **T**ranscriptome and **R**egul**O**me) is a comprehensive single-cell RNA-seq and ATAC-seq analysis suit built using [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home). MAESTRO combines several dozen tools and packages to create an integrative pipeline, which enables scRNA-seq and scATAC-seq analysis from raw sequencing data (fastq files) all the way through alignment, quality control, cell filtering, normalization, unsupervised clustering, differential expression and peak calling, celltype annotation and transcription regulation analysis. Currently, MAESTRO support [Smart-seq2](https://www.ncbi.nlm.nih.gov/pubmed/24385147), [10x-genomics](https://www.10xgenomics.com/solutions/single-cell/), [Drop-seq](https://www.cell.com/abstract/S0092-8674(15)00549-8), [SPLiT-seq](https://science.sciencemag.org/content/360/6385/176) for scRNA-seq protocols; [microfudics-based](https://www.ncbi.nlm.nih.gov/pubmed/26083756), [10x-genomics](https://www.10xgenomics.com/solutions/single-cell-atac/) and [sci-ATAC-seq](https://science.sciencemag.org/content/348/6237/910) for scATAC-seq protocols.       
        
## Change Log

### v1.0.0
* Release MAESTRO.
### v1.0.1
* Provide [docker image](https://hub.docker.com/repository/docker/winterdongqing/maestro) for easy installation. Note, the docker does not include cellranger/cellranger ATAC, as well as the corresponding genome index. Please install cellranger/cellranger ATAC following the installation instructions.
### v1.0.2
* Fix some bugs and set LISA as the default method to predict transcription factors for scRNA-seq. Note, the docker includes the LISA conda environment, but does not include required pre-computed genome datasets. Please download hg38 or mm10 datasets and update the configuration following the installation instructions.
### v1.1.0
* Change the default alignment method of MAESTRO from cellranger to [starsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) and [minimap2](https://github.com/lh3/minimap2) for accelerating the mapping time.
* Improve the memory efficiency of scATAC gene score calculation.
* Incorporate the installation of giggle into MAESTRO, and add web API for LISA function. All the core MAESTRO function can be installed through the conda environment now!
* Provide more documents for the QC parameters and add flexibility for other parameters in the workflow.
### v1.2.0
* Modify the regulatory potential model by removing the interfering peaks from adjacent genes and adjusting the weight of exon peaks. The "adjusted RP model" is set as the default gene activity scoring model with original "simple RP model" as a option.
* Modify the integration function of MAESTRO. The new function can output more intermediate figures and log files for diagnosing the possible failure in integrating rare populations.
* Add the function for annotating cell-types for scATAC-seq clusters based on public bulk chromatin accessibility data from Cistrome database.
* Provide the function of generating genome browser tracks at cluster level for scATAC-seq dataset visualization. 
* Support peak calling at the cluster level now!

## System requirements
* Linux/Unix
* Python (>= 3.0) for MAESTRO snakemake workflow
* R (>= 3.6.1) for MAESTRO R package

## Installation

### Install MAESTRO

There are two ways to install MAESTRO -- to install the full workflow through [Anaconda cloud](https://www.anaconda.com/); or to install only the R codes for exploring the processed data.

#### Installing the full solution of MAESTRO workflow through conda

MAESTRO uses the [Miniconda3](http://conda.pydata.org/miniconda.html) package management system to harmonize all of the software packages. Users can install the full solution of MAESTRO using the conda environment.

Use the following commands to install Minicoda3：
``` bash
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```
And then users can create an isolated environment for MAESTRO and install through the following commands:
``` bash
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
# To make the installation faster, we recommend using mamba
$ conda install mamba -c conda-forge
$ mamba create -n MAESTRO maestro=1.2.0 -c liulab-dfci
# Activate the environment 
$ conda activate MAESTRO
```

#### Installing the MAESTRO R package from source code

If users already have the processed datasets, like cell by gene or cell by peak matrix generate by Cell Ranger. Users can install the MAESTRO R package alone to perform the analysis from processed datasets.
``` bash
$ R
> library(devtools)
> install_github("liulab-dfci/MAESTRO")
```

### Required annotations for MAESTRO workflow

The full MAESTRO workflow requires extra annotation files and tools:

* MAESTRO depends on [starsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) and [minimap2](https://github.com/lh3/minimap2) for mapping scRNA-seq and scATAC-seq dataset. Users need to generate **the reference files** for the alignment software and specify the path of the annotations to MAESTRO through command line options.

* MAESTRO utilizes **LISA** to evaluate the enrichment of transcription factors based on the marker genes from scRNA-seq clusters. By default, users can choose the "web" option, which will use the API function in MAESTRO to perform LISA analysis. However, if users want to use the local version of LISA, they need to install [LISA](https://github.com/qinqian/lisa) locally, build the annotation files according to the LISA document, and provide the path of LISA to MAESTRO when using the RNAAnnotateTranscriptionFactor function.

* MAESTRO utilizes **giggle** to identify enrichment of transcription factor peaks in scATAC-seq cluster-specific peaks. By default [giggle](https://github.com/ryanlayer/giggle) is installed in MAESTRO environment. The giggle index for Cistrome database can be downloaded [here](http://cistrome.org/~chenfei/MAESTRO/giggle.all.tar.gz). Users need to download the file and provide the location of the giggle annotation to MAESTRO when using the ATACAnnotateTranscriptionFactor function.

## Usage
```
usage: MAESTRO [-h] [-v]
               {scrna-init,scatac-init,integrate-init,mtx-to-h5,count-to-h5,merge-h5,scrna-qc,scatac-qc,scatac-peakcount,scatac-genescore}
```

There are ten functions available in MAESTRO serving as sub-commands.

Subcommand | Description
-----------|----------
`scrna-init` | Initialize the MAESTRO scRNA-seq workflow.
`scatac-init` | Initialize the MAESTRO scATAC-seq workflow.
`integrate-init` | Initialize the MAESTRO integration workflow.
`mtx-to-h5` | Convert 10X mtx format matrix to HDF5 format.
`count-to-h5` | Convert plain text count table to HDF5 format.
`merge-h5` | Merge multiple HDF5 files, e.g. different replicates.
`scrna-qc` | Perform quality control for scRNA-seq gene-cell count matrix.
`scatac-qc` | Perform quality control for scATAC-seq peak-cell count matrix.
`scatac-peakcount` | Generate peak-cell binary count matrix.
`scatac-genescore` | Calculate gene score based on the binarized scATAC peak count.

Example of running MAESTRO can be found at the following galleries. Please use `MAESTRO
COMMAND -h` to see the detail description for each option of each module.

## Galleries & Tutorials (click on the image for details)

[<img src="./image/RNA.10x.png" width="297" height="378" />](./example/RNA_infrastructure_10x/RNA_infrastructure_10x.md)
[<img src="./image/ATAC.10x.png" width="297" height="378" />](./example/ATAC_infrastructure_10x/ATAC_infrastructure_10x.md)  
[<img src="./image/integration.10x.png" width="297" height="378" />](./example/Integration/Integration.md)
[<img src="./image/integration.large.png" width="297" height="378" />](./example/Integration_large/Integration_large.md)

## Citation


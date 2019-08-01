# STRAP RNA infrastructure

10X PBMCs 5K
In this example, we will be analyzing a dataset of 5K Peripheral blood mononuclear cells (PBMCs) of a healthy donor freely available from 10X.

## Table of Contents
[Step 0. Activate the STRAP environment](#system-requirements)             
[Step 1. Prepare your working directory](#annotation)        
[Step 2: Configure the workflow](#)      
[Step 3. Running STRAP](#SettingUpForProject)     
[Step 4. Output Files](#Output Files)              

### **Step 1. Prepare your working directory**

All work in STRAP is done in a PROJECT directory, which is simply a directory to contain a single STRAP analysis run. PROJECT directories can be named anything. Here we name it "10x.pbmc.5k".

You can initialize the workflow with
```
strap init -d 10X_PBMC_5k -m scRNA
```      

The raw data can be downloaded from 10X genomics:
```
cd 10x.pbmc.5k
mkdir data
cd data
wget http:http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
tar xvf 5k_pbmc_v3_fastqs.tar
```               

Now we have a workflow directory, and a set of FASTQ files for analysis. 

Now the PROJECT directory(`10x.pbmc.5k`) is that you fill them with the following core components: (We first lay out the directory structure and explain each element below) 
>10x.pbmc.5k/
>>strap/    
>>data/        
>>ref_files/     
>>config.yaml                       

The `strap` directory contains all of the STRAP code.The `data` directory contains all of your raw data.The `ref_files/` folder includes the genome annotation files.The `config.yaml` and metasheet.csv are configurations for your STRAP run (explained further in next section). After a successful STRAP run, another 'analysis' folder is generated which contains all of the resulting output files. 

### Step 2. Configure the workflow                             

Open the `config.yaml` file and edit it to your needs. Especially, define your single-cell platform for use. 

```
cd /root/strap/Snakemake
vi config.yaml
```

Here is an example for `config.yaml` file.

```
# Directory where fastq files are stored
fastqdir: /root/strap/Data/pbmc_k_v2_fastqs
# Sample name of fastq file (only for platform of "10xGenomics", for example, 
# if there is a file named pbmc_1k_v2_S1_L001_I1_001.fastq.gz, the sample name is "pbmc_1k_v2". )
fastqprefix: pbmc_1k_v2
# Species to use [GRCh38, mmu] (GRCh38 for human and mmu for mouse)
species: GRCh38
# Method to use [Seurat, Pagoda, scMCA, RCA, SSCC]
method: Seurat
# Platform of single cell RNA-seq [Smartseq2, 10xGenomics, Dropseq]
platform: 10xGenomics
# The prefix of output files
outprefix: pbmc_1k_v2
# Number of cores to use
cores: 8

# Reference genome 
genome:
  # # Genome index directory for STAR
  # mapindex: /mnt/Storage/home/sundongqing/RefGenome/hg38/STAR_index
  # .gtf format genome annotation file
  gtf: /root/strap/RefGenome/hg38/Homo_sapiens.GRCh38.92.gtf
  # .bed format genome annotation file 
  bed: /root/strap/RefGenome/hg38/hg38_gencode.v28.bed
  # # .txt format genome annotation file (only for platform of "Dropseq")
  # anno: /mnt/Storage/home/sundongqing/RefGenome/hg38/hg38_gencode_annotation.V28.txt
  # genome annotation file from 10xGenomics required for Cell Ranger
  cellranger: /root/strap/RefGenome/hg38/refdata-cellranger-GRCh38-3.0.0
```

If you start from count table of 10xGenomics, you need to move the data file to the Result folder and all the file and folder names except `pbmc_1k_v2` must be the same as the example here.

```
knitr::include_graphics("/Users/dongqing/Documents/Project/SingleCell/scRNA/scr/Figure/folder_structure2.png")
```

If you start from count table of other platforms (e.g., Dropseq, Smartseq2 and so on), you need to organize you file as the following. For example, if you have a count table `GSE98638_HCC.TCell.S5063.count.txt`, you need to rename it `GSE98638_HCC.TCell.S5063_expmat.txt` and move it to the Result folder. And then you need to define the platform as 'Smartseq2' (even if the platform is not 'Smartseq2') and the outprefix as "GSE98638_HCC.TCell.S5063".

```
knitr::include_graphics("/Users/dongqing/Documents/Project/SingleCell/scRNA/scr/Figure/folder_structure3.png")
```           

Once configured, the workflow can be executed with Snakemake.                

### **Step 3. Running STRAP**

To start, we must activate the STRAP CONDA ENVIRONMENT.If successful, you will see "(strap)" prepended to your command prompt.

Next we will perform a DRY-RUN to make sure that we setup the STRAP PROJECT directory correctly. In your PROJECT folder run the following command:
```
nohup snakemake --cores 8 --use-conda >10x.pbmc.1k.out &
```

### **Step 4. Output Files**

Here, we assume you've run STRAP successfully. An output directory is specified in the run() call, and will contain several useful outputs as described below.            



which will display all jobs that will be executed. If there are errors in your config file, you will be notified by Snakemake. The actual execution of the workflow can be started with

```
snakemake
```

The output files are as follows.


Here we show an annotated tSNE plot.



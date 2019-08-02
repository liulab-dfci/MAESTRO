# STRAP ATAC infrastructure

## 10X PBMCs 5K

In this example, we will be analyzing a dataset of 5K Peripheral blood mononuclear cells (PBMCs) of a healthy donor freely available from 10X. 


## Table of Contents
[Step 0. Activate the STRAP environment](#system-requirements)               
[Step 1. Prepare your working directory](#annotation)        
[Step 2: Configure the workflow](#)      
[Step 3. Running STRAP](#SettingUpForProject)     
[Step 4. Output Files](#Output Files)



### **Step 0. Activate the STRAP environment**
If you have installed **STRAP** as shown above, make sure to activate the corresponding conda environment before conducting any mageck-vispr related command via

```
source activate strap

```

### **Step 1. Prepare your working directory**

All work in STRAP is done in a PROJECT directory, which is simply a directory to contain a single STRAP analysis run. PROJECT directories can be named anything. Here we name it "10x.pbmc.5k".

You can initialize the workflow with
```
strap init 10x.pbmc.5k -m scATAC
```

The raw data can be downloaded from 10X genomics:
```
cd 10x.pbmc.5k
mkdir data
cd data
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_5k_v1/atac_pbmc_5k_v1_fastqs.tar
```

Now we have a workflow directory, and a set of FASTQ files for analysis. 

Now the PROJECT directory(`10x.pbmc.5k`) is that you fill them with the following core components: (We first lay out the directory structure and explain each element below) 
>10x.pbmc.5k/
>>strap/    
>>data/        
>>ref_files/     
>>config.yaml

The `strap` directory contains all of the STRAP code.The `data` directory contains all of your raw data.The `ref_files/` folder includes the genome annotation files.The `config.yaml` and metasheet.csv are configurations for your STRAP run (explained further in next section). After a successful STRAP run, another 'analysis' folder is generated which contains all of the resulting output files. 

### **Step 2. Configure the workflow**

Now you are in the STRAP container and all your files have been added into the container.Open the `config.yaml` file and edit it to your needs. The `config.yaml` file has several main sections : PATH,REFERENCE,PLATFORM.

```
# Directory where fastq files are stored
fastqdir: /root/strap/Data/pbmc_5k_v2_fastqs
# Sample name of fastq file (only for platform of "10xGenomics", for example, 
# if there is a file named pbmc_5k_v2_S1_L001_I1_001.fastq.gz, the sample name is "pbmc_5k_v2". )
fastqprefix: pbmc_5k_v2
# Species to use [GRCh38, mmu] (GRCh38 for human and mmu for mouse)
species: GRCh38
# Method to use [Seurat, Pagoda, scMCA, RCA, SSCC]
method: Seurat
# Platform of single cell RNA-seq [Smartseq2, 10xGenomics, Dropseq]
platform: 10xGenomics
# The prefix of output files
outprefix: pbmc_5k_v2
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

Now that we have setup our PROJECT directory (downloading the 'strap' code directory, creating our 'data' directory, and configuring our 'config.yaml'). Finally we are ready to run STRAP.     

### **Step 3. Running STRAP**

To start, we must activate the STRAP CONDA ENVIRONMENT.If successful, you will see "(strap)" prepended to your command prompt.

Next we will perform a DRY-RUN to make sure that we setup the STRAP PROJECT directory correctly. In your PROJECT folder run the following command:
```
nohup snakemake --cores 8 --use-conda >10x.pbmc.1k.out &
```


### **Step 4. Output Files**

Here, we assume you've run STRAP successfully. An output directory is specified in the run() call, and will contain several useful outputs as described below.            

#### 1. Alignment and Quality Control               
**10X_PBMC_frag.pdf**: The fragment size density of your libraries at bulk level.In a normal ATAC-seq library, you should expect to see a sharp peak at the <100 bp region (open chromatin), and a peak at ~200bp region (mono-nucleosome), and other larger peaks (multi-nucleosomes).            
![](/example/STRAP_ATAC_infrastructure/10X_PBMC_5k_frag.png)                               
 
**10X_PBMC_map.pdf**: This file include the frequency of total fragments, mapping raito, the frequency of duplicate aligned reads, and the frequency of mitochondrial aligned reads.
![](/example/STRAP_ATAC_infrastructure/10X_PBMC_5k_map.png)               

**10X_PBMC_frip.pdf**: The fraction of promoter reads in each cells(including and filtered samples).     
![](/example/STRAP_ATAC_infrastructure/10X_PBMC_5k_frip.png)               

**10X_PBMC_barcodes.txt**: the barcodes of each cell.       

#### 2. Clustering                    

**10X_PBMC_strap_peak_umap_cluster.pdf**: The UMAP plot of clustering.    
![](/example/STRAP_ATAC_infrastructure/10X_PBMC_5k_UMAP_cluster.png)              

**10X_PBMC_strap_peak_SeuratObj.rds**:The peak file.     

**10X_PBMC_strap_peak_ClusterSpecificPeaks.rds**:The specific peak file of each cluster.     


#### 3. Annotating   

**UMAP_assignIdent_10X_PBMC_strap_peak_wilcox_annotated.png**：Celltype annotation based on differential expression genes.  
![](/example/STRAP_ATAC_infrastructure/UMAP_assignIdent_10X_PBMC_5k_RPannotated.png)  

**10X_PBMC_strap_peak_ClusterSpecificGenes_wilcox.rds**: The file of specific expressed genes of each cluster.

# MAESTRO v1.0.2 docker

If it is difficult to install MAESTRO because of the system configuration, users can adopt an alternative way to get MAESTRO running through a docker container.

## Get Started

### Install docker

To start a docker container, please install [docker](https://docs.docker.com/install/) on your platform first. 

**Note:** On OS X, Docker Machine has Memory and CPU limit. We recommend to install it on linux.

### Pull the image

The docker distribution includes the latest MAESTRO conda package, as well as RABIT, GIGGLE and corresponding genome index. The [LISA](https://github.com/qinqian/lisa) conda environment has been created in the docker image, but users still need to download required hg38 or mm10 datasets and update the configuration according to the instructions.

**Note:** The docker image does not include cellranger/cellrangerATAC and required reference genome. Please install cellranger/cellranger ATAC locally following the installation instructions.

Please use the following commands to the install Minicoda3:

``` bash
$ sudo docker pull winterdongqing/maestro:1.0.2
$ sudo docker image list
REPOSITORY               TAG                 IMAGE ID            CREATED             SIZE
winterdongqing/maestro   1.0.2               d7012fe3925d        3 hours ago         11.4GB
```

### Run in the container

To run the MAESTRO container, please use the following commands. The ```-v``` flag mounts the current directory ```/home1/wangchenfei``` into ```/root/MAESTRO``` in the container.

```bash
$ pwd
/home1/wangchenfei
$ ls
annotations  docker  miniconda3  ncbi  Project  R  tmp  Tool
$ sudo docker run -v /home1/wangchenfei:/root/MAESTRO -it winterdongqing/maestro:1.0.2
```

```#``` means that users have been in the container. There are four folders in the  ```root``` directory of the container. The ```Annotation``` directory stores the RABIT index and GIGGLE index. Users can access local data and locally installed softwares through the ```MAESTRO```.

```bash
root@901b0b615e2e:~# ls
Annotation  MAESTRO  Software  miniconda3
root@901b0b615e2e:~# conda env list
# conda environments:
#
base                  *  /root/miniconda3
lisa                     /root/miniconda3/envs/lisa

root@901b0b615e2e:~# cd MAESTRO/ && ls
Project  R  Tool  annotations  docker  miniconda3  ncbi  tmp
```

Before running the MAESTRO pipeline, please prepend the Cell Ranger directory to the ```$PATH``` in the container. This will allow you to invoke the ```cellranger``` command.

**Note:** The ```cellranger-3.1.0``` in the local directory ```/home1/wangchenfei/Tool/``` is now in the ```/root/MAESTRO/Tool/``` of the container.

```bash
root@901b0b615e2e:~/MAESTRO# export PATH=/root/MAESTRO/Tool/cellranger-3.1.0:$PATH
```

Then users can run MAESTRO following the [tutorial](../example/RNA_infrastructure_10x/RNA_infrastructure_10x.md).
```bash
root@901b0b615e2e:~/MAESTRO# cd ~/MAESTRO/Project/SingleCell/scRNA/Analysis/
root@901b0b615e2e:~/MAESTRO/Project/SingleCell/scRNA/Analysis# MAESTRO init -m scRNA -d 10X_PBMC_1k_docker
root@901b0b615e2e:~/MAESTRO/Project/SingleCell/scRNA/Analysis# cd 10X_PBMC_1k_docker
root@901b0b615e2e:~/MAESTRO/Project/SingleCell/scRNA/Analysis# ls
Snakefile  config.yaml
```

Here is an example of the ```config.yaml``` file needed to run MAESTRO scRNA module in the container.
```bash
root@901b0b615e2e:~/MAESTRO/Project/SingleCell/scRNA/Analysis/10X_PBMC_1k_docker# vi config.yaml
# Directory where fastq files are stored
fastqdir: /root/MAESTRO/Project/SingleCell/scRNA/Analysis/10xPBMC_1k/Data/pbmc_1k_v3_fastqs
# Sample name of fastq file (only for platform of "10x-genomics", for example,
# if there is a file named pbmc_1k_v2_S1_L001_I1_001.fastq.gz, the sample name is "pbmc_1k_v2". )
fastqprefix: pbmc_1k_v3
# Species to use [GRCh38, GRCm38] (GRCh38 for human and GRCm38 for mouse)
species: GRCh38
# Platform of single cell RNA-seq [Smartseq2, 10x-genomics, Dropseq]
platform: 10x-genomics
# The prefix of output files
outprefix: pbmc_1k_v3
# Run RSeQC or not [True, False]
rseqc: False
# Number of cores to use
cores: 8
# annotation to run rabit
rabitlib: /root/MAESTRO/Project/SingleCell/scATAC/Code/MAESTRO/MAESTRO/annotations/Rabit_lib

# Reference genome
genome:
  # Genome index directory for STAR
  mapindex: /root/MAESTRO/annotations/refdata-cellranger-GRCh38-3.0.0/star
  # .gtf format genome annotation file
  gtf: /root/MAESTRO/annotations/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
  # genome annotation file from 10x-genomics required for Cell Ranger
  cellranger: /root/MAESTRO/annotations/refdata-cellranger-GRCh38-3.0.0
  # the prefix of transcript references for RSEM used by rsem-prepare-reference
  rsem: /root/MAESTRO/annotations/hg38/RSEM_ref/GRCh38
```

Once configured, users can use snakemake to run the workflow as they do in their local system.

```bash
root@901b0b615e2e:~/MAESTRO/Project/SingleCell/scRNA/Analysis/10X_PBMC_1k_docker# snakemake -np
Building DAG of jobs...
Job counts:
	count	jobs
	1	all
	1	scrna_analysis
	1	scrna_cellranger
	1	scrna_qc
	1	scrna_report
	1	scrna_rseqc_plot
	6

...
...
...

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
root@901b0b615e2e:~/MAESTRO/Project/SingleCell/scRNA/Analysis/10X_PBMC_1k_docker# nohup snakemake --cores 10 > 10X_PBMC_1k_docker.out &
```

## Built with
* Ubuntu 18.04 docker image
* miniconda3
* R 3.5.1

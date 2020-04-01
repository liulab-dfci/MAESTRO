# MAESTRO v1.1.0 docker

If it is difficult to install MAESTRO because of the system configuration, users can adopt an alternative way to get MAESTRO running through a docker container.

## Get Started

### Install docker

To start a docker container, please install [docker](https://docs.docker.com/install/) on your platform first. 

**Note:** On OS X, Docker Machine has Memory and CPU limit. We recommend to install it on linux.

### Pull the image

The docker distribution includes the latest MAESTRO conda package, as well as RABIT, GIGGLE and corresponding genome index. The [LISA](https://github.com/qinqian/lisa) conda environment has been created in the docker image, but users still need to download required hg38 or mm10 datasets and update the configuration according to the instructions.

**Note:** The docker image does not include required reference genome. Please download them following the help description of `MAESTRO scatac-init` or `MAESTRO scrna-init`.

Please use the following commands to the pull the MAESTRO docker image.

``` bash
$ sudo docker pull winterdongqing/maestro:1.1.0
$ sudo docker image list
REPOSITORY               TAG                 IMAGE ID            CREATED             SIZE
winterdongqing/maestro   1.1.0               831bfd245e63        3 hours ago         8.7GB
```

### Run in the container

To run the MAESTRO container, please use the following commands. The ```-v``` flag mounts the current directory ```/home1/wangchenfei``` into ```/root/MAESTRO``` in the container.

```bash
$ pwd
/home1/wangchenfei
$ ls
annotations  docker  miniconda3  ncbi  Project  R  tmp  Tool
$ sudo docker run -v /home1/wangchenfei:/root/MAESTRO -it winterdongqing/maestro:1.1.0
```

```#``` means that users have been in the container. There are four folders in the  ```root``` directory of the container. The ```Annotation``` directory stores the RABIT index and GIGGLE index. Users can access local data and locally installed softwares through the ```MAESTRO```.

```bash
root@45c5df8b6635:~# ls
MAESTRO  miniconda3
root@45c5df8b6635:~# conda env list
# conda environments:
#
base                  *  /root/miniconda3
lisa                     /root/miniconda3/envs/lisa

root@45c5df8b6635:~# cd MAESTRO/ && ls
Project  R  Tool  annotations  docker  miniconda3  ncbi  tmp
```

Then users can run MAESTRO following the [tutorial](../example/ATAC_infrastructure_10x/ATAC_infrastructure_10x.md).
```bash
root@45c5df8b6635:~/MAESTRO# cd ~/MAESTRO/Project/SingleCell/scATAC/Analysis/
root@45c5df8b6635:~/MAESTRO/Project/SingleCell/scATAC/Analysis# MAESTRO scatac-init -d 10X_PBMC_1k_MAESTRO_V110_docker \
--fastq-dir /root/MAESTRO/Project/SingleCell/scATAC/Analysis/10X_PBMC_1k/Data/atac_pbmc_1k_v1_fastqs --fastq-prefix atac_pbmc_1k_v1 \
--outprefix 10X_PBMC_1k_MAESTRO_V110 \
--whitelist /root/MAESTRO/Tool/cellranger-atac-1.1.0/cellranger-atac-cs/1.1.0/lib/python/barcodes/737K-cratac-v1.txt \
--giggleannotation /root/MAESTRO/Project/SingleCell/scATAC/Code/MAESTRO_V110/MAESTRO/annotations/giggle \
--fasta /root/MAESTRO/annotations/MAESTRO/Refdata_scATAC_MAESTRO_GRCh38_1.1.0/GRCh38_genome.fa
root@45c5df8b6635:~/MAESTRO/Project/SingleCell/scATAC/Analysis# cd 10X_PBMC_1k_MAESTRO_V110_docker
root@45c5df8b6635:~/MAESTRO/Project/SingleCell/scATAC/Analysis# ls
Snakefile  config.yaml
```

Once configured, users can use snakemake to run the workflow as they do in their local system.

```bash
root@45c5df8b6635:~/MAESTRO/Project/SingleCell/scATAC/Analysis/10X_PBMC_1k_MAESTRO_V110_docker# snakemake -np
Building DAG of jobs...
Job counts:
	count	jobs
	1	all
	1	scatac_allpeakcall
	1	scatac_analysis
	1	scatac_barcodecorrect
	1	scatac_countpeak
	1	scatac_fragmentcorrect
	1	scatac_fragmentgenerate
	1	scatac_genescore
	1	scatac_map
	1	scatac_mergepeak
	1	scatac_qcfilter
	1	scatac_qcplot
	1	scatac_qcstat_bulk
	1	scatac_qcstat_mapped
	1	scatac_qcstat_promoter
	1	scatac_qcstat_singlecell
	1	scatac_report
	1	scatac_rmdp
	18

...
...
...

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
root@45c5df8b6635:~/MAESTRO/Project/SingleCell/scATAC/Analysis/10X_PBMC_1k_MAESTRO_V110_docker# nohup snakemake --cores 10 > 10X_PBMC_1k_MAESTRO_V110_docker.out &
```

## Built with
* Ubuntu 18.04 docker image
* miniconda3
* R 3.6.2

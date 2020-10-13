# Snakemake pipeline for processing multi-sample scATACseq data


### Initiate a folder containing Snakefile and config.ymal file for multi-sample processing

Install MAESTRO following the [instruction](https://github.com/liulab-dfci/MAESTRO#installation) on the main page.


```bash
conda activate MAESTRO

MAESTRO mulit-scatac-init

```


### Create a Json file containing the fastq path for all samples

```bash
makeJsonFromFastq scatac /path/to/fastq/dir/
```

A `samples.json` will be generated. Inspect it to make sure it looks correct.


### Run the workflow

```bash
# dry run
snakemake -j 8 -np
```
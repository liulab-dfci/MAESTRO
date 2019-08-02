## 10x-genomics based scRNA-seq from human PBMC samples

In this example, we will be analyzing the scRNA-seq dataset of 8K human peripheral blood mononuclear cells (PBMCs) freely available from 10X. The raw data can be downloaded from [here](http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_fastqs.tar). 

**Step 0. Download the data and prepare the working directory**

You can initialize the workflow with
```bash
$ 
```      

The raw data can be downloaded from 10X genomics:
``` bash
$ cd 10x.pbmc.5k
$ mkdir data
$ cd data
$ wget http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_fastqs.tar
$ tar xvf pbmc8k_fastqs.tar
```               

**Step 1. Configure the workflow**                           

Open the `config.yaml` file and edit it to your needs. Especially, define your single-cell platform for use. 

```bash
```

Here is an example for `config.yaml` file.

```bash
```      

Once configured, the workflow can be executed with Snakemake.                

**Step 2. Run MAESTRO**

To start, we must activate the STRAP CONDA ENVIRONMENT.If successful, you will see "(MAESTRO)" prepended to your command prompt.

Next we will perform a DRY-RUN to make sure that we setup the STRAP PROJECT directory correctly. In your PROJECT folder run the following command:

```bash
```

**Step 3. Visualize the output files**

Here, we assume you've run STRAP successfully. An output directory is specified in the run() call, and will contain several useful outputs as described below.            

```bash
```

The output files are as follows.


**Step 4. Understand the quality control results**


```bash
```

**Step 5. Understand the analysis results**

```bash
```



















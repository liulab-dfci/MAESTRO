argue = commandArgs(T)
cwd = getwd()
setwd(argue[7])
source("scRNAseq_function.R")
setwd(file.path(cwd,argue[6]))

if(argue[3] == "10xGenomics"){
  count_file = paste0(argue[1],"matrix.mtx")
  barcode_file = paste0(argue[1],"barcodes.tsv")
  gene_file = paste0(argue[1],"features.tsv")
  exp.dat = readMM(count_file)
  gene_list = read.table(gene_file , sep = '\t', stringsAsFactors = FALSE, header = FALSE)
  barcode_list = read.table(barcode_file, sep = '\t', stringsAsFactors = FALSE, header = FALSE)
  colnames(exp.dat) = barcode_list$V1
  row.names(exp.dat) = gene_list$V1
  exp.dat = as.matrix(exp.dat)
  if(is.na(argue[9])){
    idtype = "ENSEMBL"
  } else{
    idtype = argue[9]
  }
} else if(argue[3] == "Smartseq2"){
    if(is.na(argue[9])){
      idtype = "ENSEMBL"
    } else{
      idtype = argue[9]
    }
    exp.dat <- read.table(argue[1], header = TRUE, row.names = 1, sep = '\t')
    exp.dat = exp.dat[!duplicated(row.names(exp.dat)),]
} else{
    if(is.na(argue[9])){
      idtype = "SYMBOL"
    } else{
      idtype = argue[9]
    }
    exp.dat <- read.table(argue[1], header = TRUE, row.names = 1, sep = '\t')
    exp.dat = exp.dat[!duplicated(row.names(exp.dat)),]
}

if(argue[4] == "RCA"){
  exp.dat <- Count2FPKM(exp.dat, idType = idtype, organism = argue[2])
  PipelineRCA(exp.dat, proj = argue[5])
} else{
  exp.dat <- Count2TPM(exp.dat, idType = idtype, organism = argue[2])
  exp.dat <- log2(exp.dat/10+1)
  if(argue[4] == "Seurat"){
    Seurat <- PipelineSeurat(exp.dat, proj = argue[5], min.c = 10, min.g = 1000)
  }
  if(argue[4] == "Pagoda"){
    Pagoda <- PipelinePagoda(exp.dat, proj = argue[5])
  }
  if(argue[4] == "SSCC"){
    PipelineSSCC(exp.dat, proj = argue[5])
  }
  if(argue[4] == "scMCA"){
    PipelinescMCA(exp.dat, proj = argue[5])
  }
}

rabitlibdir = argue[8]
RunRabit(Seurat$markers,Seurat$SeuratObj, rabitlibdir = rabitlibdir)



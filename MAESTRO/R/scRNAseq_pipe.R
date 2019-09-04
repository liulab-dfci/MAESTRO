library(MAESTRO)
library(Seurat)
argue = commandArgs(T)

setwd(argue[6])

if(argue[3] == "10x-genomics"){
  exp.dat = Read10X(data.dir = argue[1], gene.column = 2, unique.features = TRUE)
  exp.dat = as.matrix(exp.dat)
} else if(argue[3] == "Smartseq2"){
  exp.dat <- read.table(argue[1], header = TRUE, row.names = NULL, sep = '\t')
  exp.dat = exp.dat[!duplicated(exp.dat[,1]),]
  rownames(exp.dat) = exp.dat[,1]
  exp.dat = exp.dat[,-1]
  # exp.dat <- RNACountToTPM(exp.dat, idType = "Ensembl", organism = argue[2])
} else{
  exp.dat <- read.table(argue[1], header = TRUE, row.names = 1, sep = '\t')
  exp.dat = exp.dat[!duplicated(exp.dat[1,]),]
  # exp.dat <- RNACountToTPM(exp.dat, idType = "Symbol", organism = argue[2])
}

RNA.res <- RNARunSeurat(inputMat = exp.dat, project = argue[4], min.c = 10, min.g = 100)
RNA.res$RNA <- RNAAnnotateCelltype(RNA = RNA.res$RNA, genes = RNA.res$genes, signatures = human.immune.CIBERSORT, min.score = 0.05)
saveRDS(RNA.res$RNA, paste0(argue[4], "_scRNA_Object.rds"))
RNA.tfs <- RNAAnnotateTranscriptionFactor(RNA = RNA.res$RNA, genes = RNA.res$genes, project = argue[4], rabit.path = argue[5], organism = argue[2], top.tf = 10)

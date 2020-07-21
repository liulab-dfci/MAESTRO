library(MAESTRO)
library(Seurat)
library(future)
library(future.apply)
library(pbapply)
library(optparse)


option_list = list(
  make_option(c("--prefix"), type = "character", default = "MAESTRO",
              action = "store", help = "The prefix of the output files."
  ),
  make_option(c("--expression"), type = "character", default = "",
              action = "store", help = "The count matrix h5 file."
  ),
  make_option(c("--outdir"), type = "character", default = "Result/Analysis",
              action = "store", help = "The directory where the output files are stored."
  ),
  make_option(c("--species"), type = "character", default = "GRCh38",
              action = "store", help = "The platform of scRNA-seq."
  ),
  make_option(c("--method"), type = "character", default = "LISA",
              action = "store", help = "The method to identify driver regulators. [LISA, RABIT]"
  ),
  make_option(c("--signature"), type = "character", default = "",
              action = "store", help = "The cell signature file for celltype annotation. Default is built-in CIBERSORT immune cell signature."
  ),
  make_option(c("--lisamode"), type = "character", default = "",
              action = "store", help = "Mode to run LISA (web or local)."
  ),
  make_option(c("--condadir"), type = "character", default = "",
              action = "store", help = "Directory where miniconda or anaconda is installed (only if method is set to lisa)."
  ),
  make_option(c("--lisaenv"), type = "character", default = "lisa",
              action = "store", help = "Name of lisa environment (only if method is set to lisa)."
  ),
  make_option(c("--thread"), type = "integer", default = 1,
              action = "store", help = "Number of cores to use."
  )
)

argue = parse_args(OptionParser(option_list = option_list, usage = "Run scRNA-seq analysis pipeline"))

setwd(argue$outdir)
count_mat = argue$expression
prefix = argue$prefix
thread = argue$thread
method = argue$method
sigfile = argue$signature
lisamode = argue$lisamode
condadir = argue$condadir
lisaenv = argue$lisaenv
species = argue$species


# if(argue[3] == "10x-genomics"){
#   exp.dat = Read10X(data.dir = argue[1], gene.column = 2, unique.features = TRUE)
#   exp.dat = as.matrix(exp.dat)
# } else if(argue[3] == "Smartseq2"){
#   exp.dat <- read.table(argue[1], header = TRUE, row.names = NULL, sep = '\t')
#   exp.dat = exp.dat[!duplicated(exp.dat[,1]),]
#   rownames(exp.dat) = exp.dat[,1]
#   exp.dat = exp.dat[,-1]
#   # exp.dat <- RNACountToTPM(exp.dat, idType = "Ensembl", organism = argue[2])
# } else{
#   exp.dat <- read.table(argue[1], header = TRUE, row.names = 1, sep = '\t')
#   exp.dat = exp.dat[!duplicated(exp.dat[1,]),]
#   # exp.dat <- RNACountToTPM(exp.dat, idType = "Symbol", organism = argue[2])
# }

exp.dat = Read10X_h5(filename = count_mat, unique.features = TRUE)

plan("multiprocess", workers = as.integer(thread))
options(future.globals.maxSize = 10*1024^3)

if(sigfile %in% c("human.immune.CIBERSORT", "mouce.brain.ALLEN", "mouse.all.facs.TabulaMuris", "mouse.all.droplet.TabulaMuris")){
  signatures = sigfile
}else{
  signatures = read.table(sigfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
}
RNA.res <- RNARunSeurat(inputMat = exp.dat, project = prefix, min.c = 10, min.g = 100)
RNA.res$RNA <- RNAAnnotateCelltype(RNA = RNA.res$RNA, genes = RNA.res$genes, 
                                   signatures = signatures, min.score = 0.05)
saveRDS(RNA.res, paste0(prefix, "_scRNA_Object.rds"))
RNA.tfs <- RNAAnnotateTranscriptionFactor(RNA = RNA.res$RNA, genes = RNA.res$genes, project = prefix, 
                                          method = method, lisa.mode = lisamode, 
                                          conda.dir = condadir, lisa.envname = lisaenv, 
                                          organism = species, top.tf = 10)

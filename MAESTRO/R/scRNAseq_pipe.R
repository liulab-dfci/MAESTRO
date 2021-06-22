library(MAESTRO)
library(Seurat)
library(future)
library(future.apply)
library(pbapply)
library(optparse)
library(ggplot2)


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
  #make_option(c("--lisamode"), type = "character", default = "multi",
              #action = "store", help = "Mode to run LISA (multi or one-vs-rest)."
  #),
  #make_option(c("--method"), type = "character", default = "LISA",
              #action = "store", help = "The method to identify driver regulators. [LISA, RABIT]"
  #),
  make_option(c("--signature"), type = "character", default = "",
              action = "store", help = "The cell signature file for celltype annotation. Default is built-in CIBERSORT immune cell signature."
  ),
  make_option(c("--samplefile"), type = "character", default = NULL,
              action="store", help = "For multi-sample, file with the sample origin of each sample"), 
  make_option(c("--thread"), type = "integer", default = 1,
              action = "store", help = "Number of cores to use."
  )
)

argue = parse_args(OptionParser(option_list = option_list, usage = "Run scRNA-seq analysis pipeline"))

setwd(argue$outdir)
count_mat = argue$expression
prefix = argue$prefix
thread = argue$thread
#method = argue$method
sigfile = argue$signature
species = argue$species
samplefile = argue$samplefile


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

if(sigfile %in% c("human.immune.CIBERSORT", "mouse.brain.ALLEN", "mouse.all.facs.TabulaMuris", "mouse.all.droplet.TabulaMuris")){
  signatures = sigfile
}else{
  signatures = read.table(sigfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
}
RNA.res <- RNARunSeurat(inputMat = exp.dat, project = prefix, min.c = 10, min.g = 100)
RNA.res$RNA <- RNAAnnotateCelltype(RNA = RNA.res$RNA, genes = RNA.res$genes,
                                   signatures = signatures, min.score = 0.05)
saveRDS(RNA.res, paste0(prefix, "_scRNA_Object.rds"))
RNA.tfs <- RNAAnnotateTranscriptionFactor(RNA = RNA.res$RNA, genes = RNA.res$genes, project = prefix,
                                          organism = species, top.tf = 10)

if(!is.null(samplefile)) {
    sampleAnno <- scan(samplefile, what=list("character","character"))
    names(sampleAnno[[1]]) <- sampleAnno[[2]]
    RNA.res$RNA$sample <- sampleAnno[[1]][colnames(RNA.res$RNA)]
    p = DimPlot(RNA.res$RNA, label=T, pt.size=0.2, group.by="sample", label.size=3, repel=T)
    ggsave(paste0(RNA.res$RNA@project.name, "_samples.png"), p , width=6, height=4)
}
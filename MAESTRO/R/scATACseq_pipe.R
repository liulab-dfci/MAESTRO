library(MAESTRO)
library(Seurat)
library(Matrix)
library(future)
library(future.apply)
library(pbapply)
library(optparse)


option_list = list(
  make_option(c("--prefix"), type = "character", default = "MAESTRO",
              action = "store", help = "The prefix of the output files."
  ),
  make_option(c("--peakcount"), type = "character", default = "",
              action = "store", help = "The peak count matrix h5 file."
  ),
  make_option(c("--rpmatrix"), type = "character", default = "",
              action = "store", help = "The regulatory potential matrix h5 file."
  ),
  make_option(c("--outdir"), type = "character", default = "Result/Analysis",
              action = "store", help = "The directory where the output files are stored."
  ),
  make_option(c("--species"), type = "character", default = "GRCh38",
              action = "store", help = "The platform of scRNA-seq."
  ),
  make_option(c("--signature"), type = "character", default = "",
              action = "store", help = "The cell signature file for celltype annotation. Default is built-in CIBERSORT immune cell signature."
  ),
  make_option(c("--gigglelib"), type = "character", default = "",
              action = "store", help = "Annotation to run rabit (only if method is set to rabit)."
  ),
  make_option(c("--thread"), type = "integer", default = 1,
              action = "store", help = "Number of cores to use."
  )
)

argue = parse_args(OptionParser(option_list = option_list, usage = "Run scATAC-seq analysis pipeline"))

setwd(argue$outdir)
count_mat = argue$peakcount
rp_mat = argue$rpmatrix
prefix = argue$prefix
thread = argue$thread
sigfile = argue$signature
gigglelib = argue$gigglelib
species = argue$species

# countmatrix = read.table(argue[1], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
# RPmatrix = read.table(argue[2], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
countmatrix = Read10X_h5(count_mat)
RPmatrix = Read10X_h5(rp_mat)

plan("multiprocess", workers = as.integer(thread))
options(future.globals.maxSize = 10*1024^3)

if(sigfile %in% c("human.immune.CIBERSORT", "mouce.brain.ALLEN", "mouse.all.facs.TabulaMuris", "mouse.all.droplet.TabulaMuris")){
  signatures = sigfile
}else{
  signatures = read.table(sigfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
}
result = ATACRunSeurat(inputMat = countmatrix, project = prefix, method = "LSI")
result$ATAC = ATACAttachGenescore(ATAC = result$ATAC, RPmatrix = RPmatrix)
result$ATAC = ATACAnnotateCelltype(result$ATAC, signatures = signatures)
saveRDS(result, paste0(prefix, "_scATAC_Object.rds"))
result.tfs = ATACAnnotateTranscriptionFactor(ATAC = result$ATAC, peaks = result$peaks, project = prefix, giggle.path = gigglelib, organism = species)

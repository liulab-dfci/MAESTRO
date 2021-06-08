library(MAESTRO)
library(Seurat)
library(optparse)


option_list = list(
  make_option(c("--prefix"), type = "character", default = "MAESTRO",
              action = "store", help = "The prefix of the output files."
  ),
  make_option(c("--outdir"), type = "character", default = "Result/Analysis",
              action = "store", help = "The directory where the output files are stored."
  ),
  make_option(c("--atacobj"), type = "character", default = "",
              action = "store", help = "scATAC Seurat object."
  ),
  make_option(c("--rnaobj"), type = "character", default = "",
              action = "store", help = "scRNA Seurat object."
  )
)

argue = parse_args(OptionParser(option_list = option_list, usage = "Run integration pipeline"))

atacobj = argue$atacobj
rnaobj = argue$rnaobj
outprefix = argue$prefix
setwd(argue$outdir)

Seurat1 = readRDS(atacobj)
Seurat2 = readRDS(rnaobj)

Combine.res = Incorporate(ATAC = Seurat1$ATAC, RNA = Seurat2$RNA, RPmatrix = NULL, project = outprefix, dims.use = 1:30, RNA.res = 0.6, ATAC.res = 0.6)

Seurat.combined = Combine.res$CombinedObj
saveRDS(Seurat.combined, paste0(outprefix, "_integrate_Object.rds"))

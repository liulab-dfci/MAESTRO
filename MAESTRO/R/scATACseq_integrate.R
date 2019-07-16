argue = commandArgs(T)
cwd = getwd()
setwd(argue[5])
source("scATACseq_function.R")
setwd(file.path(cwd,argue[4]))

genescore = read.table(argue[1], header = TRUE, row.names = 1, sep = '\t')
SeuratRP <- PipelineSeurat(genescore, proj = paste0(argue[3], "_RP"), min.c = 10, min.g = 1000)
Seurat1 = paste0(argue[3], "_RP_SeuratObj.rds")
Seurat2 = argue[2]
proj.comb = paste0(argue[3], "_combined")
Seurat.combined = Incorparate(Seurat1 = Seurat1, Seurat2 = Seurat2, proj = proj.comb, type1 = "ATAC", type2 = "RNA", pcs.use = 30, dims.use = 1:20, res = 0.6, anchor = 2000)

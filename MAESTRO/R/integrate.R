library(MAESTRO)
library(Seurat)

argue = commandArgs(T)


Seurat1 = readRDS(argue[1])
Seurat2 = readRDS(argue[2])
setwd(argue[4])
Seurat.combined = Incorporate(ATAC = Seurat1, RNA = Seurat2, RPmatrix = NULL, project = argue[3], dims.use = 1:30, RNA.res = 0.6, ATAC.res = 0.6)
saveRDS(Seurat.combined, paste0(argue[3], "_integrate_Object.rds"))
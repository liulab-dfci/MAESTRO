#!/usr/bin/env Rscript

# example to merge three conditions
MAESTRO.rdsfiles <- c("1_scRNA_Object.rds",  "2_scRNA_Object.rds", "3_scRNA_Object.rds") 
MAESTRO.labels <- c("1","2","3")

# configuration
project.name <- "scRNA.merged"
dim.use <- 1:15
cluster.res <- 0.6
min.c <- 10
min.g <- 200
genes.test.use <- "presto"
genes.cutoff <- 1e-5
genes.pct <- 0.1
genes.logfc <- 0.25
only.pos <- FALSE

#=======code==========
library(MAESTRO)
library(Seurat)
library(Gmisc)

META <- data.frame( labels=MAESTRO.labels, rds=MAESTRO.rdsfiles )

scRNA.objects <- c( readRDS(MAESTRO.rdsfiles[1]), readRDS(MAESTRO.rdsfiles[2]), readRDS(MAESTRO.rdsfiles[3]) )

scRNA.merged <- merge(scRNA.objects[1]$RNA, y=c(scRNA.objects[2]$RNA, scRNA.objects[3]$RNA), add.cell.ids=META$labels, project=project.name)

# save the merged object
saveRDS(scRNA.merged,file=paste0(project.name,".rds"))

# tweak the parameters if necessary
scRNA.merged.res <- RNARunSeurat(inputMat=GetAssayData(object = scRNA.merged, slot = 'counts'), project = project.name, min.c = min.c, min.g = min.g, mito = FALSE, variable.genes = 2000, organism = "GRCh38", dims.use = dim.use, cluster.res = cluster.res, only.pos = FALSE, genes.test.use = genes.test.use, genes.cutoff = genes.cutoff, genes.pct = genes.pct, genes.logfc = genes.logfc)

# fix the orig.ident as the new labels
a <- scRNA.merged.res$RNA$orig.ident
n <- names(a)
t <- as.factor(gsub("_.*","",n))
t <- factor( t, levels=META$labels )
names(t) <- n
scRNA.merged.res$RNA$orig.ident <- t

# if cell type annotation is not needed, make `scRNA.merged.res$RNA@meta.data$assign.ident` as ``scRNA.merged.res$RNA@meta.data$seurat_clusters`
# Because the integration pipeline needs assign.ident values
scRNA.merged.res$RNA@meta.data$assign.ident <- scRNA.merged.res$RNA@meta.data$seurat_clusters

# save result before removing cell cycle effects
saveRDS(scRNA.merged.res,file=paste0(project.name,".res.before.cell.cycle.regression.rds"))

# optionally remove the cell cycle effects
#===========cell cycle===========
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scRNA.merged.res$RNA <- CellCycleScoring(scRNA.merged.res$RNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# will add a new group named Phase, and feature named S.Score and G2M.Score
scRNA.merged.res$RNA <- ScaleData(scRNA.merged.res$RNA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA.merged.res$RNA))

scRNA.merged.res$RNA <- fastDoCall("RunPCA", c(object = scRNA.merged.res$RNA, features = VariableFeatures(scRNA.merged.res$RNA)))
scRNA.merged.res$RNA <- RunUMAP(object = scRNA.merged.res$RNA, reduction = "pca", dims = dim.use)
scRNA.merged.res$RNA <- fastDoCall("FindNeighbors", c(object = scRNA.merged.res$RNA, reduction = "pca", dims = dim.use))
scRNA.merged.res$RNA <- fastDoCall("FindClusters", c(object = scRNA.merged.res$RNA, resolution = cluster.res))

scRNA.merged.res$RNA@meta.data$assign.ident <- scRNA.merged.res$RNA@meta.data$seurat_clusters

#=========DE analysis===========
cluster.genes <- FindAllMarkersMAESTRO(object = scRNA.merged.res$RNA, min.pct = genes.pct, logfc.threshold = genes.logfc, test.use = genes.test.use, only.pos = only.pos)
cluster.genes <- cluster.genes[cluster.genes$p_val_adj<genes.cutoff, ]
write.table(cluster.genes, paste0(scRNA.merged.res$RNA@project.name, "_DiffGenes_after_removing_cell_cycle.tsv"), quote = F, sep = "\t")

scRNA.merged.res$genes <- cluster.genes

# save result after removing cell cycle effects
saveRDS(scRNA.merged.res,file=paste0(project.name,".res.after.cell.cycle.regression.rds"))

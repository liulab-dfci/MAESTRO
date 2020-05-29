#!/usr/bin/env Rscript

# This script is intended to be an executable script
# Work is needed to make it into R methods

# example to merge three conditions
MAESTRO.rdsfiles <- c("1_scRNA_Object.rds",  "2_scRNA_Object.rds", "3_scRNA_Object.rds") 
MAESTRO.labels <- c("1","2","3")

# configuration
project.name <- "scRNA.merged"
organism <- "GRCh38"
dim.use <- 1:15
cluster.res <- 0.6
min.c <- 10
min.g <- 200
genes.test.use <- "presto"
genes.cutoff <- 1e-5
genes.pct <- 0.1
genes.logfc <- 0.25
only.pos <- FALSE
condadir <- "/mnt/lustre/users/tliu/miniconda3"
lisaenv <- "lisa"

#=======code==========
library(MAESTRO)
library(Seurat)
library(Gmisc)
library(purrr)

META <- data.frame( labels=MAESTRO.labels, rds=MAESTRO.rdsfiles )

scRNA.objects <- purrr::map( MAESTRO.rdsfiles, readRDS )

scRNA.merged <- merge(scRNA.objects[[1]]$RNA, y=purrr::map(scRNA.objects[2:length(scRNA.objects)],function(x){x$RNA}), add.cell.ids=MAESTRO.labels, project=project.name)

# save the merged object
saveRDS(scRNA.merged,file=paste0(project.name,".rds"))

# tweak the parameters if necessary
scRNA.merged.res <- RNARunSeurat(inputMat=GetAssayData(object = scRNA.merged, slot = 'counts'), project = project.name, min.c = min.c, min.g = min.g, mito = FALSE, variable.genes = 2000, organism = organism, dims.use = dim.use, cluster.res = cluster.res, only.pos = FALSE, genes.test.use = genes.test.use, genes.cutoff = genes.cutoff, genes.pct = genes.pct, genes.logfc = genes.logfc)

# fix the orig.ident as the new labels
a <- scRNA.merged$orig.ident
n <- names(a)
t <- as.factor(gsub("\d_.*","",n))
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

#========average expression for each cluster=========
cluster.averages <- AverageExpression(rna.aftercc$RNA)
write.table(x=cluster.averages[["RNA"]][, ], file="cluster.average.expression.tsv")

#=======Lisa analysis==========
scRNA.merged.tfs <- RNAAnnotateTranscriptionFactor(RNA = scRNA.merged.res$RNA, genes = scRNA.merged.res$genes, project = project.name, method = "LISA", lisa.mode = "local", conda.dir = condadir, lisa.envname = lisaenv, organism = organism, top.tf = 10)
saveRDS(scRNA.merged.tfs, file = paste0(project.name, ".tfs.rds"))

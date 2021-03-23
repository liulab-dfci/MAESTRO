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

ATAC.res = readRDS(atacobj)
RNA.res = readRDS(rnaobj)

setwd(argue$outdir)

cells.intersect = intersect(colnames(ATAC.res$ATAC), colnames(RNA.res$RNA))

Combined.Seurat = subset(RNA.res$RNA, cells = cells.intersect)
ATAC.subset = subset(ATAC.res$ATAC, cells = cells.intersect)
Combined.Seurat[["ATAC"]] = ATAC.subset[["ATAC"]]
Combined.Seurat[["ACTIVITY"]] = ATAC.subset[["ACTIVITY"]]

Combined.Seurat@graphs = c(Combined.Seurat@graphs, ATAC.subset@graphs)

for (i in names(Combined.Seurat@reductions)){
  reduc.colnames = colnames(Combined.Seurat@reductions[[i]]@cell.embeddings)
  reduc.colnames = paste0("rna", reduc.colnames)
  colnames(Combined.Seurat@reductions[[i]]@cell.embeddings) = reduc.colnames
  if (ncol(Combined.Seurat@reductions[[i]]@feature.loadings) > 0){
    colnames(Combined.Seurat@reductions[[i]]@feature.loadings) = reduc.colnames
  }
  Combined.Seurat@reductions[[i]]@key = paste0("rna", Combined.Seurat@reductions[[i]]@key)
}
names(Combined.Seurat@reductions) = paste0(names(Combined.Seurat@reductions), ".rna")

for (i in names(ATAC.subset@reductions)){
  reduc.colnames = colnames(ATAC.subset@reductions[[i]]@cell.embeddings)
  reduc.colnames = paste0("atac", reduc.colnames)
  colnames(ATAC.subset@reductions[[i]]@cell.embeddings) = reduc.colnames
  if (ncol(ATAC.subset@reductions[[i]]@feature.loadings) > 0){
    colnames(ATAC.subset@reductions[[i]]@feature.loadings) = reduc.colnames
  }
  ATAC.subset@reductions[[i]]@key = paste0("rna", ATAC.subset@reductions[[i]]@key)
}
names(ATAC.subset@reductions) = paste0(names(ATAC.subset@reductions), ".atac")

Combined.Seurat@reductions = c(Combined.Seurat@reductions, ATAC.subset@reductions)

colnames(ATAC.subset@meta.data) = gsub("assign.","ATAC_assign.",colnames(ATAC.subset@meta.data))
colnames(Combined.Seurat@meta.data) = gsub("assign.","RNA_assign.",colnames(Combined.Seurat@meta.data))

atac.meta.col = colnames(ATAC.subset@meta.data)[!(colnames(ATAC.subset@meta.data) %in% c("orig.ident", "assign.ident"))]
Combined.Seurat@meta.data = cbind(Combined.Seurat@meta.data, ATAC.subset@meta.data[rownames(Combined.Seurat@meta.data),atac.meta.col])

DefaultAssay(Combined.Seurat) = "ACTIVITY"
Combined.Seurat = RunPCA(Combined.Seurat, reduction.name = 'pca.activity')
Combined.Seurat = RunUMAP(Combined.Seurat, dims = 1:50, reduction = "pca.activity",reduction.name = 'umap.activity', reduction.key = 'activityUMAP_')

Combined.Seurat <- FindMultiModalNeighbors(Combined.Seurat, reduction.list = list("pca.rna", "pca.activity"), 
                                           dims.list = list(1:50, 1:50),
                                           knn.graph.name = "wknn",
                                           snn.graph.name = "wsnn",
                                           weighted.nn.name = "weighted.nn")
Combined.Seurat <- RunUMAP(Combined.Seurat, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_")
Combined.Seurat <- FindClusters(Combined.Seurat, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
saveRDS(Combined.Seurat, file.path(paste0(Combined.Seurat@project.name, "_multiome_Object.rds")))

p <- DimPlot(Combined.Seurat, reduction = "umap.wnn", group.by = "wsnn_res.0.8", pt.size = 0.2,label = TRUE, label.size = 3, repel = TRUE) + ggtitle("Joint clustering")
ggsave(file.path(paste0(Combined.Seurat@project.name, "_cluster_wsnn.png")), p, width=5.3, height=4.5)
p <- DimPlot(Combined.Seurat, reduction = "umap.wnn", group.by = "RNA_assign.ident", pt.size = 0.2,label = TRUE, label.size = 3, repel = TRUE) + ggtitle("Celltype annotation from scRNA")
ggsave(file.path(paste0(Combined.Seurat@project.name, "_annotated_wsnn.png")), p, width=6, height=4.5)



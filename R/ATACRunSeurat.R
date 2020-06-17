#' Clustering analysis for scATAC-seq data using Seurat
#'
#' Clustering analysis for scATAC-seq dataset using Seurat(version >=3.0.1). Including normalization, LSI/PCA dimension reduction, clustering and UMAP visualization. To run UMAP analysis, you must first install the umap-learn python package (e.g. via \code{pip install umap-learn}). 
#'
#' @docType methods
#' @name ATACRunSeurat
#' @rdname ATACRunSeurat
#'
#' @param inputMat Input unnormalized count matrix, with peaks as rows and cells as columns. Rownames should be like "chromosome_peakstart_peakend",
#' for example "chr10_100020591_100020841". 
#' @param type Type of the input matrix. Default is "matrix". Set to "object" if the input is Seurat object.
#' @param project Output project name. Default is "MAESTRO.scATAC.Seurat".
#' @param orign.ident Orginal identity for the input cells. If supplied, should keep the same order with the column name of the peak x cell matrix.
#' @param method Methods for dimension reduction, available options are LSI and PCA. Default is "LSI".
#' @param min.c Minimum number of cells required for a peak. Will exclude the peaks from input matrix if they only identified in 
#' less than \code{min.c} cells. Default is 10. See \code{\link{CreateSeuratObject}} for details.
#' @param min.p Minimum number of peaks required for a cell. Will exclude the cells from input matrix if less than \code{min.p}
#' peaks are deteced in one cell. Default is 100. See \code{\link{CreateSeuratObject}} for details.
#' @param dims.use Number of dimensions used for UMAP analysis. Default is 1:30, use the first 30 PCs.
#' @param cluster.res Value of the clustering resolution parameter. Please use a value above (below) 1.0 
#' if users want to obtain a larger (smaller) number of communities. Default is 0.6.
#' @param only.pos If seting true, only positive peaks will be output. Default is False.
#' @param peaks.test.use Denotes which test to use to identify differnetial peaks. Default is "presto", a fast version of Wilcoxon Rank Sum test. 
#' Available options are "wilcox" and "t". See \code{\link{FindAllMarkersMAESTRO}} for details.
#' @param peaks.cutoff Identify differential peaks with adjusted p.value less than \code{peaks.cutoff} as cluster specific peaks
#' @param peaks.pct Only test peaks that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing peaks that are very infrequently detected Default is 0.1
#' @param peaks.logfc Limit testing to peaks which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.2 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' for each cluster. Default cutoff is 1E-5.
#' @param runlsi.args Extra arguments passed to \code{\link{RunLSI}}.
#' @param runpca.args Extra arguments passed to \code{\link{RunPCA}}.
#' @param findneighbors.args Extra arguments passed to \code{\link{FindNeighbors}}.
#' @param findclusters.args Extra arguments passed to \code{\link{FindClusters}}.
#' @param \dots Extra arguments passed to \code{\link{RunUMAP}}.
#'
#' @author Chenfei Wang
#'
#' @return A list contains a Seurat object for ATAC clusters and a data frame of cluster specific peaks.
#'
#'
#' @examples
#' data(pbmc.ATAC)
#' pbmc.ATAC.res <- ATACRunSeurat(inputMat = pbmc.ATAC, project = "PBMC.scATAC.Seurat", method = "LSI")
#' str(pbmc.ATAC.res$ATAC)
#' head(pbmc.ATAC.res$peaks)
#'
#' @importFrom Seurat CreateSeuratObject DimPlot ElbowPlot FindClusters FindNeighbors NormalizeData RunLSI RunPCA RunUMAP ScaleData VariableFeatures
#' @importFrom ggplot2 ggsave
#' @importFrom Gmisc fastDoCall
#' @export

ATACRunSeurat <- function(inputMat, type = "matrix", project = "MAESTRO.scATAC.Seurat", orign.ident = NULL, 
                          min.c = 10, min.p = 100, method = "LSI", dims.use = 1:30, 
                          cluster.res = 0.6, only.pos = FALSE, peaks.test.use = "presto", 
                          peaks.cutoff = 1E-5, peaks.pct = 0.1, peaks.logfc = 0.2, 
                          runlsi.args = list(), runpca.args = list(), 
                          findneighbors.args = list(), findclusters.args = list(),...)
{
  if(type == "matrix"){SeuratObj <- CreateSeuratObject(inputMat, project = project, min.cells = min.c, min.features = min.p, assay = "ATAC")}
  if(type == "object"){SeuratObj <- inputMat}

  if(method == "LSI"){    
  #============ LSI ============
  message("LSI analysis ...")
  VariableFeatures(SeuratObj) <- names(which(Matrix::rowSums(SeuratObj) > min.c))
  # SeuratObj <- RunLSI(object = SeuratObj, scale.max = NULL) 
  SeuratObj <- fastDoCall("RunLSI", c(object = SeuratObj, runlsi.args)) 
  
  #============ UMAP ============
  message("UMAP analysis ...")
  #SeuratObj <- fastDoCall("RunUMAP", c(object = SeuratObj, dims = dims.use, reduction = "lsi", runumap.args))
  SeuratObj <- RunUMAP(object = SeuratObj, reduction = "lsi", dims = dims.use, ...)
  SeuratObj <- fastDoCall("FindNeighbors", c(object = SeuratObj, reduction = "lsi", dims = dims.use, findneighbors.args))
  SeuratObj <- fastDoCall("FindClusters", c(object = SeuratObj, resolution = cluster.res, findclusters.args))
  # SeuratObj <- FindNeighbors(object = SeuratObj, reduction = "lsi", dims = dims.use)
  # SeuratObj <- FindClusters(object = SeuratObj, resolution = cluster.res)
  p1 <- DimPlot(object = SeuratObj, pt.size = 0.5, label = TRUE)
  ggsave(file.path(paste0(project, "_cluster.png")), p1, width=5, height=4)
 
  #============ DE analysis ============
  message("Identify cluster specific peaks ...")
  SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  cluster.peaks <- NULL
  cluster.peaks <- FindAllMarkersMAESTRO(object = SeuratObj, min.pct = peaks.pct, logfc.threshold = peaks.logfc, test.use = peaks.test.use, only.pos = only.pos)
  cluster.peaks <- cluster.peaks[cluster.peaks$p_val_adj<peaks.cutoff, ]
  colnames(cluster.peaks)[7] <- "peak"
  write.table(cluster.peaks, paste0(project, "_DiffPeaks.tsv"), quote=F, sep="\t")}
  
  if(method == "PCA"){
  #============ PCA ============
  message("PCA analysis ...")
  SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  SeuratObj <- ScaleData(object = SeuratObj, var.to.regress="nCount_RNA")
  SeuratObj <- fastDoCall("RunPCA", c(object = SeuratObj, features = rownames(SeuratObj), runpca.args))
  
  # SeuratObj <- RunPCA(object = SeuratObj, features = rownames(SeuratObj))
  p2 = ElbowPlot(object = SeuratObj)
  ggsave(file.path(paste0(project,"_PCElbowPlot.png")), p2, width = 5, height = 4)
  
  #============ UMAP ============
  message("UMAP analysis ...")
  SeuratObj <- RunUMAP(object = SeuratObj, reduction = "pca", dims = dims.use, ...)
  SeuratObj <- fastDoCall("FindNeighbors", c(object = SeuratObj, reduction = "pca", dims = dims.use, findneighbors.args))
  SeuratObj <- fastDoCall("FindClusters", c(object = SeuratObj, resolution = cluster.res, findclusters.args))
  
  # SeuratObj <- FindNeighbors(object = SeuratObj, reduction = "pca", dims = dims.use)
  # SeuratObj <- FindClusters(object = SeuratObj, resolution = res)
  p3 <- DimPlot(object = SeuratObj, pt.size = 0.5, label = TRUE)
  ggsave(file.path(paste0(project, "_cluster.png")), p3, width=5, height=4)

  #============ DE analysis ============
  message("Identify cluster specific peaks ...")
  cluster.peaks <- NULL
  cluster.peaks <- FindAllMarkersMAESTRO(object = SeuratObj, min.pct = peaks.pct, logfc.threshold = peaks.logfc, test.use = peaks.test.use, only.pos = only.pos)
  cluster.peaks <- cluster.peaks[cluster.peaks$p_val_adj<peaks.cutoff, ]
  colnames(cluster.peaks)[7] <- "peak"
  write.table(cluster.peaks, paste0(proj, "_DiffPeaks.tsv"), quote=F, sep="\t")   
  }
  return(list(ATAC=SeuratObj, peaks=cluster.peaks))
}

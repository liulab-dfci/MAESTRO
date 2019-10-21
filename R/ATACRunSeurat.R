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
#' @param project Output project name. Default is "MAESTRO.scATAC.Seurat".
#' @param orign.ident Orginal identity for the input cells. If supplied, should keep the same order with the column name of the peak x cell matrix.
#' @param method Methods for dimension reduction, available options are LSI and PCA. Default is "LSI".
#' @param min.c Minimum number of cells required for a peak. Will exclude the peaks from input matrix if they only identified in 
#' less than \code{min.c} cells. Default is 50. See \code{link{CreateSeuratObject}} for details.
#' @param min.p Minimum number of peaks required for a cell. Will exclude the cells from input matrix if less than \code{min.p}
#' peaks are deteced in one cell. Default is 500. See \code{link{CreateSeuratObject}} for details.
#' @param dims.use Number of dimensions used for UMAP analysis. Default is 1:30, use the first 30 PCs.
#' @param cluster.res Value of the clustering resolution parameter. Default is 0.6.
#' @param peaks.test.use Denotes which test to use to identify differnetial peaks. Default is "wilcox". Available options are "bimod", "roc" and "t".
#' @param peaks.cutoff Identify differential peaks with adjusted p.value less than \code{peaks.cutoff} as cluster specific peaks
#' for each cluster. Default cutoff is 1E-5.
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
#' @export

ATACRunSeurat <- function(inputMat, project = "MAESTRO.scATAC.Seurat", orign.ident = NULL, method = "LSI", min.c = 50, min.p = 500, 
                          dims.use = 1:30, cluster.res = 0.6, peaks.test.use = "wilcox", peaks.cutoff = 1E-5)
{
  require(Seurat)
  require(ggplot2)
  SeuratObj <- CreateSeuratObject(inputMat, project = project, min.cells = min.c, min.features = min.p, assay = "ATAC")
  if(method == "LSI"){    
  #============ LSI ============
  message("LSI analysis ...")
  SeuratObj <- RunLSI(object = SeuratObj, n = 50, scale.max = NULL) 
  
  #============ UMAP ============
  message("UMAP analysis ...")
  SeuratObj <- RunUMAP(object = SeuratObj, reduction = "lsi", dims = dims.use)
  SeuratObj <- FindNeighbors(object = SeuratObj, reduction = "lsi", dims = dims.use)
  SeuratObj <- FindClusters(object = SeuratObj, resolution = cluster.res)
  p1 <- DimPlot(object = SeuratObj, pt.size = 0.5, label = TRUE)
  ggsave(file.path(paste0(project, "_cluster.png")), p1, width=5, height=4)
 
  #============ DE analysis ============
  message("Identify cluster specific peaks ...")
  SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  cluster.peaks <- NULL
  cluster.peaks <- FindAllMarkersMAESTRO(object = SeuratObj, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, test.use = peaks.test.use)
  cluster.peaks <- cluster.peaks[cluster.peaks$p_val_adj<peaks.cutoff, ]
  write.table(cluster.peaks, paste0(project, "_DiffPeaks.tsv"), quote=F, sep="\t")}
  
  if(method == "PCA"){
  #============ PCA ============
  message("PCA analysis ...")
  SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  SeuratObj <- ScaleData(object = SeuratObj, var.to.regress="nCount_RNA")
  SeuratObj <- RunPCA(object = SeuratObj, features = rownames(SeuratObj))
  p2 = ElbowPlot(object = SeuratObj)
  ggsave(file.path(paste0(project,"_PCElbowPlot.png")), p2, width = 5, height = 4)
  
  #============ UMAP ============
  message("UMAP analysis ...")
  SeuratObj <- RunUMAP(object = SeuratObj, reduction = "pca", dims = dims.use)
  SeuratObj <- FindNeighbors(object = SeuratObj, reduction = "pca", dims = dims.use)
  SeuratObj <- FindClusters(object = SeuratObj, resolution = res)
  p3 <- DimPlot(object = SeuratObj, pt.size = 0.5, label = TRUE)
  ggsave(file.path(paste0(project, "_cluster.png")), p3, width=5, height=4)

  #============ DE analysis ============
  message("Identify cluster specific peaks ...")
  cluster.peaks <- NULL
  cluster.peaks <- FindAllMarkersMAESTRO(object = SeuratObj, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, test.use = peaks.test.use)
  cluster.peaks <- cluster.peaks[cluster.peaks$p_val_adj<peaks.cutoff, ]
  colnames(cluster.peaks)[7] <- "peak"
  write.table(cluster.peaks, paste0(proj, "_DiffPeaks.tsv"), quote=F, sep="\t")   
  }
  return(list(ATAC=SeuratObj, peaks=cluster.peaks))
}

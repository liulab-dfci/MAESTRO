#' Clustering analysis for scATAC-seq data using scABC
#'
#' Clustering analysis for scATAC-seq dataset using scABC. Including normalization, scABC clustering and UMAP visualization. To run UMAP analysis, you must first install the umap-learn python package (e.g. via \code{pip install umap-learn}). 
#'
#' @docType methods
#' @name ATACRunscABC
#' @rdname ATACRunscABC
#'
#' @param inputMat Input unnormalized count matrix, with peaks as rows and cells as columns. Rownames should be like "chromosome_peakstart_peakend",
#' for example "chr10_100020591_100020841".
#' @param project Output project name. Default is "MAESTRO.scATAC.scABC".
#' @param min.c Minimum number of cells required for a peak. Will exclude the peaks from input matrix if they only identified in 
#' less than \code{min.c} cells. Default is 50.
#' @param min.p Minimum number of peaks required for a cell. Will exclude the cells from input matrix if less than \code{min.p}
#' peaks are deteced in one cell. Default is 500.
#' @param cluster.number Number of K for K-medoid clustering analysis. Default is 10.
#'
#' @author Chenfei Wang
#'
#' @return A list with clustering informations.
#'
#'
#' @examples
#' data(pbmc.ATAC)
#' pbmc.ATAC.res <- ATACRunscABC(inputMat = pbmc.ATAC, project = "PBMC.scATAC.scABC")
#' head(pbmc.ATAC.res)
#'
#' @import gplots
#' @importFrom scales hue_pal
#' @importFrom devtools source_url
#' @importFrom uwot umap
#' @export

ATACRunscABC <- function(inputMat, project = "MAESTRO.scATAC.scABC", min.c = 50, min.p = 500, cluster.number = 10)
{
  library(scABC)
  #============ Peak filtering and normalization ============
  message("Peak filtering and normalization ...")
  inputMat <- as.matrix(inputMat)
  peaks <- data.frame(chr=unlist(strsplit(rownames(inputMat),'_'))[seq(1,nrow(inputMat)*3,3)], 
                 start=as.numeric(unlist(strsplit(rownames(inputMat),'_'))[seq(2,nrow(inputMat)*3,3)]),
                 end=as.numeric(unlist(strsplit(rownames(inputMat),'_'))[seq(3,nrow(inputMat)*3,3)]))  
  scABCPeaksFiltered <- filterPeaks(inputMat, peaks, nreads_thresh = 1, ncells_thresh = min.c)
  scABCSamplesFiltered <- filterSamples(ForeGround = scABCPeaksFiltered$ForeGroundMatrix,
                         BackGround = matrix(nrow = dim(scABCPeaksFiltered$ForeGroundMatrix)[1], ncol = dim(scABCPeaksFiltered$ForeGroundMatrix)[2]), readsFGthresh = min.p)
  scABCForeGroundMatrix <- scABCSamplesFiltered$ForeGroundMatrix
  scABCPeaks <- scABCPeaksFiltered$peaks
  weights <- apply(scABCForeGroundMatrix, 2, mean)
  
  #============ Clustering ============
  message("Clustering ...")  
  scABCLandmarks <- computeLandmarks(scABCForeGroundMatrix, weights = weights, nCluster = cluster.number, nTop = 5000)
  scABCLandmarksAssignments <- assign2landmarks(scABCForeGroundMatrix, scABCLandmarks)
  ncols <- rev(hue_pal()(cluster.number))
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
  rowcols <- t(as.matrix(ncols[scABCLandmarksAssignments[order(scABCLandmarksAssignments)]]))
  rownames(rowcols) <- "cluster"
  scABCCell2LandmarkCorrelation <- NULL
  for(k in 1:cluster.number) scABCCell2LandmarkCorrelation <- cbind(scABCCell2LandmarkCorrelation, apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,k], method = 'spearman')))
  Normalization <- rowMeans(abs(scABCCell2LandmarkCorrelation))
  scABCCell2LandmarkCorrelationNormalized <- scABCCell2LandmarkCorrelation/Normalization
  scABCCell2LandmarkCorrelationNormalized[scABCCell2LandmarkCorrelationNormalized<0] <- 0
  scABCCell2LandmarkCorrelationNormalized <- cbind(scABCCell2LandmarkCorrelationNormalized, scABCLandmarksAssignments)
  scABCCell2LandmarkCorrelationNormalized <- scABCCell2LandmarkCorrelationNormalized[order(scABCCell2LandmarkCorrelationNormalized[,cluster.number+1]),1:cluster.number]
  png(paste0(project,"_heatmap.png"),width=8,height=6)
  heatmap.3(scABCCell2LandmarkCorrelationNormalized, dendrogram='none', Rowv=FALSE, Colv=FALSE,
           trace='none', col = scalered, margin = c(5, 5), density.info = "none",
           RowSideColors = rowcols, RowSideColorsSize=2, symm=F,symkey=F,
           symbreaks=F, scale="none", main = project)
  legend("bottomleft", legend = paste0("Cluster ", 1:cluster.number), col = ncols, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7, pch = 15)
  dev.off()
 
  #============ UMAP visualization ============
  message("UMAP visualization ...")  
  d_umap <- as.dist(1 - cor(scABCForeGroundMatrix, method = "spearman"))
  p_umap <- umap(d_umap, pca = 50, n_threads = 8)
  scABC_col <- t(as.matrix(ncols[scABCLandmarksAssignments]))
  png(paste0(project,"_umap.png"),width=6,height=6)
  plot(p_umap[,1], p_umap[,2], col=scABC_col, type = "p", pch = 20, xlab = "UMAP_1", ylab = "UMAP_2", main = project)
  legend("bottomright", legend = 1:cluster.number, col=ncols, border=FALSE, bty="n", y.intersp = 1, cex= 1, pch = 16)
  dev.off()
  
  return(scABCLandmarksAssignments)
}

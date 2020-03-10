#' Clustering analysis for scATAC-seq data using cisTopic
#'
#' Clustering analysis for scATAC-seq dataset using cisTopic(version >= 0.2.0). Including topic modeling, clustering based on cisTopics and UMAP visualization. To run UMAP analysis, you must first install the umap-learn python package (e.g. via \code{pip install umap-learn}). 
#'
#' @docType methods
#' @name ATACRuncisTopic
#' @rdname ATACRuncisTopic
#'
#' @param inputMat Input unnormalized count matrix, with peaks as rows and cells as columns. Rownames should be like "chromosome_peakstart_peakend",
#' for example "chr10_100020591_100020841".
#' @param project Output project name for the cisTopic object. Default is "MAESTRO.scATAC.cisTopic".
#' @param min.c Minimum number of cells required for a peak. Will exclude the peaks from input matrix if they only identified in 
#' less than \code{min.c} cells. Default is 50.
#' @param min.p Minimum number of peaks required for a cell. Will exclude the cells from input matrix if less than \code{min.p}
#' peaks are deteced in one cell. Default is 500.
#' @param topics.number List number of topics to calculate. Default is 5,10,20,30. 
#' @param topics.select Seleted number of topics. Default is 30.
#'
#' @author Chenfei Wang
#'
#' @return A list contains a cisTopic object for ATAC clusters and a list with clustering informations.
#'
#'
#' @examples
#' data(pbmc.ATAC)
#' pbmc.ATAC.res <- ATACRuncisTopic(inputMat = pbmc.ATAC, project = "PBMC.scATAC.cisTopic")
#' str(pbmc.ATAC.res$ATAC)
#' head(pbmc.ATAC.res$cluster)
#'
#' @importFrom scales hue_pal
#' @export

ATACRuncisTopic <- function(inputMat, project = "MAESTRO.scATAC.cisTopic", min.c = 50, min.p = 500, topics.number = c(5,10,20,30), topic.select = 30)
{
  library(cisTopic)
  library(fpc)
  
  #============ Peak filtering ============
  message("Peak filtering ...")
  inputMat = rbind(inputMat, apply(inputMat,2,sum))
  inputMat = inputMat[1:(nrow(inputMat)-1),which(inputMat[nrow(inputMat),]>min.p)]
  inputMat = cbind(inputMat, apply(inputMat,1,sum))
  inputMat = inputMat[which(inputMat[,ncol(inputMat)]>min.c),1:(ncol(inputMat)-1)]
  rownames(inputMat) <- sub("\\_", "\\:", rownames(inputMat))
  rownames(inputMat) <- sub("\\_", "\\-", rownames(inputMat))
  
  #============ Identify topics ============
  message("Identify topics ...")
  cisTopicObject <- createcisTopicObject(as.matrix(inputMat), project.name=project)
  cisTopicObject <- runModels(cisTopicObject, topic=topics.number, seed=1000, nCores=10, burnin = 120, iterations = 150, addModels=FALSE)
  png(paste0(project,"_topics_select.png"),width=10,height=5)
  par(mfrow=c(1,2))
  cisTopicObject <- selectModel(cisTopicObject, select = topic.select)
  logLikelihoodByIter(cisTopicObject, select=topics.number)
  dev.off()
  png(paste0(project,"_topics_heatmap.png"),width=12,height=5)
  cellTopicHeatmap(cisTopicObject, method='Probability')
  dev.off()  
  cisTopicObject <- runUmap(cisTopicObject, target='cell', perplexity=200, check_duplicates=FALSE)
  png(paste0(project,"_topics_umap.png"),width=10,height=5)
  par(mfrow=c(2,5))
  plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
  dev.off()    
  
  #============ UMAP visualization ============
  message("UMAP visualization ...")
  dr <- cisTopicObject@dr$cell$Umap
  dbclust <- dbscan(dr, MinPts=10, showplot=F, method="raw", eps=0.75)
  cellassign <- dbclust$cluster
  ncols <- rev(hue_pal()(max(cellassign)+1))
  cisTopic_col <- ncols[cellassign]
  png(paste0(project,"_cluster.png"),width=6,height=6)
  plot(dr[,1], dr[,2], col=cisTopic_col, type = "p", pch = 20, xlab = "UMAP_1", ylab = "UMAP_2", main = project)
  legend("bottomright", legend = 1:(max(cluster)+1), col=ncols, border=FALSE, bty="n", y.intersp = 1, cex= 1, pch = 16)
  dev.off()    
  names(cellassign) <- cisTopicObject@cell.names
  metadata <- data.frame(cluster=as.character(cellassign), row.names=names(cellassign))
  cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = metadata)
  return(list(ATAC=cisTopicObject, cluster=cellassign)) 
}
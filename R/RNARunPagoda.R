#' Clustering analysis for scRNA-seq data using Pagoda2
#'
#' Clustering analysis for scRNA-seq dataset using Pagoda2 (version 0.1.0). Including feature selection, dimension reduction, clustering and tSNE visualization. To run this function, you must first install Pagoda2 (see from https://github.com/hms-dbmi/pagoda2).
#'
#' @docType methods
#' @name RNARunPagoda
#' @rdname RNARunPagoda
#'
#' @param inputMat Input unnormalized matrix, with genes as rows and cells as columns. Could be count matrix or TPM, FPKM matrix.
#' @param project Output project name. Default is "MAESTRO.scRNA.Pagoda".
#' @param min.c Minimum number of cells required for a gene. Will exclude the genes from input matrix if they only expressed in 
#' less than \code{min.c} cells. Default is 10.
#' @param min.g Minimum number of genes required for a cell. Will exclude the cells from input matrix if less than \code{min.g}
#' genes are deteced in one cell. Default is 200.
#' @param nPCs Number of PCs considered in the clustering analysis. Default is 100. 
#' @param nKs Number of K for the KNN graph. Default is 30.
#'
#' @author Chenfei Wang
#'
#' @return A list contains a Pagoda2 object and a list for the cluster specific genes.
#'
#'
#' @examples
#' data(pbmc.RNA)
#' pbmc <- RNARunPagoda(inputMat = pbmc.RNA, project = "PBMC.scRNA.Pagoda2")
#' pbmc$RNA
#' head(pbmc$genes)
#'
#' @export

RNARunPagoda <- function(inputMat, project = "MAESTRO.scRNA.Pagoda", min.c = 10, min.g = 200, nPCs = 100, nKs = 30)
{
  library(pagoda2)
  library(igraph)
  
  PagodaObj <- Pagoda2$new(inputMat,modelType='plain',trim=10,log.scale=T)
  PagodaObj$adjustVariance(plot=T,do.par=T,gam.k=10)
  PagodaObj$calculatePcaReduction(nPcs=nPCs,n.odgenes=3e3,maxit=1000)
  PagodaObj$makeKnnGraph(k=nKs,type='PCA',center=T,distance='cosine');
  PagodaObj$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
  PagodaObj$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
  png(paste0(project, "_cluster.png"), width=4.8, height=4)
  PagodaObj$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main=project)
  dev.off()
  deSets <- get.de.geneset(PagodaObj, groups = PagodaObj$clusters$PCA[[1]], prefix = 'de_')
  saveRDS(deSets, paste0(project, "_DiffGenes.rds"))
  return(list(RNA=PagodaObj, genes=deSets))
}


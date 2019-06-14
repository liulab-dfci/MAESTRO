#' Preprocess and analysis scATAC-seq data
#' Function declare
#' @Run_scABC
#' @Run_cisTopic
#' @Run_snapATAC
#' @Run_cicero
#' @Run_chromVAR

Run_scABC <- function(peakLocFile, peakCountFile, outPrefix, minPeak = 500)
{
  library(scABC)
  start_time <- Sys.time()
  peakLoc <- read.table(peakLocFile)
  peakCount <- as.matrix(read.table(peakCountFile))
  scABCPeaksFiltered <- filterPeaks(peakCount, peakLoc, nreads_thresh = 1, ncells_thresh = 10)
  scABCSamplesFiltered <- filterSamples(ForeGround = scABCPeaksFiltered$ForeGroundMatrix,
                                       BackGround = matrix(nrow = dim(scABCPeaksFiltered$ForeGroundMatrix)[1], ncol = dim(scABCPeaksFiltered$ForeGroundMatrix)[2]), readsFGthresh = minPeak)
  scABCForeGroundMatrix <- scABCSamplesFiltered$ForeGroundMatrix
  scABCPeaks <- scABCPeaksFiltered$peaks
  weights <- apply(scABCForeGroundMatrix, 2, mean)
  scABCLandmarks <- computeLandmarks(scABCForeGroundMatrix, weights = weights, nCluster = 10, nTop = 5000)
  scABCLandmarksAssignments <- assign2landmarks(scABCForeGroundMatrix, scABCLandmarks)
  
  # Heatmap visualization
  library(gplots)
  library(devtools)
  library(RColorBrewer)
  ncols <- brewer.pal(10, "Set3")
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
  rowcols <- t(as.matrix(ncols[scABCLandmarksAssignments[order(scABCLandmarksAssignments)]]))
  rownames(rowcols) <- "cluster"
  scABCCell2LandmarkCorrelation <- NULL
  scABCCell2LandmarkCorrelation <- cbind(apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,1], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,2], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,3], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,4], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,5], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,6], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,7], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,8], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,9], method = 'spearman')),
                                         apply(scABCForeGroundMatrix, 2, function(x) cor(x, scABCLandmarks[,10], method = 'spearman')))
  Normalization <- rowMeans(abs(scABCCell2LandmarkCorrelation))
  scABCCell2LandmarkCorrelationNormalized <- scABCCell2LandmarkCorrelation/Normalization
  scABCCell2LandmarkCorrelationNormalized[scABCCell2LandmarkCorrelationNormalized<0] <- 0
  scABCCell2LandmarkCorrelationNormalized <- cbind(scABCCell2LandmarkCorrelationNormalized, scABCLandmarksAssignments)
  scABCCell2LandmarkCorrelationNormalized <- scABCCell2LandmarkCorrelationNormalized[order(scABCCell2LandmarkCorrelationNormalized[,11]),1:10]
  pdf(paste0(outPrefix,"_heatmap.pdf"),width=8,height=6)
  heatmap.3(scABCCell2LandmarkCorrelationNormalized, dendrogram='none', Rowv=FALSE, Colv=FALSE,
           trace='none', col = scalered, margin = c(5, 5), density.info = "none",
           RowSideColors = rowcols, RowSideColorsSize=2, symm=F,symkey=F,
           symbreaks=F, scale="none", main = outPrefix)
  legend("bottomleft", legend = paste0("cluster ", 1:10), col = ncols, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7, pch = 15)
  dev.off()
 
  # Umap visualization
  library(uwot)
  d_umap <- as.dist(1 - cor(scABCForeGroundMatrix, method = "spearman"))
  p_umap <- umap(d_umap, pca = 50, n_threads = 8)
  scABC_col <- t(as.matrix(ncols[scABCLandmarksAssignments]))
  pdf(paste0(outPrefix,"_umap.pdf"),width=6,height=6)
  plot(p_umap[,1], p_umap[,2], col=scABC_col, type = "p", pch = 20, xlab = "Dim 1", ylab = "Dim 2", main = outPrefix)
  legend("bottomright", legend = c(paste0("cluster ", 1:10)), col=ncols, border=FALSE, bty="n", y.intersp = 1, cex= 1, pch = 16)
  dev.off()

  saveRDS(p_umap, paste0(outPrefix,"_umap.rds"))
  saveRDS(scABCLandmarksAssignments, paste0(outPrefix,"_cluster.rds"))
  end_time <- Sys.time()
  end_time - start_time
}

Run_cisTopic <- function(peakCountFile, outPrefix, minPeak = 500)
{
  library(cisTopic)
  start_time <- Sys.time()
  peakCount = read.table(peakCountFile)
  peakCount = rbind(peakCount, apply(peakCount,2,sum))
  peakCount = peakCount[1:(nrow(peakCount)-1),which(peakCount[nrow(peakCount),]>minPeak)]
  rownames(peakCount) <- sub("\\_", "\\:", rownames(peakCount))
  rownames(peakCount) <- sub("\\_", "\\-", rownames(peakCount))
  cisTopicObject <- createcisTopicObject(as.matrix(peakCount), project.name=outPrefix)
  cisTopicObject <- runModels(cisTopicObject, topic=c(5,10,20,30), seed=987, nCores=10, burnin = 120, iterations = 150, addModels=FALSE)
  
  # Topic selection and visualization
  pdf(paste0(outPrefix,"_topics_select.pdf"),width=10,height=5)
  par(mfrow=c(1,2))
  cisTopicObject <- selectModel(cisTopicObject, select = 10)
  logLikelihoodByIter(cisTopicObject, select=c(5,10,20,30))
  dev.off()
  pdf(paste0(outPrefix,"_topics_heatmap.pdf"),width=12,height=5)
  cellTopicHeatmap(cisTopicObject, method='Probability')
  dev.off()  
  cisTopicObject <- runUmap(cisTopicObject, target='cell', perplexity=200, check_duplicates=FALSE)
  pdf(paste0(outPrefix,"_topics_umap.pdf"),width=10,height=5)
  par(mfrow=c(2,5))
  plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
  dev.off()    
  saveRDS(cisTopicObject, paste0(outPrefix,".rds")) 
  
  # Umap visualization
  cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
  cellassign <- t(scale(t(cellassign)))
  cellassignmax <- apply(cellassign, 2, which.max)
# metadata <- data.frame(celltype=sapply(strsplit(names(cellassignmax),"\\."),function(x) x[2]), cluster=as.character(cellassignmax), row.names=names(cellassignmax))
  metadata <- data.frame(cluster=as.character(cellassignmax), row.names=names(cellassignmax))
  cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = metadata)
  pdf(paste0(outPrefix,"_cluster_umap.pdf"),width=6,height=5)
  plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy='cluster', cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)
  dev.off()
    
  saveRDS(cellassignmax, paste0(outPrefix,"_cluster.rds")) 
  end_time <- Sys.time()
  end_time - start_time
}

Run_cicero <- function(peakCountFile, outPrefix)
{
  library(cicero)
  start_time <- Sys.time()
  peakCount = read.table(peakCountFile)
  colnames(peakCount) = c("Peak","Cell","Count")
  input_cds = make_atac_cds(peakCount, binarize = TRUE)
  
  set.seed(2017)
  input_cds = detectGenes(input_cds)
  input_cds = estimateSizeFactors(input_cds)
  input_cds = reduceDimension(input_cds, max_components = 2, num_dim=10, reduction_method = 'tSNE', norm_method = "none")
  tsne_coords = t(reducedDimA(input_cds))
  row.names(tsne_coords) = row.names(pData(input_cds))
  cicero_cds = make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)
   
  hg38.genome = read.table('/homes/cwang/annotations/hg38/hg38.len')
  conns = run_cicero(cicero_cds, hg38.genome,  window = 5e+05, sample_num = 100) # Takes a few minutes to run
  gene = read.table('/homes/cwang/annotations/hg38/hg38.refseq.bed')
  gene_annotation_sub = unique(gene[,c(1:3,5)])
  names(gene_annotation_sub)[4] = "gene"
  input_cds = annotate_cds_by_site(input_cds, gene_annotation_sub)
  unnorm_ga = build_gene_activity_matrix(input_cds, conns)
  num_genes = pData(input_cds)$num_genes_expressed
  names(num_genes) = row.names(pData(input_cds))
  cicero_gene_activities = normalize_gene_activities(unnorm_ga, num_genes)
  write.table(as.matrix(cicero_gene_activities),paste0(outPrefix,".txt"),sep="\t",quote=F)
  end_time <- Sys.time()
  end_time - start_time
}

Run_chromVAR <- function(peakLocFile, peakCountFile, outPrefix, minPeak = 500)
{
  library(chromVAR)
  library(motifmatchr)
  library(Matrix)
  library(SummarizedExperiment)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(chromVARmotifs)
  data("human_pwms_v2")
  
  start_time <- Sys.time()
  peakLoc <- read.table(peakLocFile)
  rownames(peakLoc) <- paste0(peakLoc[,1],'_',peakLoc[,2],'_',peakLoc[,3])
  peakCount <- read.table(peakCountFile)
  colData <- DataFrame(celltype=colnames(peakCount), depth=apply(peakCount,2,sum))
  colData <- colData[which(colData[,2]>minPeak),]
  rowData <- DataFrame(peaks=rownames(peakCount), depth=apply(peakCount[,colData[,1]],1,sum))
  rowData = rowData[which(rowData[,2]>0),]
  finalPeak = as.matrix(peakLoc[rowData[,1],])
  finalCount = as.matrix(peakCount[rowData[,1],colData[,1]])
  rowRanges = GRanges(as.character(finalPeak[,1]),IRanges(as.integer(finalPeak[,2]),as.integer(finalPeak[,3])),strand="*",score=as.integer(5),qval=1,name=rownames(finalPeak))
    
  chromVAR_counts  = SummarizedExperiment(assays=SimpleList(counts=finalCount),rowRanges=rowRanges, colData=colData)
  chromVAR_counts = addGCBias(chromVAR_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
  motif_ix = matchMotifs(human_pwms_v2, chromVAR_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
  dev = computeDeviations(object = chromVAR_counts, annotations = motif_ix)
  saveRDS(dev, paste0(outPrefix,"_dev.rds"))

  variability = computeVariability(dev)
  zscore = assays(dev)$z
  rownames(zscore) = variability[,1]
  write.table(zscore,paste0(outPrefix,"_zscore.txt"),sep="\t",quote=F)
  end_time <- Sys.time()
  end_time - start_time
}

#' Preprocess and analysis scATAC-seq data
#' Function declare
#' @FragPlot
#' @MapPlot_micro
#' @FripPlot_micro
#' @MapPlot_10x
#' @FripPlot_10x
#' @PipelineSeurat

#' Quality Control Function 
FragPlot <- function(matrix,prefix)
{
  pdf(paste0(prefix,"_frag.pdf"),width=6,height=6)
  D <- density(abs(matrix[abs(matrix)<=1000&abs(matrix)>0]))
  plot(D$x,D$y*D$n/1000,col="blue",lwd=2,type="l",main=prefix,xlab="Fragment size",ylab="Fragment count (K)")
  dev.off()
}

MapPlot_micro <- function(matrix, prefix)
{
  pdf(paste0(prefix,"_map.pdf"),width=6,height=6)
  par(mfrow=c(2,2))
  hist(log10(matrix[which(matrix[,1]>1000),2]), border=NA, col="blue",main="Total Fragments",xlim=c(3,5.5),xlab="log10(Total Fragments)",ylab="Frequency")
  hist(matrix[,2]/matrix[,1], border=NA, col="blue",main="Mapped Ratio",xlim=c(0,1),xlab="Mapped Ratio",ylab="Frequency")
  hist(matrix[,3]/matrix[,2], border=NA, col="blue",main="Duplicate Ratio",xlim=c(0,1),xlab="Duplicate Ratio",ylab="Frequency")
  hist(matrix[,4]/matrix[,2], border=NA, col="blue",main="Mitochondria Ratio",xlim=c(0,1),xlab="Mitochondria Ratio",ylab="Frequency")
  dev.off()
}
FripPlot_micro <- function(matrix, prefix, reads_cutoff = 1000, frip_cutoff = 0.05)
{
  pdf(paste0(prefix,"_frip.pdf"),width=6,height=6)
  plot(log10(matrix[which(matrix[,5]<reads_cutoff|(matrix[,6]/matrix[,5])<frip_cutoff),5]+1),matrix[which(matrix[,5]<reads_cutoff|(matrix[,6]/matrix[,5])<frip_cutoff),6]/matrix[which(matrix[,5]<reads_cutoff|(matrix[,6]/matrix[,5])<frip_cutoff),5],
   	   xlim=c(0,5),ylim=c(0,1),pch='.',col='blue',ylab='Fraction of promoter reads',xlab='Reads passed filters (log10)',main=prefix)
       points(log10(matrix[which(matrix[,5]>=reads_cutoff&(matrix[,6]/matrix[,5])>=frip_cutoff),5]+1),matrix[which(matrix[,5]>=reads_cutoff&(matrix[,6]/matrix[,5])>=frip_cutoff),6]/matrix[which(matrix[,5]>=reads_cutoff&(matrix[,6]/matrix[,5])>=frip_cutoff),5],
       pch='.',col='red')
  legend("topright",c("cells","non-cells"),col=c("red","blue"),pch=20,box.lty=0)
  dev.off()
  write.table(rownames(matrix[which(matrix[,5]>=reads_cutoff&(matrix[,6]/matrix[,5])>=frip_cutoff),]), paste0(prefix,"_barcodes.txt"), sep = "\n", quote=F, row.names=F, col.names=F)
}
MapPlot_10x <- function(matrix, prefix)
{
  pdf(paste0(prefix,"_map.pdf"),width=6,height=6)
  par(mfrow=c(2,2))
  hist(log10(matrix[which(matrix[,2]>1000),2]), border=NA, col="blue",main="Total Fragments",xlim=c(3,5.5),xlab="log10(Total Fragments)",ylab="Frequency")
  hist(1-matrix[2:nrow(matrix),5]/matrix[2:nrow(matrix),2], border=NA, col="blue",main="Mapped Ratio",xlim=c(0,1),xlab="Mapped Ratio",ylab="Frequency")
  hist(matrix[2:nrow(matrix),3]/matrix[2:nrow(matrix),2], border=NA, col="blue",main="Duplicate Ratio",xlim=c(0,1),xlab="Duplicate Ratio",ylab="Frequency")
  hist(matrix[2:nrow(matrix),7]/matrix[2:nrow(matrix),2], border=NA, col="blue",main="Mitochondria Ratio",xlim=c(0,1),xlab="Mitochondria Ratio",ylab="Frequency")
  dev.off()
}

FripPlot_10x <- function(matrix, prefix, reads_cutoff = 1000, frip_cutoff = 0.05)
{
  pdf(paste0(prefix,"_frip.pdf"),width=6,height=6)
  plot(log10(matrix[which(matrix[,8]<reads_cutoff|(matrix[,14]/matrix[,8])<frip_cutoff),8]+1),matrix[which(matrix[,8]<reads_cutoff|(matrix[,14]/matrix[,8])<frip_cutoff),14]/matrix[which(matrix[,8]<reads_cutoff|(matrix[,14]/matrix[,8])<frip_cutoff),8],
   	   xlim=c(0,5),ylim=c(0,1),pch='.',col='blue',ylab='Fraction of promoter reads',xlab='Reads passed filters (log10)',main=prefix)
       points(log10(matrix[which(matrix[,8]>=reads_cutoff&(matrix[,14]/matrix[,8])>=frip_cutoff),8]+1),matrix[which(matrix[,8]>=reads_cutoff&(matrix[,14]/matrix[,8])>=frip_cutoff),14]/matrix[which(matrix[,8]>=reads_cutoff&(matrix[,14]/matrix[,8])>=frip_cutoff),8],
       pch='.',col='red')
  legend("topright",c("cells","non-cells"),col=c("red","blue"),pch=20,box.lty=0)
  dev.off()
  write.table(as.character(matrix[which(matrix[,8]>=reads_cutoff&(matrix[,14]/matrix[,8])>=frip_cutoff),1]), paste0(prefix,"_barcodes.txt"), sep = "\n", quote=F, row.names=F, col.names=F)
}

#' Analysis Function
library(Seurat)
library(dplyr)

PipelineSeurat = function(countMat, proj, min.c = 10, min.p = 200, max.p = 100000, nfeatures = 25000, org="GRCh38",
  do.norm = T, normalization.method = "LogNormalize", dims.use = 1:15, res = 0.6, diff.p = F, diff.co = 1E-5)
{
  start_time <- Sys.time()  
  #=========QC========
  message("Check gene and cell coverage ...")
  nGene = apply(countMat, 2, function(x) length(x[x>0]))
  nCell = apply(countMat, 1, function(x) length(x[x>0]))
  pdf(paste0(proj,"_Coverage.pdf"),width=8,height=4.5)
  par(mfrow=c(1,2))
  plot(1:ncol(countMat),sort(nGene),pch=16,col="blue",ylab="Number of Peaks Detected",xlab="Cells",main="Cell Filter")
  abline(h=min.p,lwd=2,lty=2);text(ncol(countMat)/2,min.p+max(nGene)*0.05,paste0("n = ",min.p));legend("topleft",paste0("ave Peaks = ",round(mean(nGene))),box.lty=0)
  plot(1:nrow(countMat),sort(nCell),pch=16,col="blue",ylab="Number of Cells Detected",xlab="Peaks",main="Peak Filter")
  abline(h=min.c,lwd=2,lty=2);text(nrow(countMat)/2,min.c+max(nCell)*0.05,paste0("n = ",min.c));legend("topleft",paste0("ave Cells = ",round(mean(nCell))),box.lty=0)
  dev.off()
  SeuratObj = CreateSeuratObject(countMat, project = proj, min.cells = min.c, min.penes = min.p)
  
  #=========Filter========
  message("Filter cells and find variable peaks ...")  
  SeuratObj = FilterCells(object = SeuratObj, subset.names = c("nGene"), low.thresholds = min.p, high.thresholds = max.p) 
  if (do.norm) SeuratObj = NormalizeData(object = SeuratObj, normalization.method = normalization.method, scale.factor = 10000)
  SeuratObj = FindVariableGenes(object = SeuratObj, mean.function = ExpMean, dispersion.function = LogVMR,
               x.low.cutoff = 0.01, x.high.cutoff = 8, y.cutoff = 0.9, do.plot = FALSE, top.genes = nfeatures)
  SeuratObj = ScaleData(object = SeuratObj)
  
  #=========PCA===========
  message("PCA analysis ...")
  SeuratObj = RunPCA(object = SeuratObj, pc.genes = SeuratObj@var.genes, do.print = FALSE, rev.pca = TRUE)
  pdf(file.path(paste0(proj,"_PCElbowPlot.pdf")), width = 5, height = 4)
  PCElbowPlot(object = SeuratObj)
  dev.off()
  
  #=========tSNE===========
  message("t-SNE analysis ...")
  SeuratObj = FindClusters(object = SeuratObj, reduction.type = "pca", dims.use = dims.use,
                resolution = res, print.output = 0, save.SNN = TRUE)
  SeuratObj = RunTSNE(object = SeuratObj, dims.use = dims.use, do.fast = TRUE)
  pdf(file.path(paste0(proj, "_tSNE_cluster.pdf")), width=5, height=4)
  TSNEPlot(object = SeuratObj, do.label = TRUE, pt.size = 0.5, group.by = paste0("res.",res))
  dev.off()

  #=========identify marker===========
  message("Find cluster specific peaks ...")
  cluster.specific.peaks = NULL
  if (diff.p){
  cluster.specific.peaks = FindAllMarkers(object = SeuratObj, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
  cluster.specific.peaks = cluster.specific.peaks[cluster.specific.peaks$p_val_adj<diff.co, ]
  saveRDS(cluster.specific.peaks, file.path(paste0(proj, "_ClusterSpecificPeaks.rds")))}
  
  saveRDS(SeuratObj, file.path(paste0(proj, "_SeuratObj.rds")))
  end_time <- Sys.time()
  end_time - start_time  
}

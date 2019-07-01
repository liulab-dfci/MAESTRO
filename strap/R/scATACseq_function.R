#' Preprocess and analysis scATAC-seq data
#' Seurat 3.0.0
#' Function declare
#' @FragPlot
#' @MapPlot_micro
#' @FripPlot_micro
#' @MapPlot_10x
#' @FripPlot_10x
#' @PipelineSeurat
#' @Incorparate

srcdir = "../"
Rdir = "./"


#' Quality Control Function 
FragPlot <- function(matrix,prefix)
{
  png(paste0(prefix,"_frag.png"),width=6,height=6, res = 300, units = "in")
  D <- density(abs(matrix[abs(matrix)<=1000&abs(matrix)>0]))
  plot(D$x,D$y*D$n/1000,col="blue",lwd=2,type="l",main=prefix,xlab="Fragment size",ylab="Fragment count (K)")
  dev.off()
}

MapPlot_micro <- function(matrix, prefix)
{
  png(paste0(prefix,"_map.png"),width=6,height=6, res = 300, units = "in")
  par(mfrow=c(2,2))
  hist(log10(matrix[which(matrix[,1]>1000),2]), border=NA, col="blue",main="Total Fragments",xlim=c(3,5.5),xlab="log10(Total Fragments)",ylab="Frequency")
  hist(matrix[,2]/matrix[,1], border=NA, col="blue",main="Mapped Ratio",xlim=c(0,1),xlab="Mapped Ratio",ylab="Frequency")
  hist(matrix[,3]/matrix[,2], border=NA, col="blue",main="Duplicate Ratio",xlim=c(0,1),xlab="Duplicate Ratio",ylab="Frequency")
  hist(matrix[,4]/matrix[,2], border=NA, col="blue",main="Mitochondria Ratio",xlim=c(0,1),xlab="Mitochondria Ratio",ylab="Frequency")
  dev.off()
}
FripPlot_micro <- function(matrix, prefix, reads_cutoff = 1000, frip_cutoff = 0.05)
{
  png(paste0(prefix,"_frip.png"),width=6,height=6, res = 300, units = "in")
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
  png(paste0(prefix,"_map.png"),width=6,height=6, res = 300, units = "in")
  par(mfrow=c(2,2))
  hist(log10(matrix[which(matrix[,2]>1000),2]), border=NA, col="blue",main="Total Fragments",xlim=c(3,5.5),xlab="log10(Total Fragments)",ylab="Frequency")
  hist(1-matrix[2:nrow(matrix),5]/matrix[2:nrow(matrix),2], border=NA, col="blue",main="Mapped Ratio",xlim=c(0,1),xlab="Mapped Ratio",ylab="Frequency")
  hist(matrix[2:nrow(matrix),3]/matrix[2:nrow(matrix),2], border=NA, col="blue",main="Duplicate Ratio",xlim=c(0,1),xlab="Duplicate Ratio",ylab="Frequency")
  hist(matrix[2:nrow(matrix),7]/matrix[2:nrow(matrix),2], border=NA, col="blue",main="Mitochondria Ratio",xlim=c(0,1),xlab="Mitochondria Ratio",ylab="Frequency")
  dev.off()
}

FripPlot_10x <- function(matrix, prefix, reads_cutoff = 1000, frip_cutoff = 0.05)
{
  png(paste0(prefix,"_frip.png"),width=6,height=6, res = 300, units = "in")
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
library(ggplot2)
library(cowplot)
library(grid)
library(pheatmap)

PipelineSeurat = function(countMat, proj, min.c = 10, min.p = 200, max.p = 100000, nfeatures = 25000, org="GRCh38",
  do.scale = T, normalization.method = "LogNormalize", dims.use = 1:15, res = 0.6, diff.p = F, diff.co = 1E-5)
{
  start_time <- Sys.time()  
  #=========QC========
  message("Check gene and cell coverage ...")
  nGene = apply(countMat, 2, function(x) length(x[x>0]))
  nCell = apply(countMat, 1, function(x) length(x[x>0]))
  png(paste0(proj,"_Coverage.png"),width=8,height=4.5, res = 300, units = "in")
  par(mfrow=c(1,2))
  plot(1:ncol(countMat),sort(nGene),pch=16,col="blue",ylab="Number of Peaks Detected",xlab="Cells",main="Cell Filter")
  abline(h=min.p,lwd=2,lty=2);text(ncol(countMat)/2,min.p+max(nGene)*0.05,paste0("n = ",min.p));legend("topleft",paste0("ave Peaks = ",round(mean(nGene))),box.lty=0)
  plot(1:nrow(countMat),sort(nCell),pch=16,col="blue",ylab="Number of Cells Detected",xlab="Peaks",main="Peak Filter")
  abline(h=min.c,lwd=2,lty=2);text(nrow(countMat)/2,min.c+max(nCell)*0.05,paste0("n = ",min.c));legend("topleft",paste0("ave Cells = ",round(mean(nCell))),box.lty=0)
  dev.off()
  SeuratObj = CreateSeuratObject(countMat, project = proj, min.cells = min.c, min.features = min.p)
  
  #=========Filter========
  message("Filter cells and find variable peaks ...")  
  SeuratObj <- subset(SeuratObj, subset.names = c("nGene"), low.thresholds = min.p, high.thresholds = max.p) 
  SeuratObj <- NormalizeData(SeuratObj, normalization.method = normalization.method, scale.factor = 10000)
  SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = nfeatures)
  if(do.scale) SeuratObj <- ScaleData(object = SeuratObj)
  
  #=========PCA===========
  message("PCA analysis ...")
  SeuratObj <- RunPCA(object = SeuratObj, features = VariableFeatures(SeuratObj))
  # pdf(file.path(paste0(proj,"_PCElbowPlot.pdf")), width = 5, height = 4)
  p1 = ElbowPlot(object = SeuratObj)
  ggsave(file.path(paste0(proj,"_PCElbowPlot.png")), p1,  width=5, height=4)
  # dev.off()
  
  #=========tSNE===========
  message("UMAP analysis ...")
  SeuratObj <- RunUMAP(object = SeuratObj, reduction = "pca", dims = dims.use)
  SeuratObj <- FindNeighbors(object = SeuratObj, reduction = "pca", dims = dims.use)
  SeuratObj <- FindClusters(object = SeuratObj, resolution = res)
  # pdf(file.path(paste0(proj, "_UMAP_cluster.pdf")), width=5, height=4)
  p2 = DimPlot(object = SeuratObj, pt.size = 0.5, label = TRUE)
  ggsave(file.path(paste0(proj, "_UMAP_cluster.png")), p2, width=5, height=4)
  # dev.off()

  #=========identify marker===========
  message("Find cluster specific peaks ...")
  cluster.specific.peaks = NULL
  if (diff.p){
  cluster.specific.peaks = FindAllMarkers(object = SeuratObj, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
  cluster.specific.peaks = cluster.specific.peaks[cluster.specific.peaks$p_val_adj<diff.co, ]
  saveRDS(cluster.specific.peaks, file.path(paste0(proj, "_ClusterSpecificPeaks.rds")))}
  write.table(cluster.specific.peaks, file.path(paste0(proj, "_ClusterSpecificPeaks.txt")), col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
  
  saveRDS(SeuratObj, file.path(paste0(proj, "_peak_SeuratObj.rds")))
  end_time <- Sys.time()
  end_time - start_time
  return(SeuratObj)
}

Incorparate <- function(Seurat1, Seurat2, proj, type1 = "ATAC", type2 = "RNA", pcs.use = 30, dims.use = 1:20, res = 0.6, anchor = 2000)
{
  SeuratObj1 = readRDS(Seurat1)
  SeuratObj2 = readRDS(Seurat2)
  SeuratObj1$batch <- type1
  SeuratObj2$batch <- type2
 
  Seurat.anchors <- FindIntegrationAnchors(object.list = list(SeuratObj1, SeuratObj2), dims = dims.use, anchor.features = anchor)
  Seurat.combined <- IntegrateData(anchorset = Seurat.anchors, dims = dims.use)
  DefaultAssay(Seurat.combined) <- "integrated"
  Seurat.combined <- ScaleData(Seurat.combined, verbose = FALSE)
  Seurat.combined <- RunPCA(Seurat.combined, npcs = pcs.use, verbose = FALSE)
  Seurat.combined <- RunUMAP(Seurat.combined, reduction = "pca", dims = dims.use)
  Seurat.combined <- FindNeighbors(Seurat.combined, reduction = "pca", dims = dims.use)
  Seurat.combined <- FindClusters(Seurat.combined, resolution = res)
 
  Seurat.combined@project.name <- proj
  p1 <- DimPlot(Seurat.combined, reduction = "umap", group.by = "batch")
  ggsave(file.path(paste0("UMAP_origIdent_", Seurat.combined@project.name, "_source.png")), p1, width=5, height=4)
  p2 <- DimPlot(Seurat.combined, reduction = "umap", label = TRUE)
  ggsave(file.path(paste0("UMAP_origIdent_", Seurat.combined@project.name, "_all.png")), p2, width=5, height=4)
  p3 <- DimPlot(Seurat.combined, reduction = "umap", group.by = "assign.ident", cells = colnames(SeuratObj2), label = TRUE, label.size = 2)
  ggsave(file.path(paste0("UMAP_origIdent_", Seurat.combined@project.name, "_RNAonly.png")), p3, width=6, height=4)
  p4 <- DimPlot(Seurat.combined, reduction = "umap", group.by = "RNA_snn_res.0.6", cells = colnames(SeuratObj1), label = TRUE, label.size = 2)
  ggsave(file.path(paste0("UMAP_origIdent_", Seurat.combined@project.name, "_ATAConly.png")), p4, width=5, height=4)
  
  saveRDS(Seurat.combined, file.path(paste0(Seurat.combined@project.name, "_SeuratObj.rds")))
  return(Seurat.combined)
}

source(paste0(Rdir,"scRNAseq_function.R"))
AnnotateRP <- function(SeuratObj, RPmat, proj){
  SeuratObjRP <- CreateSeuratObject(RPmat, project = proj, min.cells = 3, min.features = 200)
  SeuratObjRP <- subset(SeuratObjRP, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)  

  Idents(object = SeuratObjRP) <- Idents(object = SeuratObj)[colnames(x = SeuratObjRP)]
  message("Find marker genes and annotation ...")
  cluster.markers <- NULL
  cluster.markers <- FindAllMarkers(object = SeuratObjRP, only.pos = TRUE, min.pct = 0.1)

  saveRDS(cluster.markers, paste0(SeuratObj@project.name, "_ClusterSpecificGenesRP.rds"))

  cluster.markers <- cluster.markers[cluster.markers$p_val_adj<0.00001, ]

  current.cluster.ids = as.integer(levels(cluster.markers$cluster))
  new.cluster.ids = AssignCellTypeSeurat(cluster.markers)
  SeuratObj@meta.data$assign.ident = Idents(SeuratObj)[rownames(SeuratObj@meta.data)]
  SeuratObj@meta.data$assign.ident = plyr::mapvalues(x = SeuratObj@meta.data$assign.ident,
                                                     from = current.cluster.ids, to = new.cluster.ids)
  saveRDS(SeuratObj, file.path(paste0(SeuratObj@project.name, "_SeuratObj.rds")))
  # png(paste0("UMAP_assignIdent_", SeuratObj@project.name, "_RPannotated.png"),res=300, width=6.1, height=4, units = "in")
  p1 = DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2, group.by = "assign.ident", label.size = 3)
  ggsave(paste0("UMAP_assignIdent_", SeuratObj@project.name, "_RPannotated.png"),p1, width=5.8, height=4)
  # dev.off() 
}

multiplot <- function(x,y)
{
  viewport(layout.pos.row = x, layout.pos.col = y)
}

TF_UMAP <- function(SeuratObj, RPmat){
  cluster_list = unique(as.character(Idents(SeuratObj)))
  for(i in cluster_list){
    giggle_tf = read.table(paste0("SpetificTF/cluster_",i,"/giggle_res_tfs.txt"), 
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    tf_list = giggle_tf$factor
    tf_sample = paste(giggle_tf$factor, giggle_tf$sample_id, sep = "_")
    names(tf_sample) = tf_list
    for(tf in tf_list){
      tf_RP = t(RPmat[tf,colnames(x = SeuratObj)])
      umap_RP = SeuratObj@reductions$umap@cell.embeddings
      umap_RP = merge(umap_RP, tf_RP, by.x = 0, by.y = 0)
      row.names(umap_RP) = umap_RP[,1]
      umap_RP = umap_RP[,-1]
      colnames(umap_RP)[3] = "TF"
      
      tg_list = read.table(paste0("TF/giggle_res/", proj, "_strap_peak_ClusterSpecificPeaks/cluster_",i,"/",tf_sample[tf],"_target_genes_top500.txt"), 
                           header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
      tg_RP = RPmat[tg_list,colnames(x = SeuratObj)]
      tg_RP = t(na.omit(tg_RP))
      tg_RP_avg = apply(tg_RP, 1, mean)
      umap_RP[,"TG_avg"] = tg_RP_avg[row.names(umap_RP)]
      plotdir = paste0("SpetificTF/cluster_",i)
      png(file.path(plotdir,paste0("UMAP_", SeuratObj@project.name, tf, ".png")), width=10, height=4, units = "in", res = 300)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(1,2)))
      p1 = ggplot(umap_RP, aes(x = UMAP_1, y = UMAP_2, colour = TF)) + 
        geom_point(aes(colour = TF), size = 0.2) + 
        scale_color_gradient(name = tf,low="lightgrey", high="darkblue")
      p2 = ggplot(umap_RP, aes(x = UMAP_1, y = UMAP_2, colour = TG_avg)) + 
        geom_point(aes(colour = TG_avg), size = 0.2) + 
        scale_color_gradient(name = "Target",low="lightgrey", high="darkblue")
      print(p1, vp = multiplot(1,1))
      print(p2, vp = multiplot(1,2))
      dev.off()
    }
  }
}


#' Vlnplot to show top TF'
tf_cluster = readLines(paste0(srcdir, "annotations/human_HOCOMOCOv11_full_TF_cluster.txt"))
tf_cluster_list = lapply(tf_cluster, function(x){
  return(unlist(strsplit(x, "\t")))
})
names(tf_cluster_list) = sapply(tf_cluster_list, function(x){
  return(x[1])
})
tf_cluster_list = lapply(tf_cluster_list, function(x){
  return(x[-1])
})

TF_Vlnplot <- function(SeuratObj, RPmat){
  cluster_list = unlist(strsplit(list.files("SpetificTF"), "_"))[seq(2, 2*length(list.files("SpetificTF")), 2)]
  for(i in cluster_list){
    giggle_tf = read.table(paste0("SpetificTF/cluster_",i,"/giggle_res_tfs.txt"), 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    tf_list = giggle_tf$factor

    tf_expand_list = lapply(tf_list, function(x){
      if(is.null(tf_cluster_list[[x]])){
        return(c(x)) 
      }else{
        return(tf_cluster_list[[x]])
      }
    })

    cluster_cell = names(Idents(SeuratObj))[Idents(SeuratObj) == i]
    tf_expand_maxRP = sapply(tf_expand_list, function(x){
      tf_RP = RPmat[x, cluster_cell]
      tf_RP = na.omit(tf_RP)
      if(nrow(tf_RP)>0){
        tf_RP_mean = apply(tf_RP, 1, mean)
        return(names(tf_RP_mean)[which(tf_RP_mean == max(tf_RP_mean))])
      }else{
        return(NULL)
      }
    })
    tf_expand_maxRP = unlist(tf_expand_maxRP)

    j = 0
    
    plotdir = paste0("SpetificTF/cluster_",i)
    png(file.path(plotdir, paste0(SeuratObj@project.name, "_vlnplot_cluster",i,".png")), width=11.5, height=5, res = 300, units = "in")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2,5)))
    for(tf in tf_expand_maxRP){
      tf_RP = t(RPmat[tf,colnames(x = SeuratObj)])
      RP_cluster = Idents(SeuratObj)
      RP_cluster = merge(RP_cluster, tf_RP, by.x = 0, by.y = 0)
      row.names(RP_cluster) = RP_cluster[,1]
      RP_cluster = RP_cluster[,-1]
      colnames(RP_cluster) = c("Cluster","TF")
      RP_cluster$Cluster = as.factor(RP_cluster$Cluster)
      
      p = ggplot(RP_cluster, aes(Cluster, TF)) + 
        geom_violin(aes(fill = Cluster), show.legend = FALSE, scale = "width") + 
        labs(y = "RP score", title = tf) +
        theme(axis.text.x = element_text(size = 6), axis.title = element_text(size = 12), axis.text.y = element_text(size = 10))
      print(p, vp = multiplot(as.integer(j/5)+1, j%%5+1))
      j = j + 1
    }
    dev.off()
  }
}


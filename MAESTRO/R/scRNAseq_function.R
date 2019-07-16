#' Preprocess and analysis scRNA-seq data
#' Seurat 3.0.0
#' Function declare
#' @Count2TPM
#' @Count2FPKM
#' @Ensembl2Symbol
#' @PipelineSeurat
#' @PipelinePagoda
#' @AssignCellTypeSeurat
#' @AssignCellTypePagoda
#' @RemoveBatchCCA
#' @LineageDiffusionmap
#' @LineageRNAvelocity
#' @TriangleHeatmap
#' @PipelinescMCA
#' @plotMCA
#' @PipelineRCA
#' @PipelineSSCC
#' @AssignCellTypeSSCC

library(Seurat)
library(dplyr)
library(biomaRt)
library(MAST)
library(ggplot2)
library(Matrix)

srcdir = "../"
Rdir = "./"

hg38_genome = read.delim(paste0(srcdir, "annotations/GRCh38_mart_export.txt"), check.names = FALSE)
mm10_genome = read.delim(paste0(srcdir, "annotations/GRCm38_mart_export.txt"), check.names = FALSE)
Count2TPM <- function(countMat, idType = "ENSEMBL", organism="GRCh38")
{
  if(organism=="GRCh38")
    ensembl = hg38_genome
  else
    ensembl = mm10_genome
  ensembl$Length <- abs(ensembl$`Gene end (bp)` - ensembl$`Gene start (bp)`)
  if(toupper(idType) == "ENSEMBL"){
    len <- ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Length"]
    count_rowname = ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Gene name"]
    count_index = which(!duplicated(count_rowname) & !is.na(count_rowname))
    countMat = countMat[count_index,]
    rownames(countMat) = count_rowname[count_index]
    len = len[count_index]
  } else if(toupper(idType) == "SYMBOL")
    len <- ensembl[match(rownames(countMat),ensembl$`Gene name`), "Length"]
  else
    stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")

  na_idx = which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
    countMat = countMat[!is.na(len),]
    len = len[!is.na(len)]
  }
  tmp <- countMat / len
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
  TPM = TPM[!duplicated(rownames(TPM)),]
  return(TPM)
}

Count2FPKM <- function(countMat, idType = "Ensembl", organism="GRCh38")
{
  if(organism=="GRCh38")
    ensembl = hg38_genome
  else
    ensembl = mm10_genome
  ensembl$Length <- abs(ensembl$`Gene end (bp)` - ensembl$`Gene start (bp)`)
  if(toupper(idType) == "ENSEMBL"){
    len <- ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Length"]
    count_rowname = ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Gene name"]
    count_index = which(!duplicated(count_rowname) & !is.na(count_rowname))
    countMat = countMat[count_index,]
    rownames(countMat) = count_rowname[count_index]
    len = len[count_index]
} else if(toupper(idType) == "SYMBOL")
  len <- ensembl[match(rownames(countMat),ensembl$`Gene name`), "Length"]
else
  stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")

na_idx = which(is.na(len))
if(length(na_idx)>0){
  warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
  countMat = countMat[!is.na(len),]
  len = len[!is.na(len)]
}
tmp <- countMat / len
FPKM <- 1e6 * t(t(tmp) / colSums(countMat))
FPKM = FPKM[!duplicated(rownames(FPKM)),]
return(FPKM)
}

Ensembl2Symbol <- function(countMat, organism="GRCh38")
{
  if(organism=="GRCh38")
    ensembl = hg38_genome
  else
    ensembl = mm10_genome
  count_rowname = ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Gene name"]
  count_index = which(!duplicated(count_rowname) & !is.na(count_rowname))
  countMat = countMat[count_index,]
  rownames(countMat) = count_rowname[count_index]
  na_idx = which(is.na(count_rowname))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which symbol is not available !"))
  }
  return(countMat)
}

PipelineSeurat <- function(tpmMat, proj, min.c = 3, min.g = 200, max.g = 20000, nfeatures = 2000, organism="GRCh38",
  normalization.method = NULL, do.scale = TRUE, do.center = TRUE, vars.to.regress = c("nCount_RNA","percent.mito"),
  dims.use = 1:15, marker.use = "human_immune_CIBERSORT", res = 0.6, orig.ident = NULL)
{
  #=========QC========
  message("Check gene and cell coverage ...")
  nGene = apply(tpmMat, 2, function(x) length(x[x>0]))
  nCell = apply(tpmMat, 1, function(x) length(x[x>0]))
  png(paste0(proj,"_QC_Coverage.png"), res=300, width=8,height=4.5, units = "in")
  par(mfrow=c(1,2))
  plot(1:ncol(tpmMat),sort(nGene),pch=16,col="blue",ylab="Number of Genes Expressed",xlab="Cells",main="Cell Filter")
  abline(h=min.g,lwd=2,lty=2);text(ncol(tpmMat)/2,min.g+max(nGene)*0.05,paste0("n = ",min.g));legend("topleft",paste0("ave G = ",round(mean(nGene))),box.lty=0)
  plot(1:nrow(tpmMat),sort(nCell),pch=16,col="blue",ylab="Number of Cells Expressed",xlab="Genes",main="Gene Filter")
  abline(h=min.c,lwd=2,lty=2);text(nrow(tpmMat)/2,min.c+max(nCell)*0.05,paste0("n = ",min.c));legend("topleft",paste0("ave C = ",round(mean(nCell))),box.lty=0)
  dev.off()
  
  message("Check the mitochondria and spike-in percentage ...")
  SeuratObj <- CreateSeuratObject(tpmMat, project = proj, min.cells = min.c, min.features = min.g)
  if(!is.null(orig.ident)) SeuratObj$orig.ident <- orig.ident
  if(organism=="GRCh38"){
     mito.genes <- grep("^MT-", rownames(GetAssayData(object = SeuratObj)), value = TRUE)
     ercc.genes <- grep("^ERCC", rownames(GetAssayData(object = SeuratObj)), value = TRUE)}
  else{
     mito.genes <- grep("^mt-", rownames(GetAssayData(object = SeuratObj)), value = TRUE)
     ercc.genes <- grep("^ercc", rownames(GetAssayData(object = SeuratObj)), value = TRUE)}   
  percent.mito <- colSums(as.matrix(GetAssayData(object = SeuratObj)[mito.genes, ]))/colSums(as.matrix(GetAssayData(object = SeuratObj)))
  percent.ercc <- colSums(as.matrix(GetAssayData(object = SeuratObj)[ercc.genes, ]))/colSums(as.matrix(GetAssayData(object = SeuratObj)))
  SeuratObj$percent.mito <- percent.mito
  SeuratObj$percent.ercc <- percent.ercc
  # png(paste0(proj, "_QC_Spikein.png"), res=300, width=6, height=4.5, units = "in")
  p1 = VlnPlot(SeuratObj, c("percent.mito","percent.ercc"), ncol = 2)
  ggsave(paste0(proj, "_QC_Spikein.png"), p1,  width=6, height=4.5)
  # dev.off()
  
  #=========Filter========
  message("Filter cells and find variable genes ...")  
  SeuratObj <- subset(SeuratObj, subset.name = "percent.mito", high.threshold = 0.05)
  SeuratObj <- subset(SeuratObj, subset.name = "percent.ercc", high.threshold = 0.05)
  SeuratObj <- subset(SeuratObj, subset.names = "nGene", low.thresholds = min.g, high.thresholds = max.g) 
  
  SeuratObj <- NormalizeData(object = SeuratObj, normalization.method = normalization.method, scale.factor = 10000)
  SeuratObj <- FindVariableFeatures(object = SeuratObj, selection.method = "vst", nfeatures = nfeatures)
  if(do.scale) SeuratObj <- ScaleData(object = SeuratObj, vars.to.regress = vars.to.regress, do.center = do.center)
  
  #=========PCA===========
  message("PCA analysis ...")
  SeuratObj <- RunPCA(object = SeuratObj, features = VariableFeatures(SeuratObj))
  # png(file.path(paste0(proj,"_PCElbowPlot.png")), res=300, width=5, height=4, units = "in")
  p2 = ElbowPlot(object = SeuratObj)
  ggsave(file.path(paste0(proj,"_PCElbowPlot.png")), p2,  width=5, height=4)
  # dev.off()
  
  #=========tSNE===========
  message("UMAP analysis ...")
  SeuratObj <- RunUMAP(object = SeuratObj, reduction = "pca", dims = dims.use)
  SeuratObj <- FindNeighbors(object = SeuratObj, reduction = "pca", dims = dims.use)
  SeuratObj <- FindClusters(object = SeuratObj, resolution = res)
  # png(file.path(paste0("UMAP_origIdent_", SeuratObj@project.name, "_cluster.png")), res=300, width=5, height=4, units = "in")
  p3 = DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2)
  ggsave(file.path(paste0("UMAP_origIdent_", SeuratObj@project.name, "_cluster.png")), p3,  width=5, height=4)
  # dev.off()
  # png(file.path(paste0("UMAP_origIdent_", SeuratObj@project.name, "_primary.png")), res=300, width=6, height=4, units = "in")
  p4 = DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2, group.by = "orig.ident", label.size = 3)
  ggsave(file.path(paste0("UMAP_origIdent_", SeuratObj@project.name, "_primary.png")), p4,  width=6, height=4)
  # dev.off()
  saveRDS(SeuratObj, file.path(paste0(proj, "_SeuratObj.rds")))

  #=========identify marker===========
  message("Find marker genes and annotation ...")
  cluster.markers <- NULL
  cluster.markers <- FindAllMarkers(object = SeuratObj, only.pos = TRUE, min.pct = 0.1)
  cluster.markers <- cluster.markers[cluster.markers$p_val_adj<0.000001, ]

  current.cluster.ids = as.integer(levels(cluster.markers$cluster))
  new.cluster.ids = AssignCellTypeSeurat(cluster.markers, marker.use)
  SeuratObj@meta.data$assign.ident = Idents(SeuratObj)[rownames(SeuratObj@meta.data)]
  SeuratObj@meta.data$assign.ident = plyr::mapvalues(x = SeuratObj@meta.data$assign.ident,
                       from = current.cluster.ids, to = new.cluster.ids)
  # png(file.path(paste0("UMAP_assignIdent_", SeuratObj@project.name, "_annotated.png")),res=300, width=6, height=4, units = "in")
  p5 = DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2, group.by = "assign.ident", label.size = 3)
  ggsave(file.path(paste0("UMAP_assignIdent_", SeuratObj@project.name, "_annotated.png")), p5,  width=6, height=4)
  # dev.off()  
  
  saveRDS(cluster.markers, file.path(paste0(proj, "_DiffMarkers.rds")))
  saveRDS(SeuratObj, file.path(paste0(proj, "_SeuratObj.rds")))
  return(list(SeuratObj = SeuratObj, markers = cluster.markers))
}

PipelinePagoda <- function(tpmMat, proj, min.c = 3, min.g = 200, nPCs = 100, nKs = 30)
{
  library(pagoda2)

  #=========QC========
  message("Check gene and cell coverage ...")
  nGene = apply(tpmMat, 2, function(x) length(x[x>0]))
  nCell = apply(tpmMat, 1, function(x) length(x[x>0]))
  pdf(paste0(proj,"_QC_Coverage.pdf"),width=8,height=4.5)
  par(mfrow=c(1,2))
  plot(1:ncol(tpmMat),sort(nGene),pch=16,col="blue",ylab="Number of Genes Expressed",xlab="Cells",main="Cell Filter")
  abline(h=min.g,lwd=2,lty=2);text(ncol(tpmMat)/2,min.g+max(nGene)*0.05,paste0("n = ",min.g));legend("topleft",paste0("ave G = ",round(mean(nGene))),box.lty=0)
  plot(1:nrow(tpmMat),sort(nCell),pch=16,col="blue",ylab="Number of Cells Expressed",xlab="Genes",main="Gene Filter")
  abline(h=min.c,lwd=2,lty=2);text(nrow(tpmMat)/2,min.c+max(nCell)*0.05,paste0("n = ",min.c));legend("topleft",paste0("ave C = ",round(mean(nCell))),box.lty=0)
  dev.off()

  #=========PCA & tSNE=======
  message("Filter cells, variance shrinkage, PCA and tSNE analysis ...")
  tpmMat <- tpmMat[!duplicated(rownames(tpmMat)),]
  r <- Pagoda2$new(tpmMat,modelType='plain',trim=10,log.scale=T)
  r$adjustVariance(plot=T,do.par=T,gam.k=10)
  r$calculatePcaReduction(nPcs=nPCs,n.odgenes=3e3,maxit=1000)
  r$makeKnnGraph(k=nKs,type='PCA',center=T,distance='cosine');
  require(igraph)
  r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
  r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
  png(file.path(paste0("tSNE_origIdent_", proj, "_cluster.png")),res=300, width=4, height=4, units = "in")
  r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main=proj)
  dev.off()
 
  #=========identify marker===========
  message("Find marker genes ...")
  deSets <- get.de.geneset(r, groups = r$clusters$PCA[[1]], prefix = 'de_')
  saveRDS(deSets, file.path(paste0(proj, "_Pagoda_DiffMarkers.rds")))
  
  CellType <- AssignCellTypePagoda(deSets)
  r$clusters$PCA$annotation <- r$clusters$PCA$multilevel
  r$clusters$PCA$multilevel <- plyr::mapvalues(x = r$clusters$PCA$multilevel,
                     from = levels(r$clusters$PCA$multilevel), to = CellType)
  png(file.path(paste0("tSNE_origIdent_", proj, "_annotated.png")),res=300, width=4, height=4, units = "in")
  r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main=proj)
  dev.off() 
  
  saveRDS(r, file.path(paste0(proj, "_PagodaObj.rds")))
  return(list(PagodaObj=r, markers = deSets))
}

source(paste0(Rdir,"scRNAseq_annotation.R"))
AssignCellTypeSeurat <- function(cluster.markers, marker.use = "human_immune_CIBERSORT", organism="GRCh38"){
# celltype.markers include human_immune_TCIA, markers.pubmed, "human_immune_CIBERSORT", mouse_brain (from split-seq), markers.brain.adult(from science)
  if(marker.use == "human_immune_CIBERSORT" | marker.use == "human_immune_TCIA" | marker.use == "human_immune_simple"){
    if(organism != "GRCh38"){
      marker_mouse = cluster.markers$gene
      ensembl = useMart("ensembl")
      human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      marker_human = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = marker_mouse ,mart = mouse, attributesL = "hgnc_symbol", martL = human, uniqueRows=T)
      cluster.markers = merge(cluster.markers,marker_human,by.x = "gene", by.y = "MGI.symbol")
      cluster.markers = cluster.markers[,-1]
      colnames(cluster.markers)[7] = "gene"
    }
  }
  else{
    if(organism == "GRCh38"){
      marker_human = cluster.markers$gene
      ensembl = useMart("ensembl")
      human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      marker_mouse = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = marker_human ,mart = human, attributesL = "mgi_symbol", martL = mouse, uniqueRows=T)
      cluster.markers = merge(cluster.markers,marker_mouse,by.x = "gene", by.y = "HGNC.symbol")
      cluster.markers = cluster.markers[,-1]
      colnames(cluster.markers)[7] = "gene"
    }
  }
  marker.use = get(marker.use)
  CT = sapply(as.integer(levels(cluster.markers$cluster)), function(x){
  idx = cluster.markers$cluster==x
  tmpFC = cluster.markers$avg_logFC[idx]
  names(tmpFC) = toupper(cluster.markers$gene[idx])
  score_cluster = sapply(marker.use, function(y){
    score = sum(tmpFC[y], na.rm = TRUE) / log2(length(y))
    if(!("PPBP" %in% y) & length(na.omit(tmpFC[y]))==1)
      score = 0
    return(score)
    })
    if(max(score_cluster, na.rm = TRUE)>0.1)
      CT = names(score_cluster)[which.max(score_cluster)]
    else
      CT = "Others"
  })
  return(CT)
}

AssignCellTypePagoda <- function(cluster.markers, marker.use = "human_immune_CIBERSORT"){
    marker.use = get(marker.use)
    CT = NULL
    for (i in cluster.markers){
    score_cluster = sapply(marker.use,function(x){
      score = length(intersect(x,i$genes))/length(x)
      return(score)})
    if(max(score_cluster, na.rm = TRUE)>0.05)
      CT <- c(CT,names(score_cluster)[which.max(score_cluster)])
    else
      CT <- c(CT, "Others")
      }
  return(CT)
}

AssignCellTypeMannually <- function(pbmc.markers, ident = 1, celltype.markers = "human_immune_CIBERSORT"){
    celltype.markers = get(celltype.markers)
    dx = pbmc.markers$cluster==ident
    tmpFC = pbmc.markers$avg_logFC[idx]
    names(tmpFC) = pbmc.markers$gene[idx]
    score_cluster = sapply(celltype.markers, function(y){
      score = sum(tmpFC[y], na.rm = TRUE) / log2(length(y))
      if(!("PPBP" %in% y) & length(na.omit(tmpFC[y]))==1)
        score = 0
      return(score)
    })
   return(score_cluster)
}

RemoveBatchCCA <- function(Seurat1, Seurat2, gene.co = 10000)
{
  SeuratObj1 <- readRDS(Seurat1)
  SeuratObj2 <- readRDS(Seurat2)
  SeuratObj1@meta.data$batch <- "Batch1"
  SeuratObj2@meta.data$batch <- "Batch2"
  
  g.1 <- head(rownames(SeuratObj1@hvg.info), gene.co)
  g.2 <- head(rownames(SeuratObj2@hvg.info), gene.co)
  genes.use <- unique(c(g.1, g.2))
  genes.use <- intersect(genes.use, rownames(SeuratObj1@scale.data))
  genes.use <- intersect(genes.use, rownames(SeuratObj2@scale.data))
  combined <- RunCCA(SeuratObj1, SeuratObj2, genes.use = genes.use, num.cc = 30)
  combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "batch", dims.align = 1:20) 
  combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T, check_duplicates = FALSE)
  combined <- FindClusters(combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)
  
  combined@meta.data$orig.ident <- gsub("\\..*","",colnames(combined@raw.data))
  pdf(file.path(paste0("tSNE_origIdent_", combined@project.name, "_CCA.pdf")), width=6, height=4)
  TSNEPlot(object = combined, do.label = TRUE, pt.size = 0.5, group.by = "orig.ident")
  dev.off() 
  
  cluster.markers <- NULL
  cluster.markers <- FindAllMarkers(object = combined, only.pos = TRUE, min.pct = 0.1)
  cluster.markers <- cluster.markers[cluster.markers$p_val_adj<0.01, ]
  
  current.cluster.ids = as.integer(levels(cluster.markers$cluster))
  new.cluster.ids = AssignCellTypeSeurat(cluster.markers, marker.use)
  combined@meta.data$assign.ident = combined@ident[rownames(combined@meta.data)]
  combined@meta.data$assign.ident = plyr::mapvalues(x = combined@meta.data$assign.ident,
                       from = current.cluster.ids, to = new.cluster.ids)
  pdf(file.path(paste0("tSNE_assignIdent_", combined@project.name, "_annotated.pdf")), width=6, height=4)
  TSNEPlot(object = combined, do.label = TRUE, pt.size = 0.5, group.by = "assign.ident")
  dev.off()      
  
  saveRDS(cluster.markers, file.path(paste0(proj, "_DiffMarkers.rds")))
  saveRDS(combined, file.path(paste0(proj, "_SeuratObj.rds")))
  return(list(SeuratObj = combined, markers = cluster.markers))
}

RemoveBatchMutiCCA <- function(object.list, proj, genes.use, num.ccs = 30, marker.use = "markers.CIBERSORT", res = 0.6)
{
  message("CCA analysis ...")
  for(i in 1:length(object.list)){
    object.list[[i]]@meta.data$batch <- names(object.list)[i]
  }
  
  Seurat_CCA@project.name = proj
  Seurat_CCA = RunMultiCCA(object.list = object.list, genes.use = genes.use, num.ccs = num.ccs)
  Seurat_CCA = AlignSubspace(Seurat_CCA, reduction.type = "cca", grouping.var = "batch", dims.align = 1:num.ccs)

  message("t-SNE analysis ...")
  Seurat_CCA = RunTSNE(Seurat_CCA, reduction.use = "cca.aligned", dims.use = 1:num.ccs, do.fast = T, check_duplicates = FALSE)
  Seurat_CCA = FindClusters(Seurat_CCA, reduction.type = "cca.aligned", resolution = res, dims.use = 1:num.ccs)
  
  #Seurat_CCA@meta.data$orig.ident <- gsub("\\..*","",colnames(Seurat_CCA@raw.data))
  png(file.path(paste0("tSNE_origIdent_", Seurat_CCA@project.name, "_cluster_afterCCA.png")), res=300, width=6, height=4, units = "in")
  TSNEPlot(object = Seurat_CCA, do.label = TRUE, pt.size = 0.2, group.by = paste0("res.",res))
  dev.off()
  png(file.path(paste0("tSNE_origIdent_", Seurat_CCA@project.name, "_primary_afterCCA.png")), res=300, width=6, height=4, units = "in")
  TSNEPlot(object = Seurat_CCA, do.label = TRUE, pt.size = 0.2, group.by = "batch")
  dev.off() 
  
  message("Find marker genes and annotation ...")
  cluster.markers <- NULL
  cluster.markers <- FindAllMarkers(object = Seurat_CCA, only.pos = TRUE, min.pct = 0.1)
  cluster.markers <- cluster.markers[cluster.markers$p_val_adj<0.01, ]
  
  current.cluster.ids = as.integer(levels(cluster.markers$cluster))
  new.cluster.ids = AssignCellTypeSeurat(cluster.markers, marker.use)
  Seurat_CCA@meta.data$assign.ident = Seurat_CCA@ident[rownames(Seurat_CCA@meta.data)]
  Seurat_CCA@meta.data$assign.ident = plyr::mapvalues(x = Seurat_CCA@meta.data$assign.ident,
                       from = current.cluster.ids, to = new.cluster.ids)
  png(file.path(paste0("tSNE_assignIdent_", Seurat_CCA@project.name, "_annotated_afterCCA.png")), res=300, width=6, height=4, units = "in")
  TSNEPlot(object = Seurat_CCA, do.label = TRUE, pt.size = 0.2, group.by = "assign.ident")
  dev.off()      
  
  saveRDS(cluster.markers, file.path(paste0(proj, "_DiffMarkersCCA.rds")))
  saveRDS(Seurat_CCA, file.path(paste0(proj, "_SeuratObjCCA.rds")))
  return(list(SeuratObj = Seurat_CCA, markers = cluster.markers))
}


TriangleHeatmap <- function(x1, x2, titN, colN, subN, l1, l2) {
  if (is.data.frame(x1)) {
    x1 = as.matrix(x1)
    x2 = as.matrix(x2)
  }
  xlim = c(0,nrow(x1))
  ylim = c(-1,ncol(x1)+1)
  layout(matrix(c(1,2,3),nrow=1,ncol=3,byrow=T), width = c(10,1,1))
  par(mar=c(2,5,1,1),mai=c(0.8,0.5,0.6,0.3),cex=0.75)
  plot(1:nrow(x1), 1:nrow(x1), main = "", xlim=xlim, ylim=ylim, type="n", xaxs = "r", yaxs="r", axes=F, xlab="", ylab="")
  ## change subtitle name to selection criteria
  title(titN, cex.main=1.5, font.main=4, cex.sub=1.2, font.sub=3, sub=subN)  
  axis(1, at=1:nrow(x1), labels=row.names(x1), las=2, tick=F, line=NA,outer=F, hadj=0) 
  axis(2, at=1:ncol(x1), labels=colN, las=1, tick=F, line=NA,outer=F, hadj=0.5)
  Dat <- c(t(x1))
  interval = cut(Dat, pretty(Dat, n=200), right=F) ## 200 interval
  colors = colorRampPalette(c("white","red"))(210)[as.numeric(interval)]
  for (j in 1:nrow(x1)) {
    for (i in 1:ncol(x1)) {
      rect(j-0.5,i-0.5,j+0.5,i+0.5)
      polygon(c(j,j,j+1)-0.5,c(i-1,i,i)+0.5,col=colors[i+ncol(x1)*(j-1)], border="white")    
    }
  }
  Dat <- c(t(x2))
  interval = cut(Dat, pretty(Dat, n=200), right=F) ## 200 interval
  colors = colorRampPalette(c("white","blue"))(210)[as.numeric(interval)]
  for (j in 1:nrow(x2)) {
    for (i in 1:ncol(x2)) {
      polygon(c(j,j+1,j+1)-0.5,c(i-1,i-1,i)+0.5,col=colors[i+ncol(x2)*(j-1)], border="white")    
    }
  }
  image(1, 0:100, matrix(data=0:100, nrow=1,ncol=100), col=colorRampPalette(c("white","red"))(200), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
  axis(side=2,seq(0,100,20),paste0(seq(0,100,20),"%"));box();mtext(side=2,l1,line=3,cex=0.65)
  image(1, 0:100, matrix(data=0:100, nrow=1,ncol=100), col=colorRampPalette(c("white","blue"))(200), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
  axis(side=2,seq(0,100,20),paste0(seq(0,100,20),"%"));box();mtext(side=2,l2,line=3,cex=0.65)
}

# PipelinescMCA
PipelinescMCA <- function(tpmMat, proj, numbers_plot = 3){
  library(scMCA)
  mca_result = scMCA(scdata = tpmMat, numbers_plot = 3)
  saveRDS(mca_result,file.path(paste0(proj, "_scMCAResult.rds")))
  pdf(paste0(proj, "_log_MCA_heatmap.pdf"),height = 13, width = 20)
  plotMCA(mca_result, show_tree = F)
  dev.off()
}

# plot heatmap in MCA
plotMCA <- function(mca_result,interactive_plot=F, numbers_plot=3, col_font_size = 1, row_font_size=8, show_col=T,show_bar=T, show_tree = T){
  data(ref.expr)
  cors <- mca_result$cors_matrix
  cors_index <- apply(cors,2,gettissue,numbers_plot)
  cors_index <- sort(unique(as.integer(cors_index)))
  data = cors[cors_index,]
  cors_out = reshape2::melt(data)[c(2,1,3)]
  colnames(cors_out)<- c("Cell","Cell type","Score")
  cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% top_n(n=numbers_plot,wt=Score))
  mca_result$scMCA_probility <- cors_out
  mca_result$top_cors <- numbers_plot
  height=dim(data)[1]*10+230
  tree_height = 0
  if(isTRUE(show_tree)){tree_height=50}
  
  p<-pheatmap(
    data,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D",
    fontsize_col = col_font_size,
    fontsize_row = row_font_size,
    color = colorRampPalette(c("grey", "white", "red"))(100),
    cellheight = 10,
    show_colnames = show_col,
    border_color = NA,
    height = height,
    legend = show_bar,
    treeheight_col = tree_height,
    treeheight_row = tree_height
  )
  if(isTRUE(interactive_plot)){
    
    inter_data<-data[rev(p$tree_row$order),][,p$tree_col$order]
    height= length(p$tree_row$order)*10+230
    plot_ly(x=colnames(inter_data),y=rownames(inter_data),z = inter_data, colors = colorRamp(c("grey", "white","red")),height=height, type = "heatmap", showscale=show_bar) %>% layout(autosize=T,  margin=list(l=0,r=230,b=180,t=20,pad=4),font=list(size=row_font_size),xaxis=list(showticklabels=show_col),yaxis=list(side="right"))
  }
  else{
    p
  }
}

# PipelineRCA
PipelineRCA <- function(tpmMat,proj){
  library(WGCNA)
  library(flashClust)
  library(gplots)
  library(preprocessCore)
  library(RCA)

  data_obj = dataConstruct(as.data.frame(as.matrix(tpmMat)))
  data_obj = geneFilt(obj_in = data_obj)
  data_obj = cellNormalize(data_obj, method = "scQ")
  data_obj = dataTransform(data_obj,"log10")
  data_obj = featureConstruct(data_obj,method = "SelfProjection")
  data_obj = cellClust(data_obj)
  RCAPlot(data_obj)
  saveRDS(data_obj,file.path(paste0(proj, "_RCAObj.rds")))
}

PipelineSSCC <- function(tpmMat,proj){
  library(SingleCellExperiment)
  library(sscClust)
  sce <- ssc.build(as.data.frame(as.matrix(exp.dat)), display.name = rownames(exp.dat))
  sce <- ssc.run(sce, subsampling=F, do.DE = TRUE,method.clust = "SNN",method.reduction = "pca", k.batch = 1:11)
  saveRDS(sce,file.path(paste0(proj, "_SSCC_sceObj.rds")))

  CellType = AssignCellTypeSSCC(sce@metadata$ssc$de.res$L1C1$aov.out.sig, unique(colData(sce)@listData[[1]]),celltype.markers = "markers.CIBERSORT")
  colData(sce)@listData[[2]] = CellType[colData(sce)@listData[[1]]]
  names(colData(sce)@listData)[2] = "CellType"
  
  pdf("sscc_pbmc_all_pca_snn.pdf",width = 18,height = 7)
  ssc.plot.tsne(sce,columns = c("pca.SNN.kauto","CellType"),reduced.name = "pca.tsne")
  dev.off()
}

AssignCellTypeSSCC <- function(cluster.markers, cluster.all, celltype.markers = "markers.CIBERSORT"){
  celltype.markers = get(celltype.markers)
  cluster.marker.list = list()
  cluster.markers$cluster[cluster.markers$HSD.padj.min.diff>0] = unlist(strsplit(cluster.markers[cluster.markers$HSD.padj.min.diff>0,"HSD.padj.min.cmp"],
                                                                           split = "-"))[seq(1,2*length(cluster.markers[cluster.markers$HSD.padj.min.diff>0,"HSD.padj.min.cmp"]),2)]
  cluster.markers$cluster[cluster.markers$HSD.padj.min.diff<0] = unlist(strsplit(cluster.markers[cluster.markers$HSD.padj.min.diff<0,"HSD.padj.min.cmp"],
                                                                           split = "-"))[seq(2,2*length(cluster.markers[cluster.markers$HSD.padj.min.diff<0,"HSD.padj.min.cmp"]),2)]
  for(i in 1:length(cluster.all)){
    cluster.marker.list[[i]] = rownames(cluster.markers[cluster.markers$cluster == cluster.all[i],])
  }
  names(cluster.marker.list) = cluster.all
  CT = NULL
  for (i in cluster.marker.list){
      score_cluster = sapply(celltype.markers,function(x){
        score = length(intersect(x,i))/length(x)
        return(score)})
      if(max(score_cluster, na.rm = TRUE)>0.05)
        CT <- c(CT,names(score_cluster)[which.max(score_cluster)])
      else
        CT <- c(CT, "Others")
    }
    names(CT) = names(cluster.marker.list)
    return(CT)
  }

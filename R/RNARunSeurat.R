#' Clustering analysis for scRNA-seq data using Seurat
#'
#' Clustering analysis for scRNA-seq dataset using Seurat(version >=3.0.1). Including normalization, feature selection, dimension reduction, clustering and UMAP visualization. To run UMAP analysis, you must first install the umap-learn python package (e.g. via \code{pip install umap-learn}). 
#'
#' @docType methods
#' @name RNARunSeurat
#' @rdname RNARunSeurat
#'
#' @param inputMat Input unnormalized matrix, with genes as rows and cells as columns. Could be count matrix or TPM, FPKM matrix.
#' @param type Type of the input matrix. Default is "matrix". Set to "object" if the input is Seurat object.
#' @param project Output project name. Default is "MAESTRO.scRNA.Seurat".
#' @param orig.ident Orginal identity for the input cells. If supplied, should keep the same order with the column name of the gene x cell matrix.
#' @param min.c Minimum number of cells required for a gene. Will exclude the genes from input matrix if they only expressed in 
#' less than \code{min.c} cells. Default is 10. See \code{\link{CreateSeuratObject}} for details.
#' @param min.g Minimum number of genes required for a cell. Will exclude the cells from input matrix if less than \code{min.g}
#' genes are deteced in one cell. Default is 200. See \code{\link{CreateSeuratObject}} for details.
#' @param mito Whether or not to check and normalized accroding to the fraction of mitochondria reads and ercc spike-in reads.
#' If \code{mito} is True, those cells with mitochondria reads or ercc spike-in reads higher than 0.05 will be filtered.
#' Default is FALSE.
#' @param variable.genes Number of variable genes considered in the clustering analysis. Default is 2000. See \code{\link{FindVariableFeatures}} for details.
#' @param organism Organism for the dataset. Only support "GRCh38" and "GRCm38". Default is "GRCh38".
#' @param dims.use Number of dimensions used for UMAP analysis. Default is 1:15, use the first 15 PCs.
#' @param cluster.res Value of the clustering resolution parameter. Default is 0.6.
#' @param only.pos If seting true, only positive genes will be output. Default is False.
#' @param genes.test.use Denotes which test to use to identify genes. Default is "presto", a fast version of Wilcoxon Rank Sum test. See \code{\link{FindAllMarkersMAESTRO}} for details.
#' @param genes.cutoff Identify differential expressed genes with adjusted p.value less than \code{genes.cutoff} as cluster speficic genes
#' @param genes.pct Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently detected Default is 0.1
#' @param genes.logfc Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.2 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' for each cluster. Default cutoff is 1E-5.
#' @param runpca.args Extra arguments passed to \code{\link{RunPCA}}.
#' @param findneighbors.args Extra arguments passed to \code{\link{FindNeighbors}}.
#' @param findclusters.args Extra arguments passed to \code{\link{FindClusters}}.
#' @param \dots Extra arguments passed to \code{\link{RunUMAP}}.
#'
#' @author Chenfei Wang
#'
#' @return A list contains a RNA Seurat object and a data frame for the cluster specific genes. 
#'
#'
#' @examples
#' data(pbmc.RNA)
#' pbmc.RNA.res <- RNARunSeurat(inputMat = pbmc.RNA, project = "PBMC.scRNA.Seurat", min.c = 10, min.g = 100, dims.use = 1:15)
#' str(pbmc.RNA.res$RNA)
#' head(pbmc.RNA.res$genes)
#'
#' @importFrom Seurat CreateSeuratObject DimPlot ElbowPlot FindClusters FindNeighbors FindVariableFeatures GetAssayData NormalizeData RunPCA RunUMAP ScaleData SubsetData VlnPlot
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom Gmisc fastDoCall
#' @export

RNARunSeurat <- function(inputMat, type = "matrix", project = "MAESTRO.scRNA.Seurat", orig.ident = NULL, min.c = 10, min.g = 200, mito = FALSE, 
                         variable.genes = 2000, organism = "GRCh38", dims.use = 1:15, cluster.res = 0.6, only.pos = FALSE, genes.test.use = "presto", 
                         genes.cutoff = 1E-5, genes.pct = 0.1, genes.logfc = 0.25,
                         runpca.agrs = list(), findneighbors.args = list(), 
                         findclusters.args = list(), ...)
{
  if(type == "matrix"){SeuratObj <- CreateSeuratObject(inputMat, project = project, min.cells = min.c, min.features = min.g)}
  if(type == "object"){SeuratObj <- inputMat}

  #=========Mitochondria and Spike-in========  
  if(mito){
    message("Check the mitochondria and spike-in percentage ...")
    if(organism=="GRCh38"){
       mito.genes <- grep("^MT-", rownames(GetAssayData(object = SeuratObj)), value = TRUE)
       ercc.genes <- grep("^ERCC", rownames(GetAssayData(object = SeuratObj)), value = TRUE)}
    else{
       mito.genes <- grep("^mt-", rownames(GetAssayData(object = SeuratObj)), value = TRUE)
       ercc.genes <- grep("^ercc", rownames(GetAssayData(object = SeuratObj)), value = TRUE)}   
    percent.mito <- Matrix::colSums(GetAssayData(object = SeuratObj)[mito.genes, ])/Matrix::colSums(GetAssayData(object = SeuratObj))
    percent.ercc <- Matrix::colSums(GetAssayData(object = SeuratObj)[ercc.genes, ])/Matrix::colSums(GetAssayData(object = SeuratObj))
    SeuratObj$percent.mito <- percent.mito
    SeuratObj$percent.ercc <- percent.ercc
    p1 = VlnPlot(SeuratObj, c("percent.mito","percent.ercc"), ncol = 2)
    ggsave(paste0(SeuratObj@project.name, ".spikein.png"), p1,  width=6, height=4.5)
    
    SeuratObj <- subset(x = SeuratObj, subset = percent.mito < 0.05)
    SeuratObj <- subset(x = SeuratObj, subset = percent.ercc < 0.05)
    vars.to.regress = c("nCount_RNA","percent.mito","percent.ercc")}
  else{
    vars.to.regress = "nCount_RNA"}
  
  #=========Variable genes========
  message("Normalization and identify variable genes ...")  
  SeuratObj <- NormalizeData(object = SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  SeuratObj <- FindVariableFeatures(object = SeuratObj, selection.method = "vst", nfeatures = variable.genes)
  SeuratObj <- ScaleData(object = SeuratObj, vars.to.regress = vars.to.regress)
  
  #=========PCA===========
  message("PCA analysis ...")
  SeuratObj <- fastDoCall("RunPCA", c(object = SeuratObj, features = VariableFeatures(SeuratObj), runpca.agrs))
  # SeuratObj <- RunPCA(object = SeuratObj, features = VariableFeatures(SeuratObj))
  p2 = ElbowPlot(object = SeuratObj, ndims = SeuratObj@commands$RunPCA.RNA@params$npcs)
  ggsave(file.path(paste0(SeuratObj@project.name, "_PCElbowPlot.png")), p2,  width=10, height=4)
  
  #=========UMAP===========
  message("UMAP analysis ...")
  SeuratObj <- RunUMAP(object = SeuratObj, reduction = "pca", dims = dims.use, ...)
  SeuratObj <- fastDoCall("FindNeighbors", c(object = SeuratObj, reduction = "pca", dims = dims.use, findneighbors.args))
  SeuratObj <- fastDoCall("FindClusters", c(object = SeuratObj, resolution = cluster.res, findclusters.args))
  # SeuratObj <- FindNeighbors(object = SeuratObj, reduction = "pca", dims = dims.use)
  # SeuratObj <- FindClusters(object = SeuratObj, resolution = cluster.res)
  p3 = DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2)
  ggsave(file.path(paste0(SeuratObj@project.name, "_cluster.png")), p3,  width=5, height=4)
  if(!is.null(orig.ident)){
    SeuratObj$orig.ident <- orig.ident
    p4 = DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2, group.by = "orig.ident", label.size = 3)
    ggsave(file.path(paste0(SeuratObj@project.name, "_original.png")), p4,  width=5, height=4)}

  #=========DE analysis===========
  message("Identify cluster specific genes ...")
  cluster.genes <- FindAllMarkersMAESTRO(object = SeuratObj, min.pct = genes.pct, logfc.threshold = genes.logfc, test.use = genes.test.use, only.pos = only.pos)
  cluster.genes <- cluster.genes[cluster.genes$p_val_adj<genes.cutoff, ]
  write.table(cluster.genes, paste0(SeuratObj@project.name, "_DiffGenes.tsv"), quote = F, sep = "\t")

  return(list(RNA=SeuratObj, genes=cluster.genes))
}

#' Evaluate the transcription factor enrichment in scATAC-seq using chromVAR
#'
#' Evaluate the transcription factor enrichment in scATAC-seq using chromVAR (version >= 1.1.1).
#'
#' @docType methods
#' @name ATACRunchromVAR
#' @rdname ATACRunchromVAR
#'
#' @param inputMat Input unnormalized count matrix, with peaks as rows and cells as columns. Rownames should be like "chromosome_peakstart_peakend",
#' for example "chr10_100020591_100020841".
#' @param project Output project name for the chromVAR object. Default is "MAESTRO.scATAC.chromVAR".
#' @param min.c Minimum number of cells required for a peak. Will exclude the peaks from input matrix if they only identified in 
#' less than \code{min.c} cells. Default is 50. 
#' @param min.p Minimum number of peaks required for a cell. Will exclude the cells from input matrix if less than \code{min.p}
#' peaks are deteced in one cell. Default is 500.
#' @param organism Organism for the dataset. Only support "GRCh38" and "GRCm38". Default is "GRCh38".
#'
#' @author Chenfei Wang
#'
#' @return A list contains a chromVAR deviation object, and a data frame with rows as TFs, columns as cells, and values as TF z-score in each cell. 
#'
#'
#' @examples
#' data(pbmc.ATAC)
#' pbmc.ATAC.res <- ATACRunchromVAR(inputMat = pbmc.ATAC, project = "PBMC.scATAC.chromVAR")
#' str(pbmc.ATAC.res$dev)
#' head(pbmc.ATAC.res$zscore[,1:10])
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @export

ATACRunchromVAR <- function(inputMat, project = "MAESTRO.scATAC.chromVAR", min.c = 50, min.p = 500, organism = "GRCh38")
{
  library(chromVAR)
  library(motifmatchr)
  library(chromVARmotifs)
  
  if (!requireNamespace("chromVAR", quietly = TRUE)) {
    stop("Package \"chromVAR\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("motifmatchr", quietly = TRUE)) {
    stop("Package \"motifmatchr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("chromVARmotifs", quietly = TRUE)) {
    stop("Package \"chromVARmotifs\" needed for this function to work. Please install it.",
         call. = FALSE)
  }  

  if(organism == "GRCh38") {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("Package \"BSgenome.Hsapiens.UCSC.hg38\" needed for this function to work. Please install it.",
           call. = FALSE)}
     data("human_pwms_v2")}
  if(organism == "GRCm38") {
    if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
      stop("Package \"BSgenome.Mmusculus.UCSC.mm10\" needed for this function to work. Please install it.",
           call. = FALSE)}
     data("mouse_pwms_v2")}
  
  #============ Peak filtering and normalization ============
  message("Peak filtering ...")
  inputMat <- as.matrix(inputMat)
  peaks <- data.frame(chr=unlist(strsplit(rownames(inputMat),'_'))[seq(1,nrow(inputMat)*3,3)], 
                 start=as.numeric(unlist(strsplit(rownames(inputMat),'_'))[seq(2,nrow(inputMat)*3,3)]),
                 end=as.numeric(unlist(strsplit(rownames(inputMat),'_'))[seq(3,nrow(inputMat)*3,3)]))  
  rownames(peaks) <- paste0(peaks[,1],'_',peaks[,2],'_',peaks[,3])
  colData <- DataFrame(celltype=colnames(inputMat), depth=apply(inputMat,2,sum))
  colData <- colData[which(colData[,2]>min.p),]
  rowData <- DataFrame(peaks=rownames(inputMat), depth=apply(inputMat[,colData[,1]],1,sum))
  rowData <- rowData[which(rowData[,2]>min.c),]
  finalPeak <- as.matrix(peaks[rowData[,1],])
  finalCount <- as.matrix(inputMat[rowData[,1],colData[,1]])
  rowRanges = GRanges(as.character(finalPeak[,1]),IRanges(as.integer(finalPeak[,2]),as.integer(finalPeak[,3])),strand="*",score=as.integer(5),qval=1,name=rownames(finalPeak))
  chromVAR_counts  = SummarizedExperiment(assays=SimpleList(counts=finalCount),rowRanges=rowRanges, colData=colData)
  
  #============ chromVAR analysis ============
  message("chromVAR analysis ...")
  if(organism == "GRCh38") {
     chromVAR_counts = addGCBias(chromVAR_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
     motif_ix = matchMotifs(human_pwms_v2, chromVAR_counts, genome = BSgenome.Hsapiens.UCSC.hg38)}
  if(organism == "GRCm38") {
     chromVAR_counts = addGCBias(chromVAR_counts, genome = BSgenome.Mmusculus.UCSC.mm10)
     motif_ix = matchMotifs(mouse_pwms_v2, chromVAR_counts, genome = BSgenome.Mmusculus.UCSC.mm10)}  
  dev = computeDeviations(object = chromVAR_counts, annotations = motif_ix)
  variability = computeVariability(dev)
  zscore = assays(dev)$z
  rownames(zscore) = variability[,1]
  return(list(dev = dev, zscore = zscore))
}

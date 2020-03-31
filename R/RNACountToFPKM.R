#' Convert scRNA-seq count matrix to FPKM matrix
#'
#' Convert scRNA-seq count matrix to FPKM matrix for downstream analysis.
#'
#' @docType methods
#' @name RNACountToFPKM
#' @rdname RNACountToFPKM
#'
#' @param countMat Input count matrix, with genes as rows and cells as columns.
#' @param idType Type of gene name, available options are "Ensembl" or "Symbol". Default is "Symbol".
#' @param organism Ogranism used, currently only support human "GRCh38" and mouse "GRCm38". Default is "GRCh38".
#'
#' @author Chenfei Wang, Dongqing Sun
#'
#' @return A matrix of FPKM, with genes as rows and cells as columns.
#'
#'
#' @examples
#' data(pbmc.RNA)
#' pbmc.RNA.FPKM <- RNACountToFPKM(pbmc.RNA, idType = "Symbol", organism = "GRCh38")
#' head(pbmc.RNA.FPKM[,1:10])
#'
#' @export

RNACountToFPKM <- function(countMat, idType = "Ensembl", organism = "GRCh38")
{
   if(organism=="GRCh38"){
     data(GRCh38.ensembl)
     ensembl <- GRCh38.ensembl}
   else{
     data(GRCm38.ensembl)
     ensembl <- GRCm38.ensembl}
   ensembl$Length <- abs(ensembl$`Gene.end..bp.` - ensembl$`Gene.start..bp.`)
   if(toupper(idType) == "ENSEMBL"){
      len <- ensembl[match(rownames(countMat),ensembl$`Gene.stable.ID`), "Length"]
      count_rowname = ensembl[match(rownames(countMat),ensembl$`Gene.stable.ID`), "Gene.name"]
      count_index = which(!duplicated(count_rowname) & !is.na(count_rowname))
      countMat = countMat[count_index,]
      rownames(countMat) = count_rowname[count_index]
      len = len[count_index]
   } else if(toupper(idType) == "SYMBOL")
     len <- ensembl[match(rownames(countMat),ensembl$`Gene.name`), "Length"]
   else
     stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")
   
   na_idx = which(is.na(len))
   if(length(na_idx)>0){
     warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
     countMat = countMat[!is.na(len),]
     len = len[!is.na(len)]
   }
   tmp <- countMat / len
   FPKM <- 1e6 * Matrix::t(Matrix::t(tmp) / Matrix::colSums(countMat))
   FPKM = FPKM[!duplicated(rownames(FPKM)),]
   return(FPKM)
}
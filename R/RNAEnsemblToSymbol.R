#' Convert Ensembl ID to Gene Symbol
#'
#' Convert Ensembl ID to Gene Symbol for count matrix.
#'
#' @docType methods
#' @name RNAEnsemblToSymbol
#' @rdname RNAEnsemblToSymbol
#'
#' @param countMat Input count matrix, with genes as rows and cells as columns.
#' @param organism Ogranism used, currently only support human "GRCh38" and mouse "GRCm38". Default is "GRCh38".
#'
#' @author Chenfei Wang, Dongqing Sun
#'
#' @return A matrix of count, with genes as rows and cells as columns.
#'
#'
#' @examples
#' data(pbmc.RNA)
#' pbmc.RNA.symbol <- RNAEnsemblToSymbol(pbmc.RNA, organism = "GRCh38")
#' head(pbmc.RNA.symbol)
#'
#' @export

RNAEnsemblToSymbol <- function(countMat, organism = "GRCh38")
{
  if(organism=="GRCh38"){
    data(GRCh38.ensembl)
    ensembl <- GRCh38.ensembl}
  else{
    data(GRCm38.ensembl)
    ensembl <- GRCm38.ensembl}
  count_rowname = ensembl[match(rownames(countMat),ensembl$Gene.stable.ID), "Gene.name"]
  count_index = which(!duplicated(count_rowname) & !is.na(count_rowname))
  countMat = countMat[count_index,]
  rownames(countMat) = count_rowname[count_index]
  na_idx = which(is.na(count_rowname))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which symbol is not available !"))
  }
  return(countMat)
}
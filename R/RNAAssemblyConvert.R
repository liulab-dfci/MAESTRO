#' Assembly Conversion
#'
#' Convert gene symbol between genome assemblies for count matrix.
#'
#' @docType methods
#' @name RNAAssemblyConvert
#' @rdname RNAAssemblyConvert
#'
#' @param countMat Input count matrix, with genes as rows and cells as columns.
#' @param from The genome assembly of input genes. For human, "GRCh37" (hg19) and "GRCh38" (hg38) are supported.
#' For mouse, "NCBIM37" (mm9) and "GRCm38" (mm10) are supported. Default is "GRCh37".
#' @param to The genome assembly of output genes. For human, "GRCh37" (hg19) and "GRCh38" (hg38) are supported.
#' For mouse, "NCBIM37" (mm9) and "GRCm38" (mm10)" are supported. Default is "GRCh38".
#' @param organism Ogranism used, currently only support "Human" and "Mouse". Default is "Human".
#'
#' @author Dongqing Sun
#'
#' @return A matrix of count with gene names converted.
#'
#' @export
#' 

RNAAssemblyConvert <- function(countMat, from = "GRCh37", to = "GRCh38", organism = "Human"){
  countMat = as(as.matrix(countMat), "dgCMatrix")
  genes_original = rownames(countMat)
  if (organism == "Human"){
    data(GRCh37.GRCh38)
    ensembl = GRCh37.GRCh38
    if (from == "GRCh37" & to == "GRCh38") {
      genes_convert = ensembl[match(genes_original, ensembl$Gene.name.GRCh37),"Gene.name.GRCh38"]
    }
    if (from == "GRCh38" & to == "GRCh37") {
      genes_convert = ensembl[match(genes_original, ensembl$Gene.name.GRCh38),"Gene.name.GRCh37"]
    }
  }
  
  if (organism == "Mouse"){
    data(NCBIM37.GRCm38)
    ensembl = NCBIM37.GRCm38
    if (from == "NCBIM37" & to == "GRCm38") {
      genes_convert = ensembl[match(genes_original, ensembl$Gene.name.NCBIM37),"Gene.name.GRCm38"]
    }
    if (from == "GRCm38" & to == "NCBIM37") {
      genes_convert = ensembl[match(genes_original, ensembl$Gene.name.GRCm38),"Gene.name.NCBIM37"]
    }
  }
  count_index = which((!duplicated(genes_convert)) & (!is.na(genes_convert)))
  countMat = countMat[count_index,]
  rownames(countMat) = genes_convert[count_index]
  return(countMat)
}
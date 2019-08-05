#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' Convert gene symbol to entreID, imported form MAGeckFlute.
#'
#' @docType methods
#' @name RNAEnsemblToEntrez
#' @rdname RNAEnsemblToEntrez
#' @aliases RNAEnsemblToEntrez
#'
#' @param genes A character vector, input genes to be converted.
#' @param fromType The input ID type, one of "Symbol" (default), "Entrez" and "Ensembl";
#' you can also input other valid attribute names for biomaRt.
#' @param toType The output ID type, one of "Symbol", "Entrez" (default), "Ensembl";
#' you can also input other valid attribute names for biomaRt.
#' @param organism One of "hsa"(or 'Human'), "mmu"(or 'Mouse'), "bta", "cfa", "ptr", "rno", and "ssc".
#' @param useBiomart Boolean, indicating whether use Biomart to do the transformation.
#' @param ensemblHost String, specifying ensembl host, you can use `listEnsemblArchives()`
#' to show all available Ensembl archives hosts.
#'
#' @return A character vector, named by unique input gene ids.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link[pathview]{eg2id}}
#'
#' @examples
#' data(mle.gene_summary)
#' RNAEnsemblToEntrez(mle.gene_summary$Gene[1:10], organism="hsa")
#'
#' @import biomaRt
#' @export

RNAEnsemblToEntrez <- function(genes, fromType="Symbol", toType="Entrez",
                        organism="hsa", useBiomart = FALSE,
                        ensemblHost = "www.ensembl.org"){
  requireNamespace("biomaRt")
  stopifnot(length(genes)>0)
  bods <- data.frame(package = paste0("org.", c("Hs", "Mm", "Rn", "Bt", "Cf", "Pt", "Ss"), ".eg.db"),
                     species = c("Human", "Mouse", "Rat", "Bovine", "Canine", "Chimp", "Pig"),
                     "kegg code" = c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc"),
                     check.names = FALSE, stringsAsFactors = FALSE)
  ridx = c(which(tolower(bods[,2])==tolower(organism)),
           which(tolower(bods[,3])==tolower(organism)))
  stopifnot(length(ridx)==1)
  org = bods[ridx,3]
  genes = as.character(genes)
  fromType = tolower(fromType)
  toType = tolower(toType)
  #=============Read annotation file========
  if(useBiomart){
    datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                        "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
    type = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol")
    if(org=="mmu") type[3] = "mgi_symbol"
    names(type) = c("ensembl", "entrez", "symbol")
    fromType = ifelse(fromType%in%names(type), type[fromType], fromType)
    toType = ifelse(toType%in% names(type), type[toType], toType)
    # listEnsemblArchives()
    # listMarts()
    ds = datasets[grepl(org, datasets)]
    ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds, host = ensemblHost)
    GeneID_Convert = getBM(attributes=c(fromType, toType), mart = ensembl)
    GeneID_Convert[, toType] = GeneID_Convert[, toType]
  }else{
    GeneID_Convert <- RNAEnsemblToEntrezGetOrig(org)$Symbol_Entrez
  }
  ##======Convert==========
  idx = duplicated(GeneID_Convert[, fromType])
  convert = GeneID_Convert[!idx, toType]
  names(convert) = GeneID_Convert[!idx, fromType]
  gene_after = as.character(convert[genes])
  names(gene_after) = genes
  return(gene_after)
}

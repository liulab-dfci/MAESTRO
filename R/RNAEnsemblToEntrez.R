#' Gene ID conversion between ENTREZID and SYMBOL
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
    GeneID_Convert <- getOrg(org)$Symbol_Entrez
  }
  ##======Convert==========
  idx = duplicated(GeneID_Convert[, fromType])
  convert = GeneID_Convert[!idx, toType]
  names(convert) = GeneID_Convert[!idx, fromType]
  gene_after = as.character(convert[genes])
  names(gene_after) = genes
  return(gene_after)
}

getOrg <- function(organism, update = FALSE){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  require(data.table)
  bods <- data.frame(package = paste0("org.", c("Hs", "Mm", "Rn", "Bt", "Cf", "Pt", "Ss"), ".eg.db"),
                     species = c("Human", "Mouse", "Rat", "Bovine", "Canine", "Chimp", "Pig"),
                     "kegg code" = c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc"),
                     check.names = FALSE, stringsAsFactors = FALSE)
  res=list()
  ##======
  # Get the mapping from organism to package
  ridx = c(which(tolower(bods[,2])==tolower(organism)), which(tolower(bods[,3])==tolower(organism)))
  stopifnot(length(ridx)==1)
  res$org = bods[ridx,3]
  res$pkg = bods[ridx,1]

  gzfile = c("Homo_sapiens.gene_info.gz", "Bos_taurus.gene_info.gz",
             "Canis_familiaris.gene_info.gz", "Mus_musculus.gene_info.gz",
             "Pan_troglodytes.gene_info.gz", "Rattus_norvegicus.gene_info.gz",
             "Sus_scrofa.gene_info.gz")
  names(gzfile) = c("hsa", "bta", "cfa", "mmu", "ptr", "rno", "ssc")
  locfname <- file.path(system.file("extdata", package = "MAGeCKFlute"), gzfile[res$org])

  if((!file.exists(locfname)) | update){
    ## Download gene information from NCBI ftp server
    remfname <- paste0("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/", gzfile[res$org])
    download.file(remfname, locfname, quiet = TRUE)
  }
  ## Reorder the mapping file
  data = read.csv(gzfile(locfname), sep = "\t", header = TRUE,
                  quote = "", stringsAsFactors = FALSE, comment.char = "")
  data = data[, c("GeneID", "Symbol", "Synonyms", "dbXrefs", "type_of_gene", "description")]

  data2 = apply(data, 1, function(x){
    Ensembl = "-"; HGNC = "-"
    if(grepl("HGNC", x[4])){
      HGNC = gsub(".*HGNC:", "", x[4])
      HGNC = gsub("\\|.*", "", HGNC)
    }
    if(grepl("Ensembl", x[4])){
      Ensembl = gsub(".*Ensembl:", "", x[4])
      Ensembl = gsub("\\|.*", "", Ensembl)
    }
    tmp = unlist(strsplit(x[3], "[|]"))
    as.vector(rbind(as.integer(x[1]), x[2], tmp, HGNC, Ensembl, x[5], x[6]))
  })
  data2 = matrix(unlist(data2), ncol=7, byrow = TRUE)
  colnames(data2) = c("entrez", "symbol", "synonyms", "hgnc",
                      "ensembl", "type", "fullname")
  tmp = melt(as.data.frame(data2, stringsAsFactors = FALSE),
             measure.vars=c("symbol", "synonyms"))
  mapping = tmp[, c("entrez", "value", "ensembl", "hgnc", "type", "fullname")]
  colnames(mapping)[2] = "symbol"
  idx = duplicated(paste(mapping$entrez, mapping$symbol, mapping$ensembl,
                         mapping$hgnc, mapping$type, mapping$fullname, sep = "_"))
  mapping = mapping[!idx, ]
  idx = mapping$symbol=="-" & mapping$ensembl=="-" & mapping$hgnc=="-" &
    mapping$type=="-" & mapping$fullname=="-"
  mapping = mapping[!idx, ]
  mapping$symbol = mapping$symbol
  res$Symbol_Entrez = mapping

  return(res)
}

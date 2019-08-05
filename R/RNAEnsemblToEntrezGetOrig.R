#' Determine the gene annotation package.
#'
#' Determine the gene annotation package for specific organism, imported form MAGeckFlute.
#'
#' @docType methods
#' @name RNAEnsemblToEntrezGetOrig
#' @rdname RNAEnsemblToEntrezGetOrig
#'
#' @param organism Character, KEGG species code, or the common species name, used to determine
#' the gene annotation package. For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#' @param update Boolean, indicating whether download recent annotation from NCBI.
#' @return A list containing three elements:
#' \item{organism}{species}
#' \code{pkg}{annotation package name}
#' \code{Symbol_Entrez}{a data frame, mapping between gene symbol and entrez id}
#'
#' @author Wubing Zhang
#'
#' @examples
#' ann = RNAEnsemblToEntrezGetOrig("human")
#' print(ann$pkg)
#'
#' @importFrom data.table melt
#' @export

RNAEnsemblToEntrezGetOrig <- function(organism, update = FALSE){
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

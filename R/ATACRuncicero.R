#' Generate gene score matrix for scATAC-seq dataset using cicero
#'
#' Generate gene score matrix for scATAC-seq dataset using cicero.
#'
#' @docType methods
#' @name ATACRuncicero
#' @rdname ATACRuncicero
#'
#' @param inputMat Input unnormalized count matrix, with three columns, peak, cell and count. Peaks should be like "chromosome_peakstart_peakend",
#' for example "chr10_100020591_100020841".
#' @param project Output project name for the cicero object. Default is "MAESTRO.scATAC.cicero".
#'
#' @author Chenfei Wang
#'
#' @return A list contains a cicero object for ATAC clusters and a list with clustering informations.
#'
#'
#' @examples
#' data(pbmc.ATAC)
#' pbmc.ATAC.res <- ATACRuncicero(inputMat = pbmc.ATAC, project = "PBMC.scATAC.cicero")
#' str(pbmc.ATAC.res$ATAC)
#' head(pbmc.ATAC.res$cluster)
#'
#' @export

ATACRuncicero <- function(peakCountFile, outPrefix)
{
  # cicero annotate_cds_by_site build_gene_activity_matrix detectGenes estimateSizeFactors make_atac_cds make_cicero_cds normalize_gene_activities pData reducedDimA reduceDimension run_cicero
  library(cicero)
  peakCount = read.delim(peakCountFile, header=F)
  colnames(peakCount) = c("Peak","Cell","Count")
  input_cds = make_atac_cds(peakCount, binarize = TRUE)

  set.seed(2017)
  input_cds = detectGenes(input_cds)
  input_cds = estimateSizeFactors(input_cds)
  input_cds = reduceDimension(input_cds, max_components = 2, num_dim=10, reduction_method = 'tSNE', norm_method = "none")
  tsne_coords = t(reducedDimA(input_cds))
  row.names(tsne_coords) = row.names(pData(input_cds))
  cicero_cds = make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

  hg38.genome = read.table('/homes/cwang/annotations/hg38/hg38.len')
  conns = run_cicero(cicero_cds, hg38.genome,  window = 5e+05, sample_num = 100) # Takes a few minutes to run
  gene = read.table('/homes/cwang/annotations/hg38/hg38.refseq.bed')
  gene_annotation_sub = unique(gene[,c(1:3,5)])
  names(gene_annotation_sub)[4] = "gene"
  input_cds = annotate_cds_by_site(input_cds, gene_annotation_sub)
  unnorm_ga = build_gene_activity_matrix(input_cds, conns)
  num_genes = pData(input_cds)$num_genes_expressed
  names(num_genes) = row.names(pData(input_cds))
  cicero_gene_activities = normalize_gene_activities(unnorm_ga, num_genes)
  write.table(as.matrix(cicero_gene_activities),paste0(outPrefix,".txt"),sep="\t",quote=F)
}
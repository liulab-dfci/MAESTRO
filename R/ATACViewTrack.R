#' Plot scATACseq coverage for each cluster from tabix indexed fragment.tsv.gz file
#'
#' Reads in the fragments.tsv.gz or bam file for only the given region (gene), and plot using \href{http://bioconductor.org/packages/release/bioc/html/karyoploteR.html}{karyoploteR} package
#' 
#' @docType methods
#' @name ATACViewTrack
#' @rdname ATACViewTrack
#' 
#' @param chrom chromosome
#' @param start chromosome start
#' @param end chromosome end
#' @param gene_name gene symbol for the gene one want's to plot. e.g. VEGFA. Use either chrom, start, end or
#' gene_name
#' @param upstream If use gene_name, specify the upstream bp for extending the GRanges
#' @param downstream If use gene_name, specify the downstream bp for extending the GRanges
#' @param fragment path to the fragment.tsv.gz file from 10x cellranger-atac output. indexed by tabix.
#' @param bam path to the bam file
#' @param cellbarcode_tag The cell barcode tag in the bam file, default is "CB" as in the 10x output
#' @param mapqFilter the mapq quality score filter for reads from bam file, default is 30,
#' @param grouping path to a tsv file or a dataframe object with three columns: "cell", "cluster", "depth". "cell" is the
#' cell barcode, "cluster" is the cluster id for each cell, and "depth" is the number of reads in that cell.
#' @param clusters_to_plot A character vector containing clusters to plot in the tracks from bottom to top. default is all clusters
#' @param genome hg19, hg38 or mm9, mm10
#' @param txdb A TxDb object, e.g. TxDb.Hsapiens.UCSC.hg38.knownGene
#' @param eg.db Either org.Hs.eg.db or org.Mm.eg.db
#' @param ymax ymax for each track, if not specified, max of all the tracks will be calculated. Every
#' track will use the same ymax, so it is comparable across clusters.
#' @param label_cex size of the cluster label
#' @param label_side side of the cluster label. "left" or "right".
#' @param label.margin the position of the cluster label, sepcify negative e.g., -0.6
#' to move to the center of the track body if label_side is left.
#' @param yaxis_cex size of the y-axis
#' @param track_cols a vector of colors for the tracks, should be the same length as clusters_to_plot or
#' the total number of clusters. If clusters_to_plot is NULL, and track_cols has only one element, all tracks will be plotted with
#' the same color.
#' @param tick.dist chromosome tick distance to mark, default 10kb
#' @param minor.tick.dist minor chromosome tick distance to mark, default 2000 bp
#' @param tick_label_cex size of the tick label
#'
#' @author Ming Tang, Dongqing Sun
#' 
#' @return A karyoploteR object
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @export
#'
#' @examples
#' \dontrun{
#' library(org.Hs.eg.db)
#' library(org.Mm.eg.db)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' chrom<-  "chr12"
#' start<-  69730394
#' end<- 69760971
#' ATACViewTrack(chrom = chrom, start = start, end = end, fragment = "data/atac_viz/10k_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz",
#'                     grouping = group_data, track_cols = "red")

#' ATACViewTrack(gene_name = "NKG7", fragment = "data/atac_viz/10k_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz",
#'                     grouping = "data/atac_viz/grouping.txt", tick_label_cex = 1, tick.dist = 5000,
#'                     minor.tick.dist = 1000)
#'
#' ATACViewTrack(gene_name = "NKG7", fragment = "data/atac_viz/10k_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz",
#'                     grouping = "data/atac_viz/grouping.txt", tick_label_cex = 1, tick.dist = 5000,
#'                     minor.tick.dist = 1000, track_cols = c("red", "blue", ""cadetblue2"),
#'                     clusters_to_plot = c("1", "3, "4"))
#'
#' ATACViewTrack(gene_name = "NKG7", bam = "data/atac_viz/10k_pbmc/atac_v1_pbmc_10k.bam",
#'                     cellbarcode_tag= "CB", mapqFilter = 30,
#'                     grouping = "data/atac_viz/grouping.txt", tick_label_cex = 1, tick.dist = 5000,
#'                     minor.tick.dist = 1000)
#'}

ATACViewTrack <- function(chrom = NULL, start = NULL, end =NULL, gene_name = NULL,
                               upstream = 2000, downstream = 2000, fragment = NULL,
                               bam = NULL,
                               cellbarcode_tag = "CB",
                               mapqFilter = 30,
                               grouping,
                               clusters_to_plot = NULL,
                               genome ='hg38',
                               txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                               eg.db = org.Hs.eg.db,
                               ymax = NULL, label_cex = 1, label_side = "left",
                               label.margin = 0.05,
                               yaxis_cex = 1, track_cols = "cadetblue2",
                               tick.dist = 10000, minor.tick.dist = 2000,
                               tick_label_cex = 1
){
  
  ######## read in the clustering information ##################
  if (class(grouping) == "data.frame") {
    if(! all(c("cell", "cluster", "depth") %in% colnames(grouping))) {
      stop('Grouping dataframe must have cell, cluster, and depth columns.')
    }
    grouping$cell = as.character(grouping$cell)
    grouping$cluster = as.character(grouping$cluster)
    grouping = as_tibble(grouping)
  } else {
    grouping<- readr::read_tsv(grouping)
    if(! all(c("cell", "cluster", "depth") %in% colnames(grouping))) {
      stop('Grouping dataframe must have cell, cluster, and depth columns.')
    }
  }

  ## get number of reads per group for normalization.
  ## not furthur normalize by the cell number in each group.
  message("reading in the grouping/clustering information and calculate
  a scaling factor for normalizing tracks based on 1e6/read depth")
  
  grouping<-  grouping %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(cells_in_group = dplyr::n(), total_depth_in_group = sum(depth)) %>%
    # reads per million (RPM)
    dplyr::mutate(scaling_factor = 1e6/(total_depth_in_group)) %>%
    dplyr::ungroup() %>%
    dplyr::select(cell, cluster, scaling_factor)
  
  if (!is.null(clusters_to_plot)){
    if (any(!(clusters_to_plot %in% grouping$cluster))){
      stop("clusters_to_plot should be in the grouping txt file")
    }
    if (length(track_cols) ==1){
      track_cols<- rep(track_cols, length(clusters_to_plot))
    } else if (length(clusters_to_plot) != length(track_cols)){
      stop("the track_cols should be the same length as the clusters_to_plot")
    }
  } else {
    total_number_clusters<- length(unique(grouping$cluster))
    # if only one color is specified, all tracks will be plotted with the same color
    if (length(track_cols) == 1){
      track_cols<- rep(track_cols, total_number_clusters)
    } else if (length(track_cols) != total_number_clusters) {
      stop("track_cols should be of length 1 or the same length as the total number of clusters")
    }
  }
  
  
  ########## get a GenomicRanges for the gene and extend upstream and downstream ########
  if (is.null(chrom) & is.null(start) & is.null(end) & !is.null(gene_name)){
    gene <- addGeneNameToTxdb(txdb = txdb, eg.db = eg.db)
    gr<- gene[which(gene$symbol == gene_name)]
    if (length(gr) == 0){
      stop("gene name is not found in the database")
    } else if (length(gr) > 1) {
      gr<- gr[1]
      warning("multiple GRanges found for the gene, using the first one")
      gr<- extend(gr, upstream = upstream, downstream = downstream)
    } else {
      gr<- extend(gr, upstream = upstream, downstream = downstream)
    }
    
  } else if (!is.null(chrom) & !is.null(start) & !is.null(end)){
    gr<- GRanges(seq = chrom, IRanges(start = start, end = end ))
  }
  
  gr<- GenomeInfoDb::keepStandardChromosomes(gr)
  
  #### read in the fragment.tsv.gz with "chr", "start", "end", "cell", "duplicate" columns. output from cellranger-atac ###
  #### or bam file with a CB tag for cell barcode ###############
  #### if it is a bam file, need to read in as pairs of read ######
  
  if (is.null(fragment) & !is.null(bam)) {
    ## read in only properly paired reads
    message("reading in bam file and extracting reads")
    flag<- Rsamtools::scanBamFlag(isProperPair =TRUE)
    param <- Rsamtools::ScanBamParam(what=Rsamtools::scanBamWhat(), which= gr, tag= cellbarcode_tag, flag = flag, mapqFilter = 30)
    
    # read in pairs as a list to retain meta columns
    GA <- GenomicAlignments::readGAlignmentsList(file = bam, param=param, use.names =TRUE)
    # this takes a bit longer than I want: ~1min
    meta_cols<- Reduce("rbind", lapply(GA, function(x) mcols(x)[1,]))
    meta_cols<- dplyr::left_join(as.data.frame(meta_cols), grouping, by =setNames("cell", cellbarcode_tag))
    
    ## same as: as(y, "GRanges)
    grs<- GenomicRanges::granges(GA, ignore.strand =TRUE, use.mcols = TRUE)
    grs$cell<- meta_cols[, cellbarcode_tag]
    grs$cluster<- meta_cols[, "cluster"]
    grs$scaling_factor<- meta_cols[, "scaling_factor"]
    
    grs<- grs[!is.na(grs$cell)]
    
    ## use {plyranges} in the future
    message("shifting reads +4bp and -5bp and removing PCR duplicates")
    reads<- grs %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(cluster)) %>%
      dplyr::group_by(cell) %>%
      ## remove PCR duplicates.
      dplyr::distinct(seqnames, start, end, .keep_all = TRUE) %>%
      dplyr::ungroup() %>%
      ### shift + strand 4bp and - strand -5 bp for the 9bp overhang from ATACseq
      dplyr::mutate(start = start + 4, end = end -5) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
  } else if (is.null(bam) & !is.null(fragment)){
    message("reading in fragment.tsv and extracting reads, no need to shift reads")
    reads<- Rsamtools::scanTabix(fragment, param = gr)
    reads<- reads[[1]] %>%
      tibble::enframe() %>%
      dplyr::select(-name) %>%
      tidyr::separate(value, into = c("chr", "start", "end", "cell", "duplicate"), sep = "\t") %>%
      dplyr::mutate_at(.vars = c("start", "end"), as.numeric) %>%
      # make it 1 based for R, the fragment.tsv is 0 based bed file
      dplyr::mutate(start = start + 1) %>%
      dplyr::inner_join(grouping) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  }
  
  
  #####  plotting by krayoplotR ##############
  
  message("plotting the tracks")
  pp <- karyoploteR::getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  kp <- karyoploteR::plotKaryotype(genome = genome, zoom = gr, plot.params = pp)
  kp<- karyoploteR::kpAddBaseNumbers(kp, tick.dist = tick.dist, minor.tick.dist = minor.tick.dist,
                                     add.units = TRUE, cex= tick_label_cex, digits = 6)
  ## calculate the normalized coverage
  normalized_coverage<- function(x){
    #if (!is(x, "GRangesList"))
    #  stop("'x' must be a GRangesList object")
    # specify the width to the whole chromosome to incldue the 0s
    cvgs<- lapply(x, function(x) GenomicRanges::coverage(x, width = kp$chromosome.lengths) * x$scaling_factor[1])
    return(cvgs)
  }
  
  # drop seqlevels for not used chromosomes, otherwise coverage(x, width = kp$chromosome.lengths) not working
  reads<- GenomeInfoDb::keepSeqlevels(reads, value = kp$chromosomes)
  # GRangesList object by group/cluster
  reads_by_group<- split(reads, reads$cluster)
  #print(reads_by_group)
  coverage_norm<- normalized_coverage(reads_by_group)
  
  ## calculate the max coverage if not specified
  if (is.null(ymax)) {
    yaxis_common<- ceiling(max(lapply(coverage_norm, max) %>% unlist()))
  } else {
    yaxis_common<- ymax
  }
  ## add gene information
  genes.data <- karyoploteR::makeGenesDataFromTxDb(txdb,
                                                   karyoplot=kp,
                                                   plot.transcripts = TRUE,
                                                   plot.transcripts.structure = TRUE)
  genes.data <- karyoploteR::addGeneNames(genes.data)
  genes.data <- karyoploteR::mergeTranscripts(genes.data)
  
  kp<- karyoploteR::kpPlotGenes(kp, data=genes.data, r0=0, r1=0.05, gene.name.cex = 1)
  
  ## subsetting the clusters to plot and order as the user input
  if (!is.null(clusters_to_plot)){
    coverage_norm<- coverage_norm[clusters_to_plot]
  }
  
  for(i in seq_len(length(coverage_norm))) {
    read <- coverage_norm[[i]]
    at <- karyoploteR::autotrack(i, length(coverage_norm), r0=0.1, r1=1, margin = 0.1)
    karyoploteR::kpPlotCoverage(kp, data=read,
                                r0=at$r0, r1=at$r1, col = track_cols[i], ymax = yaxis_common)
    karyoploteR::kpAxis(kp, ymin=0, ymax=yaxis_common, numticks = 2, r0=at$r0, r1=at$r1, cex = yaxis_cex, labels = c("", yaxis_common))
    karyoploteR::kpAddLabels(kp, side = label_side, labels = names(coverage_norm)[i], r0=at$r0, r1=at$r1,
                             cex=label_cex, label.margin = label.margin)
  }
}


#' Extend a GRanges object upstream and downstream
#'
#' @param x A GRanges object
#' @param upstream bp for extending upstream
#' @param downstream bp for extending downstream
#'
#' @return An extended GRanges object
#'
extend <- function(x, upstream=0, downstream=0)
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}


#' Add Gene sybmol to the GRanges object after calling genes(txdb)
#'
#' @param txdb A TxDb object, e.g. TxDb.Hsapiens.UCSC.hg19.knownGene
#' @param eg.db Either org.Hs.eg.db or org.Mm.eg.db
#'
#' @return An GRanges object for all the genes with gene symbol in the metadata column
#'
addGeneNameToTxdb<- function(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                             eg.db = org.Hs.eg.db){
  gene<- GenomicFeatures::genes(txdb)
  ## 1: 1 mapping
  ss<- AnnotationDbi::select(eg.db, keys = gene$gene_id,
                             keytype="ENTREZID", columns = "SYMBOL" )
  gene$symbol<- ss[, 2]
  return(gene)
}


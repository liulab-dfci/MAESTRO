library(MAESTRO)
library(Seurat)
library(Matrix)
library(future)
library(future.apply)
library(pbapply)
library(optparse)
library(ggplot2)


option_list = list(
  make_option(c("--prefix"), type = "character", default = "MAESTRO",
              action = "store", help = "The prefix of the output files."
  ),
  make_option(c("--peakcount"), type = "character", default = "",
              action = "store", help = "The peak count matrix h5 file."
  ),
  make_option(c("--rpmatrix"), type = "character", default = "",
              action = "store", help = "The regulatory potential matrix h5 file."
  ),
  make_option(c("--outdir"), type = "character", default = "Result/Analysis",
              action = "store", help = "The directory where the output files are stored."
  ),
  make_option(c("--species"), type = "character", default = "GRCh38",
              action = "store", help = "The platform of scRNA-seq."
  ),
  make_option(c("--annotation"), type = "character", default = "False",
              action = "store", help = "Whether or not to annotate cell types. Default is False."
  ),
  make_option(c("--method"), type = "character", default = "RP-based",
              action = "store", help = "Method to annotate cell types. ['RP-based', 'peak-based', 'both']"
  ),
  make_option(c("--signature"), type = "character", default = "",
              action = "store", help = "The cell signature file for celltype annotation. Default is built-in CIBERSORT immune cell signature."
  ),
  make_option(c("--gigglelib"), type = "character", default = "",
              action = "store", help = "Annotation to run rabit (only if method is set to rabit)."
  ),
  make_option(c("--fragment"), type = "character", default = "fragments_corrected_count.tsv.gz",
              action = "store", help = "Gzip-compressed fragment file."
  ),
  make_option(c("--thread"), type = "integer", default = 1,
              action = "store", help = "Number of cores to use."
  )
)

argue = parse_args(OptionParser(option_list = option_list, usage = "Run scATAC-seq analysis pipeline"))

setwd(argue$outdir)
count_mat = argue$peakcount
rp_mat = argue$rpmatrix
prefix = argue$prefix
thread = argue$thread
sigfile = argue$signature
gigglelib = argue$gigglelib
species = argue$species
fragment = argue$fragment
annotation = argue$annotation
method = argue$method

# countmatrix = read.table(argue[1], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
# RPmatrix = read.table(argue[2], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
countmatrix = Read10X_h5(count_mat)
RPmatrix = Read10X_h5(rp_mat)

plan("multiprocess", workers = as.integer(thread))
options(future.globals.maxSize = 10*1024^3)

if(sigfile %in% c("human.immune.CIBERSORT", "mouce.brain.ALLEN", "mouse.all.facs.TabulaMuris", "mouse.all.droplet.TabulaMuris")){
  signatures = sigfile
}else{
  signatures = read.table(sigfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
}
result = ATACRunSeurat(inputMat = countmatrix, project = prefix, method = "LSI")
result$ATAC = ATACAttachGenescore(ATAC = result$ATAC, RPmatrix = RPmatrix)
saveRDS(result, paste0(prefix, "_scATAC_Object.rds"))
if (annotation == "True") {
  if (method == "RP-based") {
    result$ATAC = ATACAnnotateCelltype(result$ATAC, signatures = signatures)

  }
  if (method == "peak-based") {
    result$ATAC <- ATACAnnotateChromatinAccessibility(ATAC = result$ATAC, 
                                                      peaks = result$peaks, 
                                                      project = prefix, 
                                                      giggle.path = gigglelib,
                                                      organism = species)
  }
  if (method == "both") {
    result$ATAC <- ATACAnnotateChromatinAccessibility(ATAC = result$ATAC, 
                                                      peaks = result$peaks, 
                                                      project = prefix, 
                                                      giggle.path = gigglelib,
                                                      organism = species)
    result$ATAC = ATACAnnotateCelltype(result$ATAC, signatures = signatures)

  }
  if ("assign.biological_resource" %in% colnames(result$ATAC@meta.data)) {
    p1 <- DimPlot(result$ATAC, label = TRUE, reduction = "umap", group.by = "assign.biological_resource", repel=T, pt.size = 0.5, label.size = 2.5)
    ggsave(file.path(paste0(result$ATAC@project.name, "_CistromeTop_annotated.png")), p1, width=7.5, height=4)
  }
}

result.tfs = ATACAnnotateTranscriptionFactor(ATAC = result$ATAC,
                                             peaks = result$peaks, 
                                             project = prefix, 
                                             giggle.path = gigglelib, 
                                             organism = species)

saveRDS(result, paste0(prefix, "_scATAC_Object.rds"))

if(species == "GRCh38"){
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  genome = "hg38"
}
if(species == "GRCm38"){
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  genome = "mm10"
}

meta_info = data.frame(cell = rownames(result[["ATAC"]]@meta.data),
                       cluster = result[["ATAC"]]@meta.data$seurat_clusters,
                       depth = result[["ATAC"]]@meta.data$nCount_ATAC)

GeneTrack = function(prefix, gene, fragment, group_file, txdb, genome){
  genetrack_file = paste(prefix, gene, "genetrack.png", sep = "_")
  tryCatch(expr = {
    png(genetrack_file, units = "in", width = 6, height = 5, res = 300)
    ATACViewTracks(gene_name = gene, downstream = 8000, 
                        yaxis_cex = 1,
                        fragment = fragment,
                        grouping = group_file,
                        tick_label_cex = 1, tick.dist = 5000,
                        track_cols = "blue", 
                        label_cex = 1,
                        minor.tick.dist = 1000, label.margin = -0.6,
                        txdb = txdb,
                        genome = genome)
    dev.off()
  }, error = function(e){
    print(e)
  })
}

GeneTrack(prefix, "MS4A1", fragment, meta_info, txdb, genome)
GeneTrack(prefix, "CD3D", fragment, meta_info, txdb, genome)


cluster_info = data.frame(Cell = rownames(result$ATAC@meta.data), Cluster = paste0("Cluster_", result$ATAC@meta.data$seurat_clusters), stringsAsFactors = FALSE)
write.table(cluster_info, paste0(prefix, "_cell_cluster.txt"), sep  = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


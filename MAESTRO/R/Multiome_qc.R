library(MAESTRO)
library(optparse)

option_list = list(
  make_option(c("--prefix"), type = "character", default = "MAESTRO",
              action = "store", help = "The prefix of the output files."
  ),
  make_option(c("--outdir"), type = "character", default = "MAESTRO",
              action = "store", help = "The directory where the output files are stored."
  ),
  make_option(c("--filtering"), type = "character", default = "",
              action = "store", help = "The result of scRNA_qc."
  ),
  make_option(c("--singlestat"), type = "character", default = "singlecell.txt",
              action = "store", help = "The result of MAESTRO QC."
  ),
  make_option(c("--count"), type = "character", default = "",
              action = "store", help = "The multiome count file."
  )
)
argue = parse_args(OptionParser(option_list = option_list, usage = "Generate QC plots."))

setwd(argue$outdir)
singlestat_file = argue$singlestat
countgene_file = argue$filtering
count_file = argue$count
prefix = argue$prefix

MultiomeFilteringPlot <- function(rna_statpath, atac_statpath, multiome_count, name){
  RCB_blue = "#2166AC"
  RCB_red = "#B2182B"
  
  UMI_gene = read.table(rna_statpath, header = TRUE, row.names = 1)
  mapping_matrix = read.table(atac_statpath,sep = "\t", row.names = 1, header = FALSE)
  count = Read10X_h5(multiome_count)
  cells = colnames(count$Peaks)
  all = intersect(rownames(UMI_gene), rownames(mapping_matrix))
  non_cells = all[!(all %in% cells)]

  png(paste0(name,"_multiome_cell_filtering.png"),width=5,height=5, res = 300, units = "in")
  par(mai = c(0.85, 0.85, 0.25, 0.25))
  plot(log10(UMI_gene[non_cells,1]+1),log10(mapping_matrix[non_cells,1]+1),
       xlim=c(0,log10(max(UMI_gene[,1])+1)+0.25), ylim=c(0,log10(max(mapping_matrix[,1])+1)+0.25), pch='.',col=RCB_red,ylab='ATAC transposition events in peaks per barcode (log10)',xlab='RNA UMIs per barcode (log10)')
  points(log10(UMI_gene[cells,1]+1),log10(mapping_matrix[cells,1]+1),
         pch='.',col=RCB_blue)
  legend("topleft",c("High-quality cells","Low-quality cells"),col=c(RCB_blue,RCB_red),pch=20, bty = "n")
  dev.off()
  # write.table(as.character(rownames(UMI_gene[which(UMI_gene[,1]>=UMI.cutoff&(UMI_gene[,2]>=gene.number.cutoff)),])), paste0(name,"_scRNA_validcells.txt"), sep = "\n", quote=F, row.names=F, col.names=F)
}
MultiomeFilteringPlot(countgene_file, singlestat_file, count_file, prefix)

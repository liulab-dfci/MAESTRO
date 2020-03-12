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
  make_option(c("--countcutoff"), type = "integer", default = 1000,
              action = "store", help = "Cutoff for the number of count in each cell."
  ),
  make_option(c("--genecutoff"), type = "integer", default = 500,
              action = "store", help = "Cutoff for the number of genes included in each cell."
  )
)
argue = parse_args(OptionParser(option_list = option_list, usage = "Generate QC plots."))

setwd(argue$outdir)
prefix = argue$prefix
countgene_file = argue$filtering
count_cutoff = argue$countcutoff
gene_cutoff = argue$genecutoff

RNAFilteringPlot(filepath = countgene_file, UMI.cutoff = count_cutoff, gene.number.cutoff = gene_cutoff, name = prefix)

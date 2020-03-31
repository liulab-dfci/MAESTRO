library(MAESTRO)
library(optparse)

option_list = list(
  make_option(c("--prefix"), type = "character", default = "MAESTRO",
              action = "store", help = "The prefix of the output files."
  ),
  make_option(c("--outdir"), type = "character", default = "MAESTRO",
              action = "store", help = "The directory where the output files are stored."
  ),
  make_option(c("--singlestat"), type = "character", default = "singlecell.txt",
              action = "store", help = "The result of MAESTRO QC."
  ),
  make_option(c("--countcutoff"), type = "integer", default = 1000,
              action = "store", help = "Cutoff for the number of count in each cell."
  ),
  make_option(c("--fripcutoff"), type = "double", default = NULL,
              action = "store", help = "Cutoff for fraction of reads in promoter."
  )
)
argue = parse_args(OptionParser(option_list = option_list, usage = "Generate QC filtering plots."))

singlestat_file = argue$singlestat
count_cutoff = argue$countcutoff
frip_cutoff = argue$fripcutoff
prefix = argue$prefix
setwd(argue$outdir)

ATACFilteringPlot(filepath = singlestat_file, name = prefix, 
                  reads.cutoff = count_cutoff, frip.cutoff = frip_cutoff)

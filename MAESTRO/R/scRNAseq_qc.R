library(MAESTRO)
library(optparse)

option_list = list(
  make_option(c("--prefix"), type = "character", default = "MAESTRO",
              action = "store", help = "The prefix of the output files."
  ),
  make_option(c("--rseqc"), type = "logical", default = FALSE,
              action = "store", help = "Do RSeQC or not."
  ),
  make_option(c("--outdir"), type = "character", default = FALSE,
              action = "store", help = "The directory where the output files are stored."
  ),
  make_option(c("--platform"), type = "character", default = FALSE,
              action = "store", help = "The platform of scRNA-seq."
  ),
  make_option(c("--bamstat"), type = "character", default = FALSE,
              action = "store", help = "The result of bam_stat."
  ),
  make_option(c("--readdistr"), type = "character", default = FALSE,
              action = "store", help = "The result of read_distribution."
  ),
  make_option(c("--qual"), type = "character", default = FALSE,
              action = "store", help = "The result of read_quality."
  ),
  make_option(c("--nvc"), type = "character", default = FALSE,
              action = "store", help = "The result of read_NVC."
  ),
  make_option(c("--gc"), type = "character", default = FALSE,
              action = "store", help = "The result of read_GC."
  ),
  make_option(c("--genecov"), type = "character", default = FALSE,
              action = "store", help = "The result of geneBody_coverage."
  ),
  make_option(c("--filtering"), type = "character", default = FALSE,
              action = "store", help = "The result of scRNA_qc."
  )
)
argue = parse_args(OptionParser(option_list = option_list, usage = "Generate QC plots."))

setwd(argue$outdir)
bamstat_file = argue$bamstat
readdistr_file = argue$readdistr
qualcode_file = argue$qual
nvc_file = argue$nvc
gc_file = argue$gc
genecov_file = argue$genecov
countgene_file = argue$filtering
prefix = argue$prefix
rseqc = argue$rseqc
platform = argue$platform

if(rseqc){
  RNAReadDistrPlot(bamstat.filepath = bamstat_file, readdistr.filepath = readdistr_file, name = prefix, platform = platform)
  RNAReadQualityPlot(filepath = qualcode_file, name = prefix)
  RNANucleotidePlot(filepath = nvc_file, name = prefix)
  RNAGCcontentPlot(filepath = gc_file, name = prefix)
  RNAGenebodyCoveragePlot(filepath = genecov_file, name = prefix)
  RNAFilteringPlot(filepath = countgene_file, UMI.cutoff = 1000, gene.number.cutoff = 500, name = prefix)
}else if(platform == "10x-genomics"){
  RNAReadDistrPlot(bamstat.filepath = bamstat_file, readdistr.filepath = readdistr_file, name = prefix, platform = platform)
  RNAFilteringPlot(filepath = countgene_file, UMI.cutoff = 1000, gene.number.cutoff = 500, name = prefix)
}else{
  RNAFilteringPlot(filepath = countgene_file, UMI.cutoff = 1000, gene.number.cutoff = 500, name = prefix)
}

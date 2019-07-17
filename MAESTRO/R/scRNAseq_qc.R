argue = commandArgs(T)
cwd = getwd()
setwd(argue[5])
source("scRNAseq_function.R")
setwd(file.path(cwd,argue[4]))

gene_cov_file = argue[1]
count_gene_file = argue[2]
prefix = argue[3]
#ReadPlot(read_distr_file, prefix)
GenebodyCoverPlot(gene_cov_file, prefix)
CountGenePlot(count_gene_file, count_cutoff = 500, gene_cutoff = 200, prefix)

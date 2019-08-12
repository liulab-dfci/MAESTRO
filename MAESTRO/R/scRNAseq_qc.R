library(MAESTRO)
argue = commandArgs(T)

setwd(argue[4])

genecov_file = argue[1]
countgene_file = argue[2]
prefix = argue[3]

RNAGenebodyCoveragePlot(filepath = genecov_file, name = prefix)
RNAFilteringPlot(filepath = countgene_file, UMI.cutoff = 1000, gene.number.cutoff = 500, name = prefix)


# argue = commandArgs(T)
# cwd = getwd()
# setwd(argue[5])
# source("scRNAseq_function.R")
# setwd(file.path(cwd,argue[4]))

# gene_cov_file = argue[1]
# count_gene_file = argue[2]
# prefix = argue[3]
# #ReadPlot(read_distr_file, prefix)
# GenebodyCoverPlot(gene_cov_file, prefix)
# CountGenePlot(count_gene_file, count_cutoff = 500, gene_cutoff = 200, prefix)
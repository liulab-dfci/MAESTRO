library(MAESTRO)
argue = commandArgs(T)

setwd(argue[7])

qualcode_file = argue[1]
nvc_file = argue[2]
gc_file = argue[3]
genecov_file = argue[4]
countgene_file = argue[5]
prefix = argue[6]

RNAReadQualityPlot(filepath = qualcode_file, name = prefix)
RNANucleotidePlot(filepath = nvc_file, name = prefix)
RNAGCcontentPlot(filepath = gc_file, name = prefix)
RNAGenebodyCoveragePlot(filepath = genecov_file, name = prefix)
RNAFilteringPlot(filepath = countgene_file, UMI.cutoff = 1000, gene.number.cutoff = 500, name = prefix)

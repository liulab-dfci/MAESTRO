library(MAESTRO)
argue = commandArgs(T)

setwd(argue[9])

bamstat_file = argue[1]
readdistr_file = argue[2]
qualcode_file = argue[3]
nvc_file = argue[4]
gc_file = argue[5]
genecov_file = argue[6]
countgene_file = argue[7]
prefix = argue[8]

RNAReadDistrPlot(bamstat.filepath = bamstat_file, readdistr.filepath = readdistr_file, name = prefix)
RNAReadQualityPlot(filepath = qualcode_file, name = prefix)
RNANucleotidePlot(filepath = nvc_file, name = prefix)
RNAGCcontentPlot(filepath = gc_file, name = prefix)
RNAGenebodyCoveragePlot(filepath = genecov_file, name = prefix)
RNAFilteringPlot(filepath = countgene_file, UMI.cutoff = 1000, gene.number.cutoff = 500, name = prefix)

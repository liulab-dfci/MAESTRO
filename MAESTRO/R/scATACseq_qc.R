library(MAESTRO)

argue = commandArgs(T)
setwd(argue[6])
ATACMapPlot(filepath = argue[1], platform = argue[4], name = argue[5])
ATACFilteringPlot(filepath = argue[1], platform = argue[4], name = argue[5])
ATACReadDistrPlot(filepath = argue[2], name = argue[5])
ATACFragmentSizePlot(filepath = argue[3], name = argue[5])

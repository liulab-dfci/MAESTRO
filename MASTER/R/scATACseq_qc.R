argue = commandArgs(T)
source(paste0(argue[6], "/scATACseq_function.R"))
setwd(argue[5])

if(argue[3] == "microfluidic"){
  stat_mtx = read.table(argue[1], header = TRUE, sep = "\t", row.names = 1)
  frag_mtx = read.table(argue[2])
  MapPlot_micro(stat_mtx,argue[4])
  FripPlot_micro(stat_mtx,argue[4])
  FragPlot(frag_mtx,argue[4])
} else{
  stat_mtx = read.csv(argue[1])
  frag_mtx = read.csv(argue[2])
  MapPlot_10x(stat_mtx,argue[4])
  FripPlot_10x(stat_mtx,argue[4])
  FragPlot(frag_mtx,argue[4])
}

argue = commandArgs(T)
source(paste0(argue[6], "/scATACseq_function.R"))
setwd(argue[5])

count_table = read.table(argue[1], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
SeuratObj = PipelineSeurat(countMat = count_table, proj = argue[4], min.c = 10, min.p = 200, max.p = 50000, 
                        org=argue[3],normalization.method = NULL, dims.use = 1:15, res = 0.6, diff.p = TRUE)
RPmat = read.table(argue[2], header = TRUE, row.names = 1, sep = '\t', check.names=FALSE)
AnnotateRP(SeuratObj, RPmat, proj = argue[4])
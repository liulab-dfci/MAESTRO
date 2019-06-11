source("../R/scATACseq_function.R")

argue = commandArgs(T)
setwd(argue[4])

count_table = read.table(argue[1], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
Seurat = PipelineSeurat(countMat = count_table, proj = argue[3], min.c = 10, min.p = 200, max.p = 50000, 
                        org=argue[2],normalization.method = "LogNormalize", dims.use = 1:15, res = 0.6)
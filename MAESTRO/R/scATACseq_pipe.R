library(MAESTRO)
library(Seurat)

argue = commandArgs(T)
setwd(argue[5])
countmatrix = read.table(argue[1], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
RPmatrix = read.table(argue[2], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

result = ATACRunSeurat(inputMat = countmatrix, proj = argue[4], method = "LSI")
result$ATAC = ATACAnnotateCelltype(result$ATAC, RPmatrix, signatures = human.immune.CIBERSORT)
result.tfs = ATACAnnotateTranscriptionFactor(ATAC = result$ATAC, peaks = result$peaks, project = argue[4], giggle.path = argue[3])

library(MAESTRO)
library(Seurat)

argue = commandArgs(T)
setwd(argue[6])
countmatrix = read.table(argue[1], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
RPmatrix = read.table(argue[2], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

result = ATACRunSeurat(inputMat = countmatrix, project = argue[5], method = "LSI")
result$ATAC = ATACAnnotateCelltype(result$ATAC, RPmatrix, signatures = human.immune.CIBERSORT)
saveRDS(result$ATAC, paste0(argue[5], "_scATAC_Object.rds"))
result.tfs = ATACAnnotateTranscriptionFactor(ATAC = result$ATAC, peaks = result$peaks, project = argue[5], giggle.path = argue[4], organism = argue[3])

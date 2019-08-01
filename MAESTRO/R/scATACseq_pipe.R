library(MAESTRO)

argue = commandArgs(T)
setwd(argue[6])
countmatrix = read.table(argue[1], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
RPmatrix = read.table(argue[2], sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
if(argue[3] == "NULL"){
   data(human.immune.CIBERSORT)
   genesignature = human.immune.CIBERSORT}
else{
   genesignature = read.table(argue[3])
result = ATACRunSeurat(inputMat = countmatrix, proj = argue[5], method = "LSI")
result$ATAC = ATACAnnotateCelltype(result$ATAC, RPmatrix, genesignature)
result.tfs = ATACAnnotateTranscriptionFactor(ATAC = result$ATAC, peaks = result$peaks, project = argue[5], giggle.path = argue[4])

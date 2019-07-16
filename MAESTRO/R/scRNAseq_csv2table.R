argue = commandArgs(T)

count.table <- read.table(argue[1], header = TRUE, row.names = 1, sep = ',')

write.table(count.table, argue[2], col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
a <-commandArgs(T)
outname <- a[1]

data <- read.table(paste(outname,'qcGC.txt',sep="_"))
pdf(file=paste(outname,'Figure3_GC.pdf',sep="_"))
gc <- rep(data[,1]/data[,3], times=data[,2])
hist(gc,probability=T,breaks=100,xlab="GC content (%)",ylab="Density of Reads",border="blue",main="")
dev.off()
 
data <- read.table(paste(outname,'qcNVC.txt',sep="_"),row.names=1)
pdf(file=paste(outname,'Figure2_NVC.pdf',sep="_"))
position <- seq(1,ncol(data))
A_count <- data["A",]
C_count <- data["C",]
G_count <- data["G",]
T_count <- data["T",]

total= A_count + C_count + G_count + T_count
ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05
yn=min(A_count/total,C_count/total,G_count/total,T_count/total)
#pdf("test5m.NVC_plot.pdf")
plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")
lines(position,T_count/total,type="o",pch=20,col="red")
lines(position,G_count/total,type="o",pch=20,col="blue")
lines(position,C_count/total,type="o",pch=20,col="cyan")
legend(40,ym,legend=c("A","T","G","C"),col=c("dark green","red","blue","cyan"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan"))
dev.off()


data <- read.table(paste(outname,'qcQul.txt',sep="_"),row.names=1)
pdf(file=paste(outname,'Figure1_quality_heatmap.pdf',sep="_"))
mat <- as.matrix(data)
Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c('spline'))
heatmap(mat,Rowv=NA,Colv=NA,xlab="Position of Read",ylab="Phred Quality Score",labRow=row.names(data),labCol=seq(ncol(data)),col = Lab.palette(256),scale="none" )
dev.off()


data <- read.table(paste(outname,'qcGBcover.txt',sep="_"))
pdf(file=paste(outname,'Figure4_GBcover.pdf',sep="_"))

plot(data[,1],data[,2],xlab="Percentile of gene body (5'->3')",ylab='Read number',type='s')
dev.off()
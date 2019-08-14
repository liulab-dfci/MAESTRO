library(cluster)
a<-commandArgs(T)

cluster_mat <- read.table(a[1],row.names=1,header=T)
qc_mat <- read.table(a[2],row.names=1,header=T)
outname <- a[3]
use_cluster_mat <- cluster_mat[which(cluster_mat[,3]>0),]
unuse_cluster_mat <- cluster_mat[which(cluster_mat[,3]==0),]
silhouette_evaluate <- function(clusmat){
    C <- clusmat[,3]
    D <- dist(clusmat[,1:2])
    return(silhouette(C,D))
}
silScore <- silhouette_evaluate(use_cluster_mat)

if (is.na(silScore)){
    use_cluster_mat_sil <- cbind(use_cluster_mat,rep('NA',nrow(use_cluster_mat)))
}else{
    use_cluster_mat_sil <- cbind(use_cluster_mat,silScore[,'sil_width'])
    pdf(file=paste(outname,"_Figure12_silhouetteScore.pdf",sep=""))
    plot(silScore,main="")
    dev.off()
}
colnames(use_cluster_mat_sil)[4] <- 'silhouette' 
if(nrow(unuse_cluster_mat) == 0){
    cluster_mat_sil  <- use_cluster_mat_sil[rownames(cluster_mat),]
}else{
    unuse_cluster_mat_sil <- cbind(unuse_cluster_mat,rep('NA',nrow(unuse_cluster_mat)))
    colnames(unuse_cluster_mat_sil)[4] <- 'silhouette' 
    cluster_mat_sil <- rbind(use_cluster_mat_sil,unuse_cluster_mat_sil)[rownames(cluster_mat),]
}



use_qc_mat <- qc_mat[rownames(cluster_mat_sil),]
pdf(file=paste(outname,"_Figure13_totalUMIcolored.pdf",sep=""))
initCol <- c("#4575B4","#91BFDB","#E0F3F8","#FFFFBF","#FEE090","#FC8D59","#D73027")
l <- log10(use_qc_mat[,2])
#if (logsc) {
#  f <- l == 0
#  l <- log(l)
#  l[f] <- NA
#}
mi <- min(l,na.rm=TRUE)
ma <- max(l,na.rm=TRUE)
ColorRamp <- colorRampPalette(initCol)(100)
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
v <- round((l - mi)/(ma - mi)*99 + 1,0)
layout(matrix(data=c(1,2), nrow=2, ncol=1), heights=c(4,1))
par(mar = c(4,4,2.5,2))
plot(cluster_mat_sil[,1:2],xlab="t-SNE 1",ylab="t-SNE 2",main="",pch=20,cex=0,col="grey")
for ( k in 1:length(v) ){
    points(cluster_mat_sil[k,1],cluster_mat_sil[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
}
par(mar = c(4,4,1.5,2))
image(ColorLevels,1,
    matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),
    col=ColorRamp,
    xlab="log10(total UMI count)",ylab="",
    yaxt="n")
dev.off()

pdf(file=paste(outname,"_Figure14_intronRate_colored.pdf",sep=""))
initCol <- c("#4575B4","#91BFDB","#E0F3F8","#FFFFBF","#FEE090","#FC8D59","#D73027")
l <- use_qc_mat[,'intron']/(use_qc_mat[,'utr3']+use_qc_mat[,'utr5']+use_qc_mat[,'cds']+use_qc_mat[,'intron'])
#if (logsc) {
#  f <- l == 0
#  l <- log(l)
#  l[f] <- NA
#}
mi <- min(l,na.rm=TRUE)
ma <- max(l,na.rm=TRUE)
ColorRamp <- colorRampPalette(initCol)(100)
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
v <- round((l - mi)/(ma - mi)*99 + 1,0)
layout(matrix(data=c(1,2), nrow=2, ncol=1), heights=c(4,1))
par(mar = c(4,4,2.5,2))
plot(cluster_mat_sil[,1:2],xlab="t-SNE 1",ylab="t-SNE 2",main="",pch=20,cex=0,col="grey")
for ( k in 1:length(v) ){
    points(cluster_mat_sil[k,1],cluster_mat_sil[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
}
par(mar = c(4,4,1.5,2))
image(ColorLevels,1,
    matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),
    col=ColorRamp,
    xlab="intron rate",ylab="",
    yaxt="n")
dev.off()


outmat <- cbind(cluster_mat_sil,use_qc_mat)
write.table(outmat,file=paste(outname,'_features_clustercell.txt',sep=""),row.names=T,col.names=T,sep="\t",quote=F)



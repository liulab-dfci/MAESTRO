a<-commandArgs(T)
### read parameter, 

### here qcmatrix and expmatrix only contain cell_barcodes with >= covergncutoff(100, default) gene covered, 
### see step1
qcdata <- read.table(a[1],row.names=1,header=T)
expdata <- read.table(a[2],row.names=1,header=T)
outname <- a[3]
measure <- as.numeric(a[4])
measure_cutoff <- as.numeric(a[5])
remove_nondup <- as.numeric(a[6])
non_dup_cutoff <- as.numeric(a[7])
total_umi <- as.numeric(a[8])
cluster_qcdata <- a[9]
cluster_expdata <- a[10]
pngplot <- as.numeric(a[11])
### individual qc1 : duplicate rate of all cell_barcodes

duprate <- (qcdata[,'allreads']-qcdata[,'umi'])/qcdata[,'allreads']
names(duprate) <- row.names(qcdata)

#nodup_cell <- names(duprate)[which(duprate < non_dup_cutoff)]

pdf(file=paste(outname,"_Figure5_duprate.pdf",sep=""))
hist(duprate,breaks=200,border="blue",xlab="Reads duplicate rate\n(#mappable reads-#UMI)/#mappable reads",freq=T,main="")
abline(v=non_dup_cutoff,col="darkblue",lwd=3)
legend("top",lwd=3,col="darkblue",legend="cutoff of cell barcodes\nwith low duplcate rate",bty="n",cex=1.2)
dev.off()


### individual qc2 : umi v.s. covered gene number in each cell_barcodes
coverGN <- qcdata[,'coveredGN']
names(coverGN) <- row.names(qcdata)
umi <- qcdata[,'umi']
names(umi) <- row.names(qcdata)

### 2 column matrix for umi v.s. covered gene number
qc2matrix <- cbind(log10(umi),log10(coverGN))
### calculate umi count cutoff based on the umi count of 1000th(default) cell 
if (remove_nondup == 1){
    nodupumi <- umi[which(duprate >= non_dup_cutoff)]
    umicutoff <- nodupumi[order(nodupumi,decreasing=T)][min(length(nodupumi),measure_cutoff)]
}else{
    umicutoff <- umi[order(umi,decreasing=T)][min(length(umi),measure_cutoff)]
}
### qc2 plot , and output cluster cell name to "cluster_cell_name"
if (pngplot == 1){
    png(file=paste(outname,"_Figure7_umi_coverGN.png",sep=""))
}else{
    pdf(file=paste(outname,"_Figure7_umi_coverGN.pdf",sep="")) 
}
plot(qc2matrix,xlab="log10 #UMI",ylab="log10 #covered gene",pch=".",col="red")
### use covered gene number as cutoff
if(measure == 1){
    ### remove non-duplicate cell_barcode
    if(remove_nondup == 1){
        points(qc2matrix[which(duprate < non_dup_cutoff),],pch=".",col="blue")
        points(qc2matrix[which(coverGN >= measure_cutoff & duprate >= non_dup_cutoff),],pch=".",col="purple")
        legend("topleft",legend=c("normal duplicate rate cell barcodes","low duplicate rate cell barcodes","cell selected for clustering\n(normal dup + high covered)"),col=c("red","blue","purple"),bty="n",lwd=3,cex=1.2)
        # here cluster_cell_name is selected clustering cells 
        cluster_cell_name <- row.names(qc2matrix[which(coverGN >= measure_cutoff & duprate >= non_dup_cutoff),])
    ### not remove non-duplicate cell_barcode
    }else{
        points(qc2matrix[which(coverGN >= measure_cutoff),],pch=".",col="purple")
        legend("topleft",legend=c("cell barcode","STAMP barcodes\n(high covered gene number)"),col=c("red","purple"),bty="n",lwd=3,cex=1.2)
        cluster_cell_name <- row.names(qc2matrix[which(coverGN >= measure_cutoff),])
    }
    abline(h=log10(measure_cutoff))
### use umi count as cutoff
}else{ 
    if(remove_nondup == 1){
        points(qc2matrix[which(duprate < non_dup_cutoff),],pch=".",col="blue")
        points(qc2matrix[which(umi >= umicutoff & duprate >= non_dup_cutoff),],pch=".",col="purple")
        legend("topleft",legend=c("normal duprate cell barcodes","low duplicate cellbarcodes","cell selected for clustering\n(normal dup + highest umi counts)"),col=c("red","blue","purple"),bty="n",lwd=3,cex=1.2)
        cluster_cell_name <- row.names(qc2matrix[which(umi >= umicutoff & duprate >= non_dup_cutoff),])
    }else{
        points(qc2matrix[which(umi >= umicutoff),],pch=".",col="purple")
        legend("topleft",legend=c("cell barcode","STAMP barcodes\n(highest umi counts)"),col=c("red","purple"),bty="n",lwd=3,cex=1.2)
        cluster_cell_name <- row.names(qc2matrix[which(umi >= umicutoff),])
    }
    abline(v=log10(umicutoff))
}
dev.off()


### individual qc3 :  cumulative umi v.s. duprate or cumulative covered GN v.s. duprate
if (pngplot == 1){
    png(file=paste(outname,"_Figure6_cumUMI_duprate.png",sep=""))
}else{
    pdf(file=paste(outname,"_Figure6_cumUMI_duprate.pdf",sep=""))
}
### use covered gene as cutoff
### generate cumulative detected gene number curve, overwrite on the scatter plot of duplicate rate
### note that the cumulate gene number is smaller than directly add up detected gene number of each cell_barcode, because different cell_barcode may cover same gene subsets. 
if(measure == 1){
    cell_coverGN_order <- names(sort(coverGN,decreasing=T))
    cs_coverGN <- c()
    current_cumulate_gn <- (expdata[,cell_coverGN_order[1]] > 0)
    cs_coverGN <- c(cs_coverGN, length(which(current_cumulate_gn)))
    for(i in 2:length(cell_coverGN_order)){
        current_cumulate_gn <- ( current_cumulate_gn | (expdata[,cell_coverGN_order[i]] > 0) ) 
        cs_coverGN <- c(cs_coverGN, length(which(current_cumulate_gn)))
    }
    names(cs_coverGN) <- cell_coverGN_order
    cs_coverGN_mat <- cbind(seq(length(cs_coverGN)),(cs_coverGN)/nrow(expdata))
    duprate_mat <- cbind(seq(length(cs_coverGN)),(duprate[names(cs_coverGN)]))
    #smoothScatter(duprate_mat,transformation = function(x) x^2,xlab="Cell barcodes ordered by covered gene number from high to low",ylab="Reads duplicate rate")
    smoothScatter(duprate_mat,xlab="Cell barcodes ordered by covered gene number from high to low",ylab="Reads duplicate rate")

    abline_cutoff <- length(which(coverGN >= measure_cutoff))
    abline(v=abline_cutoff,col="purple",lwd=2,lty=2)
    par(new=T)
    plot(cs_coverGN_mat*100,type="l",lwd=2,col="red",xlab="",ylab="",main="",axes=F)
    axis(side=4)
    mtext("Detected gene%",4)
    legend("topright",legend=c("reads duplicate rate (left axis)","detected gene% (right axis)"),col=c("blue","red"),lwd=3,bty="n",cex=1.2)

### use umi count as cutoff
}else{
    # cumulative umi number, cell_barcodes sorted by umi number
    csumi <- cumsum(sort(umi,decreasing=T))
    csumi_mat <- cbind(seq(length(csumi)),(csumi)/total_umi)
    # csumi_matlog <- cbind(seq(length(csumi)),log10(csumi))
    duprate_mat <- cbind(seq(length(csumi)),(duprate[names(csumi)]))
 
    smoothScatter(duprate_mat,xlab="Cell barcodes ordered by UMI count from high to low",ylab="Reads duplicate rate")
    abline(v=measure_cutoff,col="purple",lwd=2,lty=2)
    par(new=T)
    #plot(csumi_matlog,type="l",lwd=2,col="purple",xlab="",ylab="",main="",axes=F)
    #par(new=T)
    plot(csumi_mat*100,type="l",lwd=2,col="red",xlab="",ylab="",main="",axes=F)
    axis(side=4)
    mtext("cumulative UMI%",4)
    legend("topright",legend=c("reads duplicate rate (left axis)","cumulative UMI% (right axis)"),col=c("blue","red"),lwd=3,bty="n",cex=1.2)
}
dev.off()

### qc4 : intron rate,  following steps are all based on selected real cell  

real_qcdata <- qcdata[cluster_cell_name,]
real_expdata <- expdata[,cluster_cell_name]

intron_exon_ratio <- (real_qcdata[,'intron']+1)/(real_qcdata[,'utr3']+real_qcdata[,'utr5']+real_qcdata[,'cds']+real_qcdata[,'intron']+1)

pdf(file=paste(outname,"_Figure9_intronrate.pdf",sep=""))
hist(intron_exon_ratio,n=200,border="blue",xlab="Intron rate",main="")
dev.off()

###  qc5 cover gene number distribution
pdf(file=paste(outname,"_Figure8_coverGN.pdf",sep=""))
hist(log10(coverGN[cluster_cell_name]),border="blue",n=200,main='',xlab="log10 #covered genes")
dev.off()

### output 

write.table(real_qcdata,file=cluster_qcdata,row.names=T,col.names=T,sep="\t",quote=F)
write.table(real_expdata,file=cluster_expdata,row.names=T,col.names=T,sep="\t",quote=F)


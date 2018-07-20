#############################
# call nucleosome positions #
#############################
################################################################
# call in wild type and met1 - identify differential occupancy #
################################################################
###########################################################################
# analyze wild type and met1 SPO11 in reciprocal differential nucleosomes #
###########################################################################

library(Rsamtools)
library(GenomicRanges)
library(nucleR)
library(IRanges)

bam.lap.pwd <- "/Users/IanHenderson/Desktop/Choi/BAM_NUCS_masters/"
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
genome <- readDNAStringSet(file="arabidopsis_genome")
names(genome) <- c("1","2","3","4","5","chloroplast","mitochondria")
chr.size.at <- width(genome)
names(chr.size.at) <- names(genome)
chr.size.at <- chr.size.at[1:5]
nuc.names <- c("WT_nuc","Arp6_nuc","Met1_nuc")

j=1	
nuc.bam <- readGAlignmentPairs(paste(bam.lap.pwd,nuc.names[j],"_k10_bt2_mapped_lowmiss_unique_both_sort.bam",sep=""))
nuc.ranges <- ranges(nuc.bam)
nuc.chrs <- as.data.frame(seqnames(nuc.bam))
nuc.chrs <- as.character(nuc.chrs[,1])
nuc.strand <- as.character(strand(nuc.bam))
nuc.ranged <- RangedData(ranges=nuc.ranges,space=nuc.chrs,strand=nuc.strand)
save(nuc.ranged,file="wt.nuc.ranged")
j=2	
nuc.bam <- readGAlignmentPairs(paste(bam.lap.pwd,nuc.names[j],"_k10_bt2_mapped_lowmiss_unique_both_sort.bam",sep=""))
nuc.ranges <- ranges(nuc.bam)
nuc.chrs <- as.data.frame(seqnames(nuc.bam))
nuc.chrs <- as.character(nuc.chrs[,1])
nuc.strand <- as.character(strand(nuc.bam))
nuc.ranged <- RangedData(ranges=nuc.ranges,space=nuc.chrs,strand=nuc.strand)
save(nuc.ranged,file="arp6.nuc.ranged")
j=3	
nuc.bam <- readGAlignmentPairs(paste(bam.lap.pwd,nuc.names[j],"_k10_bt2_mapped_lowmiss_unique_both_sort.bam",sep=""))
nuc.ranges <- ranges(nuc.bam)
nuc.chrs <- as.data.frame(seqnames(nuc.bam))
nuc.chrs <- as.character(nuc.chrs[,1])
nuc.strand <- as.character(strand(nuc.bam))
nuc.ranged <- RangedData(ranges=nuc.ranges,space=nuc.chrs,strand=nuc.strand)
save(nuc.ranged,file="met1.nuc.ranged")

load("wt.nuc.ranged")
chr.size <- NULL
for(i in 1:5){
print(i)
nuc.chr.ranged <- nuc.ranged[i]
chr.size <- c(chr.size,length(space(nuc.chr.ranged)))
}
lib.size <- sum(chr.size)
for(i in 1:5){
print(i)
nuc.chr.ranged <- nuc.ranged[i]
read.trim <- processReads(nuc.chr.ranged,type="paired",fragmentLen=200,trim=40)
wt.chr.cov <- coverage(read.trim,width=chr.size.at[i])
wt.chr.cov <- as.integer(wt.chr.cov[[1]])
wt.chr.cov <- wt.chr.cov/lib.size
wt.chr.cov <- wt.chr.cov*1000000
chr.coords <- seq(i,length(wt.chr.cov),by=1)
trim.fft <- filterFFT(wt.chr.cov,pcKeepComp=0.02)
peaks <- peakDetection(trim.fft,threshold="25%",score=TRUE,width=140)
write.table(peaks,file=paste(i,"_wt_peaks.txt",sep=""))
write.table(trim.fft,file=paste(i,"_wt.cov.fft.txt",sep=""))
}

load("met1.nuc.ranged")
chr.size <- NULL
for(i in 1:5){
print(i)
nuc.chr.ranged <- nuc.ranged[i]
chr.size <- c(chr.size,length(space(nuc.chr.ranged)))
}
lib.size <- sum(chr.size)
for(i in 1:5){
print(i)
nuc.chr.ranged <- nuc.ranged[i]
read.trim <- processReads(nuc.chr.ranged,type="paired",fragmentLen=200,trim=40)
wt.chr.cov <- coverage(read.trim,width=chr.size.at[i])
wt.chr.cov <- as.integer(wt.chr.cov[[1]])
wt.chr.cov <- wt.chr.cov/lib.size
wt.chr.cov <- wt.chr.cov*1000000
chr.coords <- seq(i,length(wt.chr.cov),by=1)
trim.fft <- filterFFT(wt.chr.cov,pcKeepComp=0.02)
peaks <- peakDetection(trim.fft,threshold="25%",score=TRUE,width=140)
write.table(peaks,file=paste(i,"_met1_peaks.txt",sep=""))
write.table(trim.fft,file=paste(i,"_met1.cov.fft.txt",sep=""))
}

########################################################################################################
# Calculate Col and met1 NUC values within Col highly positioned nucleosomes - identify differentials  #
########################################################################################################

bam.lap.pwd <- "/Users/IanHenderson/Desktop/Choi/BAM_NUCS_masters/"
library(GenomicRanges)
library(Biostrings)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(segmentSeq)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(seqLogo)
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
all.fft.peaks <- NULL
for(i in 1:5){
print(i)
fft.col	<- read.table(file=paste(i,"_wt.cov.fft.txt",sep=""))
fft.col <- fft.col[,1]
fft.met	<- read.table(file=paste(i,"_met1.cov.fft.txt",sep=""))
fft.met <- fft.met[,1]
print(i)
fft.coords <- seq(1,length(fft.col),by=1)
fft.ir.coords <- IRanges(start=fft.coords,width=1)
fft.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=fft.ir.coords)	
nuccol.peaks <- read.table(file=paste(i,"_wt_peaks.txt",sep=""))
print(dim(nuccol.peaks))
nuccol.peaks[,1] <- rep(i,length(nuccol.peaks[,1]))
if(length(which(nuccol.peaks[,2]<501)>0)){
  nuccol.peaks <- nuccol.peaks[-which(nuccol.peaks[,2]<501),]
}
if(length(which(nuccol.peaks[,2]>chr.lens[i]-501))){
  nuccol.peaks <- nuccol.peaks[-which(nuccol.peaks[,2]>chr.lens[i]-501),]
}
nuccol.peaks <- nuccol.peaks[which(nuccol.peaks$score>0.7),]
print(dim(nuccol.peaks))
nuccol.ir.coords <- IRanges(start=nuccol.peaks$start,end=nuccol.peaks$end)
nuccol.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=nuccol.ir.coords)	
print(i)
nuccol.overlaps <- getOverlaps(nuccol.gr.coords,fft.gr.coords,whichOverlaps=T)
nuccol.fftcol <- sapply(nuccol.overlaps,function(x) mean(fft.col[x]))	
nuccol.fftmet <- sapply(nuccol.overlaps,function(x) mean(fft.met[x]))	
print(i)
nuccol.peaks <- cbind(nuccol.peaks,nuccol.fftcol,nuccol.fftmet)
all.fft.peaks <- rbind(all.fft.peaks,nuccol.peaks)
}
write.table(all.fft.peaks,file="all.fft.peaks.txt")
all.fft.peaks <- read.table(file="all.fft.peaks.txt")

diff.nucs <- all.fft.peaks[,9]-all.ftt.peaks[,8]

plot(all.fft.peaks[,8],all.fft.peaks[,9],pch=".")
abline(0,1,col=2)

#############################################################
# SPO11 and NUCS in Col vs met1 at differential nucleosomes #
#############################################################

spo.lap.pwd <- "/Users/IanHenderson/Desktop/Choi/BAM_SPO11_masters/"
nuc.lap.pwd <- "/Users/IanHenderson/Desktop/Choi/BAM_NUCS_masters/"
library(GenomicRanges)
library(Biostrings)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(segmentSeq)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(seqLogo)
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
all.fft.peaks <- read.table(file="all.fft.peaks.txt")
all.peaks <- NULL
for(i in 1:5){
nuccol.peaks <- all.fft.peaks[which(all.fft.peaks[,1]==i),]
nuc.diff <- nuccol.peaks[,8]-nuccol.peaks[,9]
nuccol.peaks <- nuccol.peaks[which(nuc.diff>0.05),]
nuccol.peaks <- nuccol.peaks[order(nuccol.peaks[,9]),]
if(length(which(nuccol.peaks[,2]<501)>0)){
  nuccol.peaks <- nuccol.peaks[-which(nuccol.peaks[,2]<501),]
}
if(length(which(nuccol.peaks[,2]>chr.lens[i]-501))){
  nuccol.peaks <- nuccol.peaks[-which(nuccol.peaks[,2]>chr.lens[i]-501),]
}
all.peaks <- rbind(all.peaks,nuccol.peaks)	
}
write.table(all.peaks,file="nuc.peaks.forAJT.txt")

spo.lap.pwd <- "/Users/IanHenderson/Desktop/Choi/BAM_SPO11_masters/"
nuc.lap.pwd <- "/Users/IanHenderson/Desktop/Choi/BAM_NUCS_masters/"
library(GenomicRanges)
library(Biostrings)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(segmentSeq)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(seqLogo)
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
all.nuccol.spocol <- rep(0,times=1140)
all.rancol.spocol <- rep(0,times=1140)
all.nuccol.spomet <- rep(0,times=1140)
all.rancol.spomet <- rep(0,times=1140)
all.nuccol.nuccol <- rep(0,times=1140)
all.rancol.nuccol <- rep(0,times=1140)
all.nuccol.nucmet <- rep(0,times=1140)
all.rancol.nucmet <- rep(0,times=1140)
chr.nuccol.tots <- NULL
for(i in 1:5){	
nuc.col	<- read.table(file=paste(nuc.lap.pwd,"WT_nuc_chr",i,"_cov.norm.both.txt",sep=""))
nuc.col <- nuc.col[,1]
nuc.met	<- read.table(file=paste(nuc.lap.pwd,"Met1_nuc_chr",i,"_cov.norm.both.txt",sep=""))
nuc.met <- nuc.met[,1]
print(i)
nuc.coords <- seq(1,length(nuc.col),by=1)
nuc.ir.coords <- IRanges(start=nuc.coords,width=1)
nuc.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=nuc.ir.coords)	
print(i)
rpi.col	<- read.table(file=paste(spo.lap.pwd,"col_chr",i,"_cov.norm.both.txt",sep=""))
rpi.col <- rpi.col[,1]
rpi.met	<- read.table(file=paste(spo.lap.pwd,"met1_chr",i,"_cov.norm.both.txt",sep=""))
rpi.met <- rpi.met[,1]
print(i)
rpi.coords <- seq(1,length(rpi.col),by=1)
rpi.ir.coords <- IRanges(start=rpi.coords,width=1)
rpi.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=rpi.ir.coords)	
print(i)
nuccol.peaks <- all.fft.peaks[which(all.fft.peaks[,1]==i),]
nuc.diff <- nuccol.peaks[,8]-nuccol.peaks[,9]
nuccol.peaks <- nuccol.peaks[which(nuc.diff>0.05),]
nuccol.peaks <- nuccol.peaks[order(nuccol.peaks[,9]),]
if(length(which(nuccol.peaks[,2]<501)>0)){
  nuccol.peaks <- nuccol.peaks[-which(nuccol.peaks[,2]<501),]
}
if(length(which(nuccol.peaks[,2]>chr.lens[i]-501))){
  nuccol.peaks <- nuccol.peaks[-which(nuccol.peaks[,2]>chr.lens[i]-501),]
}
nuccol.ir.coords <- IRanges(start=nuccol.peaks$start-500,end=nuccol.peaks$end+500)
nuccol.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=nuccol.ir.coords)	
rancol <- round(runif(n=length(nuccol.peaks[,1]),min=1001,max=chr.lens[i]-1001))
rancol.ir.coords <- IRanges(start=rancol-570,end=rancol+569)
rancol.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=rancol.ir.coords)	
print(i)
nuccol.overlaps <- getOverlaps(nuccol.gr.coords,rpi.gr.coords,whichOverlaps=T)
nuccol.spocol <- sapply(nuccol.overlaps,function(x) cbind(rpi.col[x]))	
nuccol.spocol.sum <- rowSums(nuccol.spocol)
rancol.overlaps <- getOverlaps(rancol.gr.coords,rpi.gr.coords,whichOverlaps=T)
rancol.spocol <- sapply(rancol.overlaps,function(x) cbind(rpi.col[x]))	
rancol.spocol.sum <- rowSums(rancol.spocol)
print(i)
nuccol.spomet <- sapply(nuccol.overlaps,function(x) cbind(rpi.met[x]))	
nuccol.spomet.sum <- rowSums(nuccol.spomet)
rancol.overlaps <- getOverlaps(rancol.gr.coords,rpi.gr.coords,whichOverlaps=T)
rancol.spomet <- sapply(rancol.overlaps,function(x) cbind(rpi.met[x]))	
rancol.spomet.sum <- rowSums(rancol.spomet)
print(i)
nuccol.overlaps <- getOverlaps(nuccol.gr.coords,nuc.gr.coords,whichOverlaps=T)
nuccol.nuccol <- sapply(nuccol.overlaps,function(x) cbind(nuc.col[x]))	
nuccol.nuccol.sum <- rowSums(nuccol.nuccol)
rancol.overlaps <- getOverlaps(rancol.gr.coords,nuc.gr.coords,whichOverlaps=T)
rancol.nuccol <- sapply(rancol.overlaps,function(x) cbind(nuc.col[x]))	
rancol.nuccol.sum <- rowSums(rancol.nuccol)
print(i)
nuccol.overlaps <- getOverlaps(nuccol.gr.coords,nuc.gr.coords,whichOverlaps=T)
nuccol.nucmet <- sapply(nuccol.overlaps,function(x) cbind(nuc.met[x]))	
nuccol.nucmet.sum <- rowSums(nuccol.nucmet)
rancol.overlaps <- getOverlaps(rancol.gr.coords,nuc.gr.coords,whichOverlaps=T)
rancol.nucmet <- sapply(rancol.overlaps,function(x) cbind(nuc.met[x]))	
rancol.nucmet.sum <- rowSums(rancol.nucmet)
print(i)
all.nuccol.nuccol <- all.nuccol.nuccol+nuccol.nuccol.sum
all.rancol.nuccol <- all.rancol.nuccol+rancol.nuccol.sum
all.nuccol.nucmet <- all.nuccol.nucmet+nuccol.nucmet.sum
all.rancol.nucmet <- all.rancol.nucmet+rancol.nucmet.sum
print(i)	
all.nuccol.spocol <- all.nuccol.spocol+nuccol.spocol.sum
all.rancol.spocol <- all.rancol.spocol+rancol.spocol.sum
all.nuccol.spomet <- all.nuccol.spomet+nuccol.spomet.sum
all.rancol.spomet <- all.rancol.spomet+rancol.spomet.sum
print(i)
chr.nuccol.tots <- c(chr.nuccol.tots,length(nuccol.peaks[,1]))
}
nuccol.tot <- sum(chr.nuccol.tots)
all.nuccol.spocol <- all.nuccol.spocol/nuccol.tot
all.rancol.spocol <- all.rancol.spocol/nuccol.tot
all.nuccol.spomet <- all.nuccol.spomet/nuccol.tot
all.rancol.spomet <- all.rancol.spomet/nuccol.tot
all.nuccol.nuccol <- all.nuccol.nuccol/nuccol.tot
all.rancol.nuccol <- all.rancol.nuccol/nuccol.tot
all.nuccol.nucmet <- all.nuccol.nucmet/nuccol.tot
all.rancol.nucmet <- all.rancol.nucmet/nuccol.tot
write.table(all.nuccol.spocol,file="nuccol.diff.spocol.txt")
write.table(all.rancol.spocol,file="rancol.diff.spocol.txt")
write.table(all.nuccol.spomet,file="nuccol.diff.spomet.txt")
write.table(all.rancol.spomet,file="rancol.diff.spomet.txt")
write.table(all.nuccol.nuccol,file="nuccol.diff.nuccol.txt")
write.table(all.rancol.nuccol,file="rancol.diff.nuccol.txt")
write.table(all.nuccol.nucmet,file="nuccol.diff.nucmet.txt")
write.table(all.rancol.nucmet,file="rancol.diff.nucmet.txt")

par(mfcol=c(3,3))
par(mar=c(1.8,1.8,1.8,1.8))

xplot <- seq(-500,639,by=1)
spo.lim=c(5e-9,8.5e-09)
nuc.lim=c(7e-9,1.4e-8)

plot(xplot,all.nuccol.spocol,type="l",main="SPO11 Col @ diff.nuc",ylim=spo.lim,col=2)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)
plot(xplot,all.nuccol.spomet,type="l",main="SPO11 met1 @ diff.nuc",ylim=spo.lim,col=4)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)
plot(xplot,all.nuccol.spocol,type="l",main="overlay",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,all.nuccol.spomet,type="l",ylim=spo.lim,col=4,yaxt="n")
axis(side=4)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,all.nuccol.nuccol,type="l",main="NUC Col @ diff.nuc",ylim=nuc.lim,col=2)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)
plot(xplot,all.nuccol.nucmet,type="l",main="NUC met1 @ diff.nuc",ylim=nuc.lim,col=4)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)
plot(xplot,all.nuccol.nuccol,type="l",main="overlay",ylim=nuc.lim,col=2)
par(new=T)
plot(xplot,all.nuccol.nucmet,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,all.rancol.spocol,type="l",main="SPO11 NUC Col @ Ran",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,all.rancol.nuccol,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)

plot(xplot,all.rancol.spomet,type="l",main="SPO11 NUC met1 @ Ran",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,all.rancol.nucmet,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)

plot(xplot,all.rancol.spocol,type="l",main="overlay",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,all.rancol.nuccol,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)
par(new=T)
plot(xplot,all.rancol.spomet,type="l",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,all.rancol.nucmet,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)

##########################################################
# Methylation in Col vs met1 at differential nucleosomes #
##########################################################

library(Rsamtools)
library(GenomicRanges)
library(nucleR)
library(IRanges)
library(GenomicRanges)
library(Biostrings)
library(seqinr)
library(GenomicAlignments)
library(parallel)
library(segmentSeq)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(seqLogo)

wt.cg <- read.table(file="GSM980986_WT_rep2.cg.txt",header=T)
wt.chg <- read.table(file="GSM980986_WT_rep2.chg.txt",header=T)
wt.chh <- read.table(file="GSM980986_WT_rep2.chh.txt",header=T)
met.cg <- read.table(file="GSM981031_met1.cg.txt",header=T)
met.chg <- read.table(file="GSM981031_met1.chg.txt",header=T)
met.chh <- read.table(file="GSM981031_met1.chh.txt",header=T)
all.fft.peaks <- read.table(file="all.fft.peaks.txt")

chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")

all.nuccol.cgcol <- rep(0,times=1140)
all.rancol.cgcol <- rep(0,times=1140)
all.nuccol.chgcol <- rep(0,times=1140)
all.rancol.chgcol <- rep(0,times=1140)
all.nuccol.chhcol <- rep(0,times=1140)
all.rancol.chhcol <- rep(0,times=1140)
all.nuccol.cgmet <- rep(0,times=1140)
all.rancol.cgmet <- rep(0,times=1140)
all.nuccol.chgmet <- rep(0,times=1140)
all.rancol.chgmet <- rep(0,times=1140)
all.nuccol.chhmet <- rep(0,times=1140)
all.rancol.chhmet <- rep(0,times=1140)
chr.nuccol.tots <- NULL
for(i in 1:5){
print(i)
wt.cg.chr <- wt.cg[which(wt.cg[,1]==i),]
wt.chg.chr <- wt.chg[which(wt.chg[,1]==i),]
wt.chh.chr <- wt.chh[which(wt.chh[,1]==i),]
print(i)
met.cg.chr <- met.cg[which(met.cg[,1]==i),]
met.chg.chr <- met.chg[which(met.chg[,1]==i),]
met.chh.chr <- met.chh[which(met.chh[,1]==i),]
print(i)
wt.cg.vals <- rep(0,times=chr.lens[i])
wt.cg.vals[wt.cg.chr[,2]] <- wt.cg.chr[,3]
wt.chg.vals <- rep(0,times=chr.lens[i])
wt.chg.vals[wt.chg.chr[,2]] <- wt.chg.chr[,3]
wt.chh.vals <- rep(0,times=chr.lens[i])
wt.chh.vals[wt.chh.chr[,2]] <- wt.chh.chr[,3]
print(i)
met.cg.vals <- rep(0,times=chr.lens[i])
met.cg.vals[met.cg.chr[,2]] <- met.cg.chr[,3]
met.chg.vals <- rep(0,times=chr.lens[i])
met.chg.vals[met.chg.chr[,2]] <- met.chg.chr[,3]
met.chh.vals <- rep(0,times=chr.lens[i])
met.chh.vals[met.chh.chr[,2]] <- met.chh.chr[,3]
print(i)
meth.coords <- seq(1,chr.lens[i],by=1)
meth.ir.coords <- IRanges(start=meth.coords,width=1)
meth.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=meth.ir.coords)	
print(i)
nuccol.peaks <- all.fft.peaks[which(all.fft.peaks[,1]==i),]
nuc.diff <- nuccol.peaks[,8]-nuccol.peaks[,9]
nuccol.peaks <- nuccol.peaks[which(nuc.diff>0.05),]
nuccol.peaks <- nuccol.peaks[order(nuccol.peaks[,9]),]
if(length(which(nuccol.peaks[,2]<501)>0)){
  nuccol.peaks <- nuccol.peaks[-which(nuccol.peaks[,2]<501),]
}
if(length(which(nuccol.peaks[,2]>chr.lens[i]-501))){
  nuccol.peaks <- nuccol.peaks[-which(nuccol.peaks[,2]>chr.lens[i]-501),]
}
nuccol.ir.coords <- IRanges(start=nuccol.peaks$start-500,end=nuccol.peaks$end+500)
nuccol.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=nuccol.ir.coords)	
rancol <- round(runif(n=length(nuccol.peaks[,1]),min=1001,max=chr.lens[i]-1001))
rancol.ir.coords <- IRanges(start=rancol-570,end=rancol+569)
rancol.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=rancol.ir.coords)	
print(i)
nuccol.overlaps <- getOverlaps(nuccol.gr.coords,meth.gr.coords,whichOverlaps=T)
nuccol.cgcol <- sapply(nuccol.overlaps,function(x) cbind(wt.cg.vals[x]))	
nuccol.cgcol.sum <- rowSums(nuccol.cgcol)
nuccol.chgcol <- sapply(nuccol.overlaps,function(x) cbind(wt.chg.vals[x]))	
nuccol.chgcol.sum <- rowSums(nuccol.chgcol)
nuccol.chhcol <- sapply(nuccol.overlaps,function(x) cbind(wt.chh.vals[x]))	
nuccol.chhcol.sum <- rowSums(nuccol.chhcol)
rancol.overlaps <- getOverlaps(rancol.gr.coords,meth.gr.coords,whichOverlaps=T)
rancol.cgcol <- sapply(rancol.overlaps,function(x) cbind(wt.cg.vals[x]))	
rancol.cgcol.sum <- rowSums(rancol.cgcol)
rancol.chgcol <- sapply(rancol.overlaps,function(x) cbind(wt.chg.vals[x]))	
rancol.chgcol.sum <- rowSums(rancol.chgcol)
rancol.chhcol <- sapply(rancol.overlaps,function(x) cbind(wt.chh.vals[x]))	
rancol.chhcol.sum <- rowSums(rancol.chhcol)
print(i)
nuccol.cgmet <- sapply(nuccol.overlaps,function(x) cbind(met.cg.vals[x]))	
nuccol.cgmet.sum <- rowSums(nuccol.cgmet)
nuccol.chgmet <- sapply(nuccol.overlaps,function(x) cbind(met.chg.vals[x]))	
nuccol.chgmet.sum <- rowSums(nuccol.chgmet)
nuccol.chhmet <- sapply(nuccol.overlaps,function(x) cbind(met.chh.vals[x]))	
nuccol.chhmet.sum <- rowSums(nuccol.chhmet)
rancol.overlaps <- getOverlaps(rancol.gr.coords,meth.gr.coords,whichOverlaps=T)
rancol.cgmet <- sapply(rancol.overlaps,function(x) cbind(met.cg.vals[x]))	
rancol.cgmet.sum <- rowSums(rancol.cgmet)
rancol.chgmet <- sapply(rancol.overlaps,function(x) cbind(met.chg.vals[x]))	
rancol.chgmet.sum <- rowSums(rancol.chgmet)
rancol.chhmet <- sapply(rancol.overlaps,function(x) cbind(met.chh.vals[x]))	
rancol.chhmet.sum <- rowSums(rancol.chhmet)
print(i)
all.nuccol.cgcol <- all.nuccol.cgcol+nuccol.cgcol.sum
all.rancol.cgcol <- all.rancol.cgcol+rancol.cgcol.sum
all.nuccol.cgmet <- all.nuccol.cgmet+nuccol.cgmet.sum
all.rancol.cgmet <- all.rancol.cgmet+nuccol.cgmet.sum
all.nuccol.chgcol <- all.nuccol.chgcol+nuccol.chgcol.sum
all.rancol.chgcol <- all.rancol.chgcol+rancol.chgcol.sum
all.nuccol.chgmet <- all.nuccol.chgmet+nuccol.chgmet.sum
all.rancol.chgmet <- all.rancol.chgmet+rancol.chgmet.sum
all.nuccol.chhcol <- all.nuccol.chhcol+nuccol.chhcol.sum
all.rancol.chhcol <- all.rancol.chhcol+rancol.chgcol.sum
all.nuccol.chhmet <- all.nuccol.chhmet+nuccol.chhmet.sum
all.rancol.chhmet <- all.rancol.chhmet+rancol.chhmet.sum
chr.nuccol.tots <- c(chr.nuccol.tots,length(nuccol.peaks[,1]))
}
nuccol.tot <- sum(chr.nuccol.tots)
all.nuccol.cgcol <- all.nuccol.cgcol/nuccol.tot
all.rancol.cgcol <- all.rancol.cgcol/nuccol.tot
all.nuccol.chgcol <- all.nuccol.chgcol/nuccol.tot
all.rancol.chgcol <- all.rancol.chgcol/nuccol.tot
all.nuccol.chhcol <- all.nuccol.chhcol/nuccol.tot
all.rancol.chhcol <- all.rancol.chhcol/nuccol.tot
all.nuccol.cgmet <- all.nuccol.cgmet/nuccol.tot
all.rancol.cgmet <- all.rancol.cgmet/nuccol.tot
all.nuccol.chgmet <- all.nuccol.chgmet/nuccol.tot
all.rancol.chgmet <- all.rancol.chgmet/nuccol.tot
all.nuccol.chhmet <- all.nuccol.chhmet/nuccol.tot
all.rancol.chhmet <- all.rancol.chhmet/nuccol.tot
write.table(all.nuccol.cgcol,file="nuccol.diff.cgcol.txt")
write.table(all.rancol.cgcol,file="rancol.diff.cgcol.txt")
write.table(all.nuccol.chgcol,file="nuccol.diff.chgcol.txt")
write.table(all.rancol.chgcol,file="rancol.diff.chgcol.txt")
write.table(all.nuccol.chhcol,file="nuccol.diff.chhcol.txt")
write.table(all.rancol.chhcol,file="rancol.diff.chhcol.txt")
write.table(all.nuccol.cgmet,file="nuccol.diff.cgmet.txt")
write.table(all.rancol.cgmet,file="rancol.diff.cgmet.txt")
write.table(all.nuccol.chgmet,file="nuccol.diff.chgmet.txt")
write.table(all.rancol.chgmet,file="rancol.diff.chgmet.txt")
write.table(all.nuccol.chhmet,file="nuccol.diff.chhmet.txt")
write.table(all.rancol.chhmet,file="rancol.diff.chhmet.txt")

test <- seq(1,500,by=1)
j=50
ma <- rep(1,test[j])/test[j]
filt.nuccol.cgcol <- filter(all.nuccol.cgcol,ma)
filt.nuccol.chgcol <- filter(all.nuccol.chgcol,ma)
filt.nuccol.chhcol <- filter(all.nuccol.chhcol,ma)
filt.nuccol.cgmet <- filter(all.nuccol.cgmet,ma)
filt.nuccol.chgmet <- filter(all.nuccol.chgmet,ma)
filt.nuccol.chhmet <- filter(all.nuccol.chhmet,ma)
filt.rancol.cgcol <- filter(all.rancol.cgcol,ma)
filt.rancol.chgcol <- filter(all.rancol.chgcol,ma)
filt.rancol.chhcol <- filter(all.rancol.chhcol,ma)
filt.rancol.cgmet <- filter(all.rancol.cgmet,ma)
filt.rancol.chgmet <- filter(all.rancol.chgmet,ma)
filt.rancol.chhmet <- filter(all.rancol.chhmet,ma)

nuccol.spocol <- read.table(file="nuccol.diff.spocol.txt")
nuccol.spocol <- nuccol.spocol[,1]
rancol.spocol <- read.table(file="rancol.diff.spocol.txt")
rancol.spocol <- rancol.spocol[,1]
nuccol.spomet <- read.table(file="nuccol.diff.spomet.txt")
nuccol.spomet <- nuccol.spomet[,1]
rancol.spomet <- read.table(file="rancol.diff.spomet.txt")
rancol.spomet <- rancol.spomet[,1]
nuccol.nuccol <- read.table(file="nuccol.diff.nuccol.txt")
nuccol.nuccol <- nuccol.nuccol[,1]
rancol.nuccol <- read.table(file="rancol.diff.nuccol.txt")
rancol.nuccol <- rancol.nuccol[,1]
nuccol.nucmet <- read.table(file="nuccol.diff.nucmet.txt")
nuccol.nucmet <- nuccol.nucmet[,1]
rancol.nucmet <- read.table(file="rancol.diff.nucmet.txt")
rancol.nucmet <- rancol.nucmet[,1]

meth.lim <- c(0,0.045)
spo.lim=c(5e-9,8.5e-09)
nuc.lim=c(7e-9,1.4e-8)

par(mfcol=c(2,4))
par(mar=c(1.8,1.8,1.8,1.8))
xplot <- seq(-500,639,by=1)

plot(xplot,filt.nuccol.cgcol,type="l",ylim=meth.lim,main="diff.nuc CG CHG CHH col",col=4)
par(new=T)
plot(xplot,filt.nuccol.chgcol,type="l",ylim=meth.lim,col="green")
par(new=T)
plot(xplot,filt.nuccol.chhcol,type="l",ylim=meth.lim,col=2)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,filt.nuccol.cgmet,type="l",ylim=meth.lim,main="diff.nuc CG CHG CHH met",col=4)
par(new=T)
plot(xplot,filt.nuccol.chgmet,type="l",ylim=meth.lim,col="green")
par(new=T)
plot(xplot,filt.nuccol.chhmet,type="l",ylim=meth.lim,col=2)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,filt.rancol.cgcol,type="l",ylim=meth.lim,main="diff.ran CG CHG CHH col",col=4)
par(new=T)
plot(xplot,filt.rancol.chgcol,type="l",ylim=meth.lim,col="green")
par(new=T)
plot(xplot,filt.rancol.chhcol,type="l",ylim=meth.lim,col=2)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,filt.rancol.cgmet,type="l",ylim=meth.lim,main="diff.ran CG CHG CHH met",col=4)
par(new=T)
plot(xplot,filt.rancol.chgmet,type="l",ylim=meth.lim,col="green")
par(new=T)
plot(xplot,filt.rancol.chhmet,type="l",ylim=meth.lim,col=2)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,nuccol.spocol,type="l",main="diff.nucs SPO11 NUCS Col",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,nuccol.nuccol,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,nuccol.spomet,type="l",main="diff.nucs SPO11 NUCS met1",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,nuccol.nucmet,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,rancol.spocol,type="l",main="diff.ran SPO11 NUCS Col",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,rancol.nuccol,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)

plot(xplot,rancol.spomet,type="l",main="diff.ran SPO11 NUCS met1",ylim=spo.lim,col=2)
par(new=T)
plot(xplot,rancol.nucmet,type="l",ylim=nuc.lim,col=4,yaxt="n")
axis(side=4)
abline(v=0,col=1,lty=2,lwd=0.5)
abline(v=140,col=1,lty=2,lwd=0.5)


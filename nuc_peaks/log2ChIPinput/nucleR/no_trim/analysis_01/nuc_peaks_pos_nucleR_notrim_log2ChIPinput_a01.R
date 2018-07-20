#############################################
# call and rank nuc peak positions         #
#############################################
# Rsamtools version 1.26.1
# GenomicRanges version 1.26.4
# nucleR version 2.6.0
# IRanges version 2.8.2

library(Rsamtools)
library(GenomicRanges)
library(nucleR)
library(IRanges)
library(doParallel)
library(parallel)

print(packageVersion("Rsamtools"))
print(packageVersion("GenomicRanges"))
print(packageVersion("nucleR"))
print(packageVersion("IRanges"))

#Definitions
bam.dir <- "/projects/ajt200/BAM_masters/nucleosomes/"
cov.dir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/"
out.dir <- "/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/no_trim/analysis_01/"
hist.dir <- "/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/no_trim/analysis_01/hist/"
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
#centromeres <- c(15086045,3607929,13587786,3956021,11725024)
#chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
#genome <- readDNAStringSet(file="/data/public_data/arabidopsis/TAIR_10/TAIR10_chr_all.fa")
#names(genome) <- c("1","2","3","4","5","chloroplast","mitochondria")
#chr.size.at <- width(genome)
#names(chr.size.at) <- names(genome)
#chr.size.at <- chr.size.at[1:5]
#lib.names <- c("REC8_ChIP", "REC8_input")

##Create REC8_ChIP ranged data object from bam
#lib.bam <- readGAlignmentPairs(paste0(bam.dir, lib.names[1], "_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
#lib.ranges <- ranges(lib.bam)
#lib.chrs <- as.data.frame(seqnames(lib.bam))
#lib.chrs <- as.character(lib.chrs[,1])
#lib.strand <- as.character(strand(lib.bam))
#REC8_ChIP.ranged <- RangedData(ranges=lib.ranges,space=lib.chrs,strand=lib.strand)
#save(REC8_ChIP.ranged, file=paste0(out.dir, "REC8_ChIP.ranged"))
#
##Create REC8_input ranged data object from bam
#lib.bam <- readGAlignmentPairs(paste0(bam.dir, lib.names[2], "_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
#lib.ranges <- ranges(lib.bam)
#lib.chrs <- as.data.frame(seqnames(lib.bam))
#lib.chrs <- as.character(lib.chrs[,1])
#lib.strand <- as.character(strand(lib.bam))
#REC8_input.ranged <- RangedData(ranges=lib.ranges,space=lib.chrs,strand=lib.strand)
#save(REC8_input.ranged, file=paste0(out.dir, "REC8_input.ranged"))
#
##Calculate REC8_ChIP library size
#load(file=paste0(out.dir, "REC8_ChIP.ranged"))
#chr.size <- NULL
#for(i in 1:5) {
#  print(i)
#  REC8_ChIP.chr.ranged <- REC8_ChIP.ranged[i]
#  chr.size <- c(chr.size,length(space(REC8_ChIP.chr.ranged)))
#}
#REC8_ChIP.lib.size <- sum(chr.size)
#
##Calculate REC8_input library size
#load(file=paste0(out.dir, "REC8_input.ranged"))
#chr.size <- NULL
#for(i in 1:5) {
#  print(i)
#  REC8_input.chr.ranged <- REC8_input.ranged[i]
#  chr.size <- c(chr.size,length(space(REC8_input.chr.ranged)))
#}
#REC8_input.lib.size <- sum(chr.size)
#
##Calculate REC8_ChIP normalised coverage
#for(i in 1:5) {
#  print(i)
#  REC8_ChIP.chr.ranged <- REC8_ChIP.ranged[i]
#  #removed "trim=40" (trim each read to 40 bp around the nucleosome dyad) as this step applies to nucleosome data 
#  REC8_ChIP.read.pair <- processReads(REC8_ChIP.chr.ranged, type="paired", fragmentLen=200)
#  REC8_ChIP.chr.cov <- coverage(REC8_ChIP.read.pair, width=chr.size.at[i])
#  REC8_ChIP.chr.cov <- as.integer(REC8_ChIP.chr.cov[[1]])
#  REC8_ChIP.chr.cov <- REC8_ChIP.chr.cov/REC8_ChIP.lib.size
#  REC8_ChIP.chr.cov <- REC8_ChIP.chr.cov*1000000
#  write.table(REC8_ChIP.chr.cov, file=paste0(out.dir, i, "_REC8_ChIP.chr.cov.txt"))
#}
#
##Calculate REC8_input normalised coverage
#for(i in 1:5) {
#  print(i)
#  REC8_input.chr.ranged <- REC8_input.ranged[i]
#  #removed "trim=40" (trim each read to 40 bp around the nucleosome dyad) as this step applies to nucleosome data 
#  REC8_input.read.pair <- processReads(REC8_input.chr.ranged, type="paired", fragmentLen=200)
#  REC8_input.chr.cov <- coverage(REC8_input.read.pair, width=chr.size.at[i])
#  REC8_input.chr.cov <- as.integer(REC8_input.chr.cov[[1]])
#  REC8_input.chr.cov <- REC8_input.chr.cov/REC8_input.lib.size
#  REC8_input.chr.cov <- REC8_input.chr.cov*1000000
#  write.table(REC8_input.chr.cov, file=paste0(out.dir, i, "_REC8_input.chr.cov.txt"))
#}

registerDoParallel(cores=5)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# load log2(nucleosomes/naked DNA) normalised coverage file and create FFT coverage and peaks files for each chromosome
cov <- read.table(file = paste0(cov.dir, "WT_log2nucNakedDNA_norm_allchrs_coverage_coord_tab.bed"))
foreach(i = 1:5, .combine = 'c') %dopar% {
  print(i)
  chr.cov <- cov[cov[,1] == chrs[i], 4]
  chr.coords <- seq(1, length(chr.cov), by = 1)
  cov.fft <- filterFFT(chr.cov, pcKeepComp = 0.02)
  #width = 140 (width = 200 for REC8 peak identification)
  peaks <- peakDetection(cov.fft, threshold="25%", score=TRUE, width=140)
  write.table(peaks, file=paste0(out.dir, i, "_log2nucNakedDNA_peaks.txt"))
  write.table(cov.fft, file=paste0(out.dir, i, "_log2nucNakedDNA.cov.fft.txt"))
}
rm(cov)

########################################################################################################
# Calculate mean FFT coverage values within peaks                                                 #
########################################################################################################

library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(segmentSeq)
library(seqLogo)
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
all.fft.peaks <- NULL
for(i in 1:5) {
  print(i)
  fft.nuc <- read.table(file=paste0(out.dir, i, "_log2nucNakedDNA.cov.fft.txt"))
  fft.nuc <- fft.nuc[,1]
  print(i)
  fft.coords <- seq(1,length(fft.nuc),by=1)
  fft.ir.coords <- IRanges(start=fft.coords,width=1)
  fft.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=fft.ir.coords)
  nuc.peaks <- read.table(file=paste0(out.dir, i, "_log2nucNakedDNA_peaks.txt"))  
  print(dim(nuc.peaks))
  nuc.peaks[,1] <- rep(i,length(nuc.peaks[,1]))
  if(length(which(nuc.peaks[,2]<1001)>0)) {
     nuc.peaks <- nuc.peaks[-which(nuc.peaks[,2]<1001),]
  }
  if(length(which(nuc.peaks[,2]>chr.lens[i]-1001))) {
     nuc.peaks <- nuc.peaks[-which(nuc.peaks[,2]>chr.lens[i]-1001),]
  }
  print(dim(nuc.peaks))
  nuc.ir.coords <- IRanges(start=nuc.peaks$start,end=nuc.peaks$end)
  nuc.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=nuc.ir.coords)
  print(i)
  nuc.overlaps <- getOverlaps(nuc.gr.coords,fft.gr.coords,whichOverlaps=T)
  nuc.fft <- sapply(nuc.overlaps,function(x) mean(fft.nuc[x]))
  print(i)
  nuc.peaks <- cbind(nuc.peaks,nuc.fft)
  all.fft.peaks <- rbind(all.fft.peaks,nuc.peaks)
}
write.table(all.fft.peaks,file=paste0(out.dir, "all.fft.peaks.txt"))
all.fft.peaks <- read.table(file=paste0(out.dir, "all.fft.peaks.txt"))

head(all.fft.peaks)
#   space start  end width      score    score_w   score_h   nuc.fft

dim(all.fft.peaks)
#

pdf(paste0(hist.dir, "hist_nuc_peaks_score_scoreh_FFTcov_20-07-17_brks1000.pdf"))
par(mfrow = c(2,2))
###Plot histogram of nuc peak combined scores ("score" column)
hist(all.fft.peaks$score, breaks=1000, col="black", ann=FALSE, xaxt = "n")
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc peak combined score"), expression("No. of peaks")))
axis(1, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))

###Plot histogram of nuc peak height scores ("score_h" column)
hist(all.fft.peaks$score_h, breaks=1000, col="black", ann=FALSE, xaxt = "n")
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc peak height score"), expression("No. of peaks")))
axis(1, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))

###Plot histogram of mean nuc FFT coverage at peaks
hist(all.fft.peaks$nuc.fft, breaks=1000, col="black", ann=FALSE)
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc FFT coverage at peaks"), expression("No. of peaks")))
dev.off()


##############################################################################################################
# Create ranking of peaks using dplyr percent_rank() function on FFT coverage values                        #
##############################################################################################################

library(dplyr)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

all.fft.peaks.arms <- NULL
all.fft.peaks.peri <- NULL
for(i in 1:5) {
  chr.all.fft.peaks <- all.fft.peaks[all.fft.peaks$space == i,]
  chr.all.fft.peaks.arms <- chr.all.fft.peaks %>%
    filter(end < pericenStart[i] | start > pericenEnd[i])
  chr.all.fft.peaks.peri <- chr.all.fft.peaks %>%
    filter(end > pericenStart[i] & start < pericenEnd[i])
  all.fft.peaks.arms <- rbind(all.fft.peaks.arms, chr.all.fft.peaks.arms)
  all.fft.peaks.peri <- rbind(all.fft.peaks.peri, chr.all.fft.peaks.peri)
}
print("Unfiltered arm peaks:")
print(dim(all.fft.peaks.arms)[[1]])
#
print("Unfiltered pericentromeric and centromeric peaks:")
print(dim(all.fft.peaks.peri)[[1]])
#
print("Unfiltered peaks:")
print(dim(all.fft.peaks.arms)[[1]] + dim(all.fft.peaks.peri)[[1]])
#

#Peaks on arms
all.fft.peaks.arms.rank <- mutate(all.fft.peaks.arms, nuc.rank = percent_rank(nuc.fft))
write.table(all.fft.peaks.arms.rank, file = paste0(out.dir, "all.fft.peaks.arms.rank_nucfft.txt"))

#Peaks in centromeric and pericentromeric regions
all.fft.peaks.peri.rank <- mutate(all.fft.peaks.peri, nuc.rank = percent_rank(nuc.fft))
write.table(all.fft.peaks.peri.rank, file = paste0(out.dir, "all.fft.peaks.peri.rank_nucfft.txt"))

pdf(paste0(hist.dir, "hist_nuc_arm_peaks_score_scoreh_FFTcov_19-07-17_brks1000.pdf"))
par(mfrow = c(2,2))
###Plot histogram of nuc peak combined scores ("score" column) at arm peaks
hist(all.fft.peaks.arms$score, breaks=1000, col="black", ann=FALSE, xaxt = "n")
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc arm peak combined score"), expression("No. of peaks")))
axis(1, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))

###Plot histogram of nuc peak height scores ("score_h" column) at arm peaks
hist(all.fft.peaks.arms$score_h, breaks=1000, col="black", ann=FALSE, xaxt = "n")
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc arm peak height score"), expression("No. of peaks")))
axis(1, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))

###Plot histogram of mean nuc FFT coverage at peaks at arm peaks
hist(all.fft.peaks.arms$nuc.fft, breaks=1000, col="black", ann=FALSE)
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc FFT coverage at arm peaks"), expression("No. of peaks")))
dev.off()

pdf(paste0(hist.dir, "hist_nuc_peri_peaks_score_scoreh_FFTcov_19-07-17_brks1000.pdf"))
par(mfrow = c(2,2))
###Plot histogram of nuc peak combined scores ("score" column) at peri peaks
hist(all.fft.peaks.peri$score, breaks=1000, col="black", ann=FALSE, xaxt = "n")
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc peri peak combined score"), expression("No. of peaks")))
axis(1, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))

###Plot histogram of nuc peak height scores ("score_h" column) at peri peaks
hist(all.fft.peaks.peri$score_h, breaks=1000, col="black", ann=FALSE, xaxt = "n")
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc peri peak height score"), expression("No. of peaks")))
axis(1, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))

###Plot histogram of mean nuc FFT coverage at peaks at peri peaks
hist(all.fft.peaks.peri$nuc.fft, breaks=1000, col="black", ann=FALSE)
mtext(side=c(1, 2), line=c(2.5, 2.5), text=c(expression("nuc FFT coverage at peri peaks"), expression("No. of peaks")))
dev.off()


###########################################################################################
# Filter arm and peri peaks based on distribution peaks over "score" and "score_h" values #
###########################################################################################

all.fft.peaks.arms.S54 <- all.fft.peaks.arms[all.fft.peaks.arms$score > 0.54,]
all.fft.peaks.peri.S54 <- all.fft.peaks.peri[all.fft.peaks.peri$score > 0.54,]
all.fft.peaks.arms.S60 <- all.fft.peaks.arms[all.fft.peaks.arms$score > 0.60,]
all.fft.peaks.peri.S60 <- all.fft.peaks.peri[all.fft.peaks.peri$score > 0.60,]
all.fft.peaks.arms.S70 <- all.fft.peaks.arms[all.fft.peaks.arms$score > 0.70,]
all.fft.peaks.peri.S70 <- all.fft.peaks.peri[all.fft.peaks.peri$score > 0.70,]
all.fft.peaks.arms.SH975 <- all.fft.peaks.arms[all.fft.peaks.arms$score_h > 0.975,]
all.fft.peaks.peri.SH975 <- all.fft.peaks.peri[all.fft.peaks.peri$score_h > 0.975,]

write.table(all.fft.peaks.arms.S54, file = paste0(out.dir, "all.fft.peaks.arms.S54.txt"))
write.table(all.fft.peaks.arms.S60, file = paste0(out.dir, "all.fft.peaks.arms.S60.txt"))
write.table(all.fft.peaks.arms.S70, file = paste0(out.dir, "all.fft.peaks.arms.S70.txt"))
write.table(all.fft.peaks.arms.SH975, file = paste0(out.dir, "all.fft.peaks.arms.SH975.txt"))
write.table(all.fft.peaks.peri.S54, file = paste0(out.dir, "all.fft.peaks.peri.S54.txt"))
write.table(all.fft.peaks.peri.S60, file = paste0(out.dir, "all.fft.peaks.peri.S60.txt"))
write.table(all.fft.peaks.peri.S70, file = paste0(out.dir, "all.fft.peaks.peri.S70.txt"))
write.table(all.fft.peaks.peri.SH975, file = paste0(out.dir, "all.fft.peaks.peri.SH975.txt"))

###Create GRanges objects to merge overlapping peaks
#score > 0.54
armPeaksS54GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.arms.S54$space), ranges = IRanges(start = all.fft.peaks.arms.S54$start, end = all.fft.peaks.arms.S54$end), strand = "+")
periPeaksS54GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.peri.S54$space), ranges = IRanges(start = all.fft.peaks.peri.S54$start, end = all.fft.peaks.peri.S54$end), strand = "+")
print("Unmerged arm peaks with score > 0.54")
print(length(armPeaksS54GR))
print("Unmerged peri peaks with score > 0.54")
print(length(periPeaksS54GR))
print("Unmerged peaks with score > 0.54")
print(length(armPeaksS54GR)+length(periPeaksS54GR))
#score > 0.60
armPeaksS60GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.arms.S60$space), ranges = IRanges(start = all.fft.peaks.arms.S60$start, end = all.fft.peaks.arms.S60$end), strand = "+")
periPeaksS60GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.peri.S60$space), ranges = IRanges(start = all.fft.peaks.peri.S60$start, end = all.fft.peaks.peri.S60$end), strand = "+")
print("Unmerged arm peaks with score > 0.60")
print(length(armPeaksS60GR))
print("Unmerged peri peaks with score > 0.60")
print(length(periPeaksS60GR))
print("Unmerged peaks with score > 0.60")
print(length(armPeaksS60GR)+length(periPeaksS60GR))
#score > 0.70
armPeaksS70GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.arms.S70$space), ranges = IRanges(start = all.fft.peaks.arms.S70$start, end = all.fft.peaks.arms.S70$end), strand = "+")
periPeaksS70GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.peri.S70$space), ranges = IRanges(start = all.fft.peaks.peri.S70$start, end = all.fft.peaks.peri.S70$end), strand = "+")
print("Unmerged arm peaks with score > 0.70")
print(length(armPeaksS70GR))
print("Unmerged peri peaks with score > 0.70")
print(length(periPeaksS70GR))
print("Unmerged peaks with score > 0.70")
print(length(armPeaksS70GR)+length(periPeaksS70GR))
#score_h > 0.975
armPeaksSH975GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.arms.SH975$space), ranges = IRanges(start = all.fft.peaks.arms.SH975$start, end = all.fft.peaks.arms.SH975$end), strand = "+")
periPeaksSH975GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.peri.SH975$space), ranges = IRanges(start = all.fft.peaks.peri.SH975$start, end = all.fft.peaks.peri.SH975$end), strand = "+")
print("Unmerged arm peaks with score_h > 0.975")
print(length(armPeaksSH975GR))
print("Unmerged peri peaks with score_h > 0.975")
print(length(periPeaksSH975GR))
print("Unmerged peaks with score_h > 0.975")
print(length(armPeaksSH975GR)+length(periPeaksSH975GR))

##Merge overlapping peaks
armPeaksS54GRmerge <- reduce(armPeaksS54GR)
armPeaksS60GRmerge <- reduce(armPeaksS60GR)
armPeaksS70GRmerge <- reduce(armPeaksS70GR)
armPeaksSH975GRmerge <- reduce(armPeaksSH975GR)
periPeaksS54GRmerge <- reduce(periPeaksS54GR)
periPeaksS60GRmerge <- reduce(periPeaksS60GR)
periPeaksS70GRmerge <- reduce(periPeaksS70GR)
periPeaksSH975GRmerge <- reduce(periPeaksSH975GR)

save(armPeaksS54GRmerge, file = paste0(out.dir, "armPeaksS54GRmerge.RData"))
save(armPeaksS60GRmerge, file = paste0(out.dir, "armPeaksS60GRmerge.RData"))
save(armPeaksS70GRmerge, file = paste0(out.dir, "armPeaksS70GRmerge.RData"))
save(armPeaksSH975GRmerge, file = paste0(out.dir, "armPeaksSH975GRmerge.RData"))
save(periPeaksS54GRmerge, file = paste0(out.dir, "periPeaksS54GRmerge.RData"))
save(periPeaksS60GRmerge, file = paste0(out.dir, "periPeaksS60GRmerge.RData"))
save(periPeaksS70GRmerge, file = paste0(out.dir, "periPeaksS70GRmerge.RData"))
save(periPeaksSH975GRmerge, file = paste0(out.dir, "periPeaksSH975GRmerge.RData"))


sessionInfo()


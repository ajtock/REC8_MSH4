#############################
# call nucleosome positions #
#############################
# Rsamtools version 1.26.1
# GenomicRanges version 1.26.4
# nucleR version 2.6.0
# IRanges version 2.8.2

library(Rsamtools)
library(GenomicRanges)
library(nucleR)
library(segmentSeq)
library(doParallel)
library(parallel)

#Definitions
bam.dir <- "/projects/ajt200/BAM_masters/nucleosomes/"
bam.dir2 <- "/projects/ajt200/BAM_masters/WT_nakedDNA/"
#cov.dir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/"
out.dir <- "/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/"
hist.dir <- "/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/hist/"
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
genome <- readDNAStringSet(file="/data/public_data/arabidopsis/TAIR_10/TAIR10_chr_all.fa")
names(genome) <- c("1","2","3","4","5","chloroplast","mitochondria")
chr.size.at <- width(genome)
names(chr.size.at) <- names(genome)
chr.size.at <- chr.size.at[1:5]
lib.names <- c("WT_nuc", "WT_nakedDNA")

# Create nucleosomes ranged data object from bam	
lib.bam <- readGAlignmentPairs(paste0(bam.dir, lib.names[1], "_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
lib.ranges <- ranges(lib.bam)
lib.chrs <- as.data.frame(seqnames(lib.bam))
lib.chrs <- as.character(lib.chrs[,1])
lib.strand <- as.character(strand(lib.bam))
nuc.ranged <- RangedData(ranges = lib.ranges, space = lib.chrs, strand = lib.strand)
save(nuc.ranged, file = paste0(out.dir, "nuc.ranged.RData"))

# Create naked DNA ranged data object from bam
lib.bam <- readGAlignmentPairs(paste0(bam.dir2, lib.names[2], "_trim51_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
lib.ranges <- ranges(lib.bam)
lib.chrs <- as.data.frame(seqnames(lib.bam))
lib.chrs <- as.character(lib.chrs[,1])
lib.strand <- as.character(strand(lib.bam))
nakedDNA.ranged <- RangedData(ranges = lib.ranges, space = lib.chrs, strand = lib.strand)
save(nakedDNA.ranged, file = paste0(out.dir, "nakedDNA.ranged.RData"))

# Calculate nucleosomes library size
chr.size <- NULL
for(i in 1:5) {
  print(i)
  nuc.chr.ranged <- nuc.ranged[i]
  chr.size <- c(chr.size, length(space(nuc.chr.ranged)))
}
nuc.lib.size <- sum(chr.size)

# Calculate naked DNA library size
chr.size <- NULL
for(i in 1:5) {
  print(i)
  nakedDNA.chr.ranged <- nakedDNA.ranged[i]
  chr.size <- c(chr.size, length(space(nakedDNA.chr.ranged)))
}
nakedDNA.lib.size <- sum(chr.size)

registerDoParallel(cores=5)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# Calculate nucleosomes and naked DNA library-size-normalised coverage,
# and compute log2 fold coverage for FFT transformation and peak calling 
foreach(i = 1:5, .combine = 'c') %dopar% {
  print(i)
  nuc.chr.ranged <- nuc.ranged[i]
  nuc.read.pair <- processReads(nuc.chr.ranged, type = "paired", fragmentLen = 200, trim = 40)
  nuc.chr.cov <- coverage(nuc.read.pair, width = chr.size.at[i])
  nuc.chr.cov <- as.integer(nuc.chr.cov[[1]])
  nuc.chr.cov <- nuc.chr.cov/nuc.lib.size
  nuc.chr.cov <- nuc.chr.cov*1000000
  write.table(nuc.chr.cov, file = paste0(out.dir, i, "_nuc.chr.cov.txt"))
  print(i)
  nakedDNA.chr.ranged <- nakedDNA.ranged[i]
  nakedDNA.read.pair <- processReads(nakedDNA.chr.ranged, type = "paired", fragmentLen = 200, trim = 40)
  nakedDNA.chr.cov <- coverage(nakedDNA.read.pair, width = chr.size.at[i])
  nakedDNA.chr.cov <- as.integer(nakedDNA.chr.cov[[1]])
  nakedDNA.chr.cov <- nakedDNA.chr.cov/nakedDNA.lib.size
  nakedDNA.chr.cov <- nakedDNA.chr.cov*1000000
  write.table(nakedDNA.chr.cov, file = paste0(out.dir, i, "_nakedDNA.chr.cov.txt"))
  print(i)
  nuc.chr.cov.offset <- nuc.chr.cov+1
  nakedDNA.chr.cov.offset <- nakedDNA.chr.cov+1
  log2nucNakedDNA.chr.cov <- log2(nuc.chr.cov.offset/nakedDNA.chr.cov.offset)
  log2nucNakedDNA.chr.cov <- (log2nucNakedDNA.chr.cov-mean(log2nucNakedDNA.chr.cov, na.rm = T))/sd(log2nucNakedDNA.chr.cov, na.rm =  T)
  cov.fft <- filterFFT(log2nucNakedDNA.chr.cov, pcKeepComp = 0.02)
  peaks <- peakDetection(cov.fft, threshold = "25%", score = TRUE, width = 140)
  write.table(peaks, file = paste0(out.dir, i, "_log2nucNakedDNA_peaks.txt"))
  write.table(cov.fft, file=paste0(out.dir, i, "_log2nucNakedDNA.cov.fft.txt"))
}

########################################################################################################
# Calculate mean FFT coverage values within peaks                                                      #
########################################################################################################
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
#9      1  1308 1447   140 0.58840304 1.00000000 0.17680608 0.02624212
#10     1  1433 1572   140 0.36625808 0.29712682 0.43538935 0.60278705
#11     1  1487 1626   140 0.59654707 0.40207235 0.79102179 0.87854635
#12     1  1586 1725   140 0.05070546 0.00155342 0.09985749 0.30804682
#13     1  1648 1787   140 0.46312876 0.51424190 0.41201562 0.27925514
#14     1  1713 1852   140 0.47813638 0.54561123 0.41066154 0.28708027

dim(all.fft.peaks)
#[1] 892654      8

pdf(paste0(hist.dir, "hist_nuc_peaks_trim_score_scoreh_FFTcov_23-08-17_brks1000.pdf"))
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
# Separate arm and pericetromeric peaks                                                                      #
# Create ranking of peaks using dplyr percent_rank() function on FFT coverage values                         #
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
#[1] 635366
print("Unfiltered pericentromeric and centromeric peaks:")
print(dim(all.fft.peaks.peri)[[1]])
#[1] 257288
print("Unfiltered peaks:")
print(dim(all.fft.peaks.arms)[[1]] + dim(all.fft.peaks.peri)[[1]])
#[1] 892654

#Peaks on arms
all.fft.peaks.arms.rank <- mutate(all.fft.peaks.arms, nuc.rank = percent_rank(nuc.fft))
write.table(all.fft.peaks.arms.rank, file = paste0(out.dir, "all.fft.peaks.arms.rank_nucfft.txt"))

#Peaks in centromeric and pericentromeric regions
all.fft.peaks.peri.rank <- mutate(all.fft.peaks.peri, nuc.rank = percent_rank(nuc.fft))
write.table(all.fft.peaks.peri.rank, file = paste0(out.dir, "all.fft.peaks.peri.rank_nucfft.txt"))

pdf(paste0(hist.dir, "hist_nuc_arm_peaks_trim_score_scoreh_FFTcov_23-08-17_brks1000.pdf"))
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

pdf(paste0(hist.dir, "hist_nuc_peri_peaks_trim_score_scoreh_FFTcov_23-07-17_brks1000.pdf"))
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

all.fft.peaks.arms.S55 <- all.fft.peaks.arms[all.fft.peaks.arms$score > 0.55,]
all.fft.peaks.peri.S55 <- all.fft.peaks.peri[all.fft.peaks.peri$score > 0.55,]
all.fft.peaks.arms.S60 <- all.fft.peaks.arms[all.fft.peaks.arms$score > 0.60,]
all.fft.peaks.peri.S60 <- all.fft.peaks.peri[all.fft.peaks.peri$score > 0.60,]
all.fft.peaks.arms.S70 <- all.fft.peaks.arms[all.fft.peaks.arms$score > 0.70,]
all.fft.peaks.peri.S70 <- all.fft.peaks.peri[all.fft.peaks.peri$score > 0.70,]
all.fft.peaks.arms.SH99 <- all.fft.peaks.arms[all.fft.peaks.arms$score_h > 0.99,]
all.fft.peaks.peri.SH99 <- all.fft.peaks.peri[all.fft.peaks.peri$score_h > 0.99,]

write.table(all.fft.peaks.arms.S55, file = paste0(out.dir, "all.fft.peaks.arms.S55.txt"))
write.table(all.fft.peaks.arms.S60, file = paste0(out.dir, "all.fft.peaks.arms.S60.txt"))
write.table(all.fft.peaks.arms.S70, file = paste0(out.dir, "all.fft.peaks.arms.S70.txt"))
write.table(all.fft.peaks.arms.SH99, file = paste0(out.dir, "all.fft.peaks.arms.SH99.txt"))
write.table(all.fft.peaks.peri.S55, file = paste0(out.dir, "all.fft.peaks.peri.S55.txt"))
write.table(all.fft.peaks.peri.S60, file = paste0(out.dir, "all.fft.peaks.peri.S60.txt"))
write.table(all.fft.peaks.peri.S70, file = paste0(out.dir, "all.fft.peaks.peri.S70.txt"))
write.table(all.fft.peaks.peri.SH99, file = paste0(out.dir, "all.fft.peaks.peri.SH99.txt"))

###Create GRanges objects to merge overlapping peaks
#score > 0.55
armPeaksS55GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.arms.S55$space), ranges = IRanges(start = all.fft.peaks.arms.S55$start, end = all.fft.peaks.arms.S55$end), strand = "+")
periPeaksS55GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.peri.S55$space), ranges = IRanges(start = all.fft.peaks.peri.S55$start, end = all.fft.peaks.peri.S55$end), strand = "+")
print("Unmerged arm peaks with score > 0.55")
print(length(armPeaksS55GR))
print("Unmerged peri peaks with score > 0.55")
print(length(periPeaksS55GR))
print("Unmerged peaks with score > 0.55")
print(length(armPeaksS55GR)+length(periPeaksS55GR))
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
#score_h > 0.99
armPeaksSH99GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.arms.SH99$space), ranges = IRanges(start = all.fft.peaks.arms.SH99$start, end = all.fft.peaks.arms.SH99$end), strand = "+")
periPeaksSH99GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.peri.SH99$space), ranges = IRanges(start = all.fft.peaks.peri.SH99$start, end = all.fft.peaks.peri.SH99$end), strand = "+")
print("Unmerged arm peaks with score_h > 0.99")
print(length(armPeaksSH99GR))
print("Unmerged peri peaks with score_h > 0.99")
print(length(periPeaksSH99GR))
print("Unmerged peaks with score_h > 0.99")
print(length(armPeaksSH99GR)+length(periPeaksSH99GR))

##Merge overlapping peaks
armPeaksS55GRmerge <- reduce(armPeaksS55GR)
armPeaksS60GRmerge <- reduce(armPeaksS60GR)
armPeaksS70GRmerge <- reduce(armPeaksS70GR)
armPeaksSH99GRmerge <- reduce(armPeaksSH99GR)
periPeaksS55GRmerge <- reduce(periPeaksS55GR)
periPeaksS60GRmerge <- reduce(periPeaksS60GR)
periPeaksS70GRmerge <- reduce(periPeaksS70GR)
periPeaksSH99GRmerge <- reduce(periPeaksSH99GR)

save(armPeaksS55GRmerge, file = paste0(out.dir, "armPeaksS55GRmerge.RData"))
save(armPeaksS60GRmerge, file = paste0(out.dir, "armPeaksS60GRmerge.RData"))
save(armPeaksS70GRmerge, file = paste0(out.dir, "armPeaksS70GRmerge.RData"))
save(armPeaksSH99GRmerge, file = paste0(out.dir, "armPeaksSH99GRmerge.RData"))
save(periPeaksS55GRmerge, file = paste0(out.dir, "periPeaksS55GRmerge.RData"))
save(periPeaksS60GRmerge, file = paste0(out.dir, "periPeaksS60GRmerge.RData"))
save(periPeaksS70GRmerge, file = paste0(out.dir, "periPeaksS70GRmerge.RData"))
save(periPeaksSH99GRmerge, file = paste0(out.dir, "periPeaksSH99GRmerge.RData"))


sessionInfo()


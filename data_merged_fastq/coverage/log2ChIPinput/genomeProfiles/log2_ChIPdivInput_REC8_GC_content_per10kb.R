##############################################################################################
# Calculate %GC and %AT in 10-kb windows and generate genome plots vs. log2-transformed REC8 #
##############################################################################################

library(segmentSeq)
library(seqinr)
#inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/"
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
#inNames <- c("WT_SPO11-oligo_meanAllReps")
#ChIPnames <- c("WT_SPO11_oligos")
#inputnames <- c("WT_H3K4me3_input", "MSH4_input")
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrIndex <- c(1:5)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
windows <- c(10000)
winNames <- c("10kb")

genome <- read.fasta(file = "/data/public_data/arabidopsis/TAIR_10/TAIR10_chr_all.fa")
for(i in 1:5) {
genome[[i]] <- genome[[i]][1:length(genome[[i]])]
}
genome <- genome[1:5]
print(summary(genome))
GCcont <- function(x) (length(x[x == "g" | x == "c" | x == "s"])/length(x))*100
ATcont <- function(x) (length(x[x == "a" | x == "t" | x == "w"])/length(x))*100
#calculate GC content per chromosome
lapply(seq_along(chrs), function(x) {
  print(paste0(chrs[x], " GC content"))
  print(GCcont(genome[[x]]))
})
#calculate AT content per chromosome
lapply(seq_along(chrs), function(x) {
  print(paste0(chrs[x], " AT content"))
  print(ATcont(genome[[x]]))
})

################################
# make cumulative genomes      #
################################

sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

for(s in 1:length(windows)) {
  print(s)
  wins <- windows[s]
  print(winNames[s])
  cumGCwinDat <- NULL
  cumATwinDat <- NULL
  for(i in 1:5) {
    print(i)
    genomeChr <- genome[[i]]
    coords <- seq(1, length(genomeChr), by = 1)
    IRcoords <- IRanges(start = coords, width = 1)
    GRcoords <- GRanges(seqnames = chrs[i], strand = "+", ranges = IRcoords)
    seqWindows <- seq(1, chrLens[i], by = wins)
    seqWindows <- c(seqWindows, chrLens[i])
    cumWindows <- seqWindows + sumchr[i]
    windowsIRanges <- IRanges(start = seqWindows, width = wins)
    windowsGRanges <- GRanges(seqnames = chrs[i], strand = "+", ranges = windowsIRanges)
    overlaps <- getOverlaps(windowsGRanges, GRcoords, whichOverlaps = TRUE)
    GCwinVals <- sapply(overlaps, function(x) GCcont(genomeChr[x]))
    ATwinVals <- sapply(overlaps, function(x) ATcont(genomeChr[x]))
    GCwinDat <- cbind(cumWindows, GCwinVals)
    ATwinDat <- cbind(cumWindows, ATwinVals)
    cumGCwinDat <- rbind(cumGCwinDat, GCwinDat)
    cumATwinDat <- rbind(cumATwinDat, ATwinDat)
  }
  write.table(cumGCwinDat, file = paste0(outDir, "GC_content_genome_", winNames[s], ".txt"))
  write.table(cumATwinDat, file = paste0(outDir, "AT_content_genome_", winNames[s], ".txt"))
}


REC8_ChIP_file <- paste0(outDir, "REC8_ChIP_genome_norm_coverage_", winNames[1], ".txt")
REC8_input_file <- paste0(outDir, "REC8_input_genome_norm_coverage_", winNames[1], ".txt")
GC_content_file <- paste0(outDir, "GC_content_genome_", winNames[s], ".txt")
AT_content_file <- paste0(outDir, "AT_content_genome_", winNames[s], ".txt")

# plot log2-normalised REC8 (genome-scale), adjusted so that the genome-wide average and SD of log2 ratio scores are 0 and 1, respectively, and overlay with %GC and %AT
  REC8_ChIP <- read.table(REC8_ChIP_file, header = TRUE)
  REC8_input <- read.table(REC8_input_file, header = TRUE)
  REC8_norm <- log2(REC8_ChIP[,2]/REC8_input[,2])
  REC8_norm <- (REC8_norm-mean(REC8_norm, na.rm = T))/sd(REC8_norm, na.rm = T)
  #REC8_norm_noNA <- REC8_norm[!is.na(REC8_norm)]
  GC_content <- read.table(GC_content_file, header = TRUE)
  AT_content <- read.table(AT_content_file, header = TRUE)
  GC <- GC_content[,2]
  AT <- AT_content[,2]

test <- seq(1, 1000, by = 1)
j = 100
ma <- rep(1, test[j])/test[j]
filt_REC8_norm <- stats::filter(REC8_norm, ma)
which_na <- which(is.na(filt_REC8_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_REC8_norm[left_na[length(left_na)]+1]
filt_REC8_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_REC8_norm[right_na[1]-1]
filt_REC8_norm[right_na] <- right_val
filt_REC8_norm_noNA <- filt_REC8_norm[!is.na(filt_REC8_norm)]

filt_GC <- stats::filter(GC, ma)
which_na <- which(is.na(filt_GC) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_GC[left_na[length(left_na)]+1]
filt_GC[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_GC[right_na[1]-1]
filt_GC[right_na] <- right_val
filt_GC_noNA <- filt_GC[!is.na(filt_GC)]

filt_AT <- stats::filter(AT, ma)
which_na <- which(is.na(filt_AT) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_AT[left_na[length(left_na)]+1]
filt_AT[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_AT[right_na[1]-1]
filt_AT[right_na] <- right_val 
filt_AT_noNA <- filt_AT[!is.na(filt_AT)]

ymin_REC8 <- min(filt_REC8_norm_noNA)
ymax_REC8 <- max(filt_REC8_norm_noNA)
ymin_GC <- min(filt_GC_noNA)
ymax_GC <- max(filt_GC_noNA)
ymin_AT <- min(filt_AT_noNA)
ymax_AT <- max(filt_AT_noNA)

pdf(file = paste0(outDir, "REC8_log2norm_mean0_sd1_GC_content_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(GC_content[,1], filt_GC, type = "l", lwd = 1.5, col = "darkgreen",
     ylim = c(ymin_GC, ymax_GC),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_GC_noNA), lty = 2, lwd = 1, col = "darkgreen")
axis(side = 4, at = pretty(c(filt_GC)))
mtext(side = 4, line = 3, cex = 1, text = expression(paste("GC content (%)")), col = "darkgreen")
par(new = T)
plot(REC8_ChIP[,1], filt_REC8_norm, type = "l", lwd = 1.5, col = "red",
     ylim = c(ymin_REC8, ymax_REC8),
     main = paste0("Spearman's rho = ", round(cor(filt_REC8_norm, filt_GC, method = "spearman"), digits = 2)),
     xlab = "",
     ylab = "log2(REC8 ChIP/REC8 input)", col.lab = "red")
abline(v = sumchr, lty = 1, lwd = 1)
abline(v = centromeres, lty = 2, lwd = 1, col = "black")
abline(v = pericenStart, lty = 3, lwd = 2, col = "black")
abline(v = pericenEnd, lty = 3, lwd = 2, col = "black")
abline(h = mean(filt_REC8_norm_noNA), lty = 2, lwd = 1, col = "red")
dev.off()

pdf(file = paste0(outDir, "REC8_log2norm_mean0_sd1_AT_content_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(AT_content[,1], filt_AT, type = "l", lwd = 1.5, col = "dodgerblue4",
     ylim = c(ymin_AT, ymax_AT),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_AT_noNA), lty = 2, lwd = 1, col = "dodgerblue4")
axis(side = 4, at = c(54, 58, 62, 66), labels = c("54", "58", "62", "66"))
#axis(side = 4, at = pretty(c(filt_AT)))
mtext(side = 4, line = 3, cex = 1, text = expression(paste("AT content (%)")), col = "dodgerblue4")
par(new = T)
plot(REC8_ChIP[,1], filt_REC8_norm, type = "l", lwd = 1.5, col = "red",
     ylim = c(ymin_REC8, ymax_REC8),
     main = paste0("Spearman's rho = ", round(cor(filt_REC8_norm, filt_AT, method = "spearman"), digits = 2)),
     xlab = "",
     ylab = "log2(REC8 ChIP/REC8 input)", col.lab = "red")
abline(v = sumchr, lty = 1, lwd = 1)
abline(v = centromeres, lty = 2, lwd = 1, col = "black")
abline(v = pericenStart, lty = 3, lwd = 2, col = "black")
abline(v = pericenEnd, lty = 3, lwd = 2, col = "black")
abline(h = mean(filt_REC8_norm_noNA), lty = 2, lwd = 1, col = "red")
dev.off()


print(cor.test(filt_REC8_norm, filt_GC, method = "pearson"))
#
#	Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_GC
#t = 54.306, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.4308640 0.4596472
#sample estimates:
#      cor 
#0.4453707 
print(cor.test(filt_REC8_norm, filt_GC, method = "spearman"))
#
#	Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_GC
#S = 1.716e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.3922288 

print(cor.test(filt_REC8_norm, filt_AT, method = "pearson"))
#
#	Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_AT
#t = -67.18, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.5369744 -0.5109301
#sample estimates:
#       cor 
#-0.5240747 
print(cor.test(filt_REC8_norm, filt_AT, method = "spearman"))
#
#	Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_AT
#S = 4.0723e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho 
#-0.4423024 


sessionInfo()

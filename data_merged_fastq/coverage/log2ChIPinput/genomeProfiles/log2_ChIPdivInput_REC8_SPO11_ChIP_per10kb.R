#############################################################################################
# Calculate log2 ratio of ChIP divided by input coverage values in 100-kb windows           #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
inDir <- "/projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/"
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
inNames <- c("WT_SPO11_ChIP4_nakedDNA", "WT_SPO11_ChIP13_nakedDNA")
ChIPnames <- c("WT_SPO11_ChIP4", "WT_SPO11_ChIP13")
#inputnames <- c("WT_nakedDNA", "WT_nakedDNA")
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrIndex <- c(1:5)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
windows <- c(10000)
winNames <- c("10kb")

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
  for(k in 1:length(inNames)) {
    print(k)
    ChIPinput <- read.table(file = paste0(inDir, inNames[k], "_norm_allchrs_coverage_coord_tab.bed"))
    cumChIPwinDat <- NULL
    #cuminputwinDat <- NULL
    for(i in 1:5) {
      print(i)
      chrChIPinput <- ChIPinput[ChIPinput$V1 == chrs[i],]
      ChIPcov <- chrChIPinput[,4]
      #inputcov <- chrChIPinput[,5]
      covCoords <- seq(1, length(ChIPcov), by = 1)
      covIRcoords <- IRanges(start = covCoords, width = 1)
      covGRcoords <- GRanges(seqnames = chrs[i], strand = "+", ranges = covIRcoords)
      seqWindows <- seq(1, chrLens[i], by = wins)
      seqWindows <- c(seqWindows, chrLens[i])
      cumWindows <- seqWindows + sumchr[i]
      windowsIRanges <- IRanges(start = seqWindows, width = wins)
      windowsGRanges <- GRanges(seqnames = chrs[i], strand = "+", ranges = windowsIRanges)
      overlaps <- getOverlaps(windowsGRanges, covGRcoords, whichOverlaps = TRUE)
      ChIPcovWinVals <- sapply(overlaps, function(x) mean(ChIPcov[x])+1) #+1 is offset to avoid infinite values
      #inputcovWinVals <- sapply(overlaps, function(x) mean(inputcov[x])+1) #+1 is offset to avoid infinite values
      ChIPwinDat <- cbind(cumWindows, ChIPcovWinVals)
      #inputwinDat <- cbind(cumWindows, inputcovWinVals)
      cumChIPwinDat <- rbind(cumChIPwinDat, ChIPwinDat)
      #cuminputwinDat <- rbind(cuminputwinDat, inputwinDat)
      #write.table(ChIPwinDat, file = paste0(outDir, ChIPnames[k], "_chr", chrIndex[i], "_norm_coverage_", winNames[s], ".txt"))
      #write.table(inputwinDat, file = paste0(outDir, inputnames[k], "_chr", chrIndex[i], "_norm_coverage_", winNames[s], ".txt"))
    }
    write.table(cumChIPwinDat, file = paste0(outDir, ChIPnames[k], "_genome_norm_coverage_", winNames[s], ".txt"))
    #write.table(cuminputwinDat, file = paste0(outDir, inputnames[k], "_genome_norm_coverage_", winNames[s], ".txt")) 
  }
}

REC8_ChIP_file <- paste0(outDir, "REC8_ChIP_genome_norm_coverage_", winNames[1], ".txt")
REC8_input_file <- paste0(outDir, "REC8_input_genome_norm_coverage_", winNames[1], ".txt")
WT_SPO11_ChIP4_file <- paste0(outDir, "WT_SPO11_ChIP4_genome_norm_coverage_", winNames[1], ".txt")
WT_SPO11_ChIP13_file <- paste0(outDir, "WT_SPO11_ChIP13_genome_norm_coverage_", winNames[1], ".txt")
WT_nakedDNA_file <- paste0(outDir, "WT_nakedDNA_genome_norm_coverage_", winNames[1], ".txt")

# plot log2-normalised REC8 and SPO11_ChIP4 coverage profiles (genome-scale), adjusted so that the genome-wide average and SD of log2 ratio scores are 0 and 1, respectively
  REC8_ChIP <- read.table(REC8_ChIP_file, header = TRUE)
  REC8_input <- read.table(REC8_input_file, header = TRUE)
  WT_SPO11_ChIP4 <- read.table(WT_SPO11_ChIP4_file, header = TRUE)
  WT_SPO11_ChIP13 <- read.table(WT_SPO11_ChIP13_file, header = TRUE)
  WT_nakedDNA <- read.table(WT_nakedDNA_file, header = TRUE)
  REC8_norm <- log2(REC8_ChIP[,2]/REC8_input[,2])
  REC8_norm <- (REC8_norm-mean(REC8_norm, na.rm = T))/sd(REC8_norm, na.rm = T)
  #REC8_norm_noNA <- REC8_norm[!is.na(REC8_norm)]
  SPO11_ChIP4_norm <- log2(WT_SPO11_ChIP4[,2]/WT_nakedDNA[,2])
  SPO11_ChIP4_norm <- (SPO11_ChIP4_norm-mean(SPO11_ChIP4_norm, na.rm = T))/sd(SPO11_ChIP4_norm, na.rm = T)
  #SPO11_ChIP4_norm_noNA <- SPO11_ChIP4_norm[!is.na(SPO11_ChIP4_norm)]
  SPO11_ChIP13_norm <- log2(WT_SPO11_ChIP13[,2]/WT_nakedDNA[,2])
  SPO11_ChIP13_norm <- (SPO11_ChIP13_norm-mean(SPO11_ChIP13_norm, na.rm = T))/sd(SPO11_ChIP13_norm, na.rm = T)
  #SPO11_ChIP13_norm_noNA <- SPO11_ChIP13_norm[!is.na(SPO11_ChIP13_norm)]

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

filt_SPO11_ChIP4_norm <- stats::filter(SPO11_ChIP4_norm, ma)
which_na <- which(is.na(filt_SPO11_ChIP4_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_SPO11_ChIP4_norm[left_na[length(left_na)]+1]
filt_SPO11_ChIP4_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_SPO11_ChIP4_norm[right_na[1]-1]
filt_SPO11_ChIP4_norm[right_na] <- right_val
filt_SPO11_ChIP4_norm_noNA <- filt_SPO11_ChIP4_norm[!is.na(filt_SPO11_ChIP4_norm)]

filt_SPO11_ChIP13_norm <- stats::filter(SPO11_ChIP13_norm, ma)
which_na <- which(is.na(filt_SPO11_ChIP13_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_SPO11_ChIP13_norm[left_na[length(left_na)]+1]
filt_SPO11_ChIP13_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_SPO11_ChIP13_norm[right_na[1]-1]
filt_SPO11_ChIP13_norm[right_na] <- right_val
filt_SPO11_ChIP13_norm_noNA <- filt_SPO11_ChIP13_norm[!is.na(filt_SPO11_ChIP13_norm)]

ymin_REC8 <- min(filt_REC8_norm_noNA)
ymax_REC8 <- max(filt_REC8_norm_noNA)
ymin_SPO11_ChIP4 <- min(filt_SPO11_ChIP4_norm_noNA)
ymax_SPO11_ChIP4 <- max(filt_SPO11_ChIP4_norm_noNA)
ymin_SPO11_ChIP13 <- min(filt_SPO11_ChIP13_norm_noNA)
ymax_SPO11_ChIP13 <- max(filt_SPO11_ChIP13_norm_noNA)

pdf(file = paste0(outDir, "REC8_SPO11_ChIP4_log2norm_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(WT_SPO11_ChIP4[,1], filt_SPO11_ChIP4_norm, type = "l", lwd = 1.5, col = "darkcyan",
     ylim = c(ymin_SPO11_ChIP4, ymax_SPO11_ChIP4),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_SPO11_ChIP4_norm_noNA), lty = 2, lwd = 1, col = "darkcyan")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(SPO11-1 ChIP/naked DNA)")), col = "darkcyan")
par(new = T)
plot(REC8_ChIP[,1], filt_REC8_norm, type = "l", lwd = 1.5, col = "red",
     ylim = c(ymin_REC8, ymax_REC8),
     main = "",
     xlab = "",
     ylab = "log2(REC8 ChIP/REC8 input)", col.lab = "red")
abline(v = sumchr, lty = 1, lwd = 1)
abline(v = centromeres, lty = 2, lwd = 1, col = "black")
abline(v = pericenStart, lty = 3, lwd = 2, col = "black")
abline(v = pericenEnd, lty = 3, lwd = 2, col = "black")
abline(h = mean(filt_REC8_norm_noNA), lty = 2, lwd = 1, col = "red")
dev.off()

pdf(file = paste0(outDir, "REC8_SPO11_ChIP13_log2norm_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(WT_SPO11_ChIP13[,1], filt_SPO11_ChIP13_norm, type = "l", lwd = 1.5, col = "darkcyan",
     ylim = c(ymin_SPO11_ChIP13, ymax_SPO11_ChIP13),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_SPO11_ChIP13_norm_noNA), lty = 2, lwd = 1, col = "darkcyan")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(SPO11-1 ChIP/naked DNA)")), col = "darkcyan")
par(new = T)
plot(REC8_ChIP[,1], filt_REC8_norm, type = "l", lwd = 1.5, col = "red",
     ylim = c(ymin_REC8, ymax_REC8),
     main = "",
     xlab = "",
     ylab = "log2(REC8 ChIP/REC8 input)", col.lab = "red")
abline(v = sumchr, lty = 1, lwd = 1)
abline(v = centromeres, lty = 2, lwd = 1, col = "black")
abline(v = pericenStart, lty = 3, lwd = 2, col = "black")
abline(v = pericenEnd, lty = 3, lwd = 2, col = "black")
abline(h = mean(filt_REC8_norm_noNA), lty = 2, lwd = 1, col = "red")
dev.off()

print(cor.test(filt_REC8_norm, filt_SPO11_ChIP4_norm, method = "pearson"))
print(cor.test(filt_REC8_norm, filt_SPO11_ChIP4_norm, method = "spearman"))

print(cor.test(filt_REC8_norm, filt_SPO11_ChIP13_norm, method = "pearson"))
print(cor.test(filt_REC8_norm, filt_SPO11_ChIP13_norm, method = "spearman"))

sessionInfo()

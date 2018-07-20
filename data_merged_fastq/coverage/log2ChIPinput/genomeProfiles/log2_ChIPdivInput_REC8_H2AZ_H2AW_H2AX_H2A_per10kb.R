#############################################################################################
# Calculate log2 ratio of ChIP divided by input coverage values in 100-kb windows           #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
inDir <- "/projects/ajt200/BAM_masters/H2A/coverage/"
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
inNames <- c("H2AZ_ChIP_input", "H2AW_ChIP_input", "H2AX_ChIP_input", "H2A_ChIP_input")
ChIPnames <- c("H2AZ_ChIP", "H2AW_ChIP", "H2AX_ChIP", "H2A_ChIP")
inputnames <- c("H2A_input", "H2A_input", "H2A_input", "H2A_input")
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
    cuminputwinDat <- NULL
    for(i in 1:5) {
      print(i)
      chrChIPinput <- ChIPinput[ChIPinput$V1 == chrs[i],]
      ChIPcov <- chrChIPinput[,4]
      inputcov <- chrChIPinput[,5]
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
      inputcovWinVals <- sapply(overlaps, function(x) mean(inputcov[x])+1) #+1 is offset to avoid infinite values
      ChIPwinDat <- cbind(cumWindows, ChIPcovWinVals)
      inputwinDat <- cbind(cumWindows, inputcovWinVals)
      cumChIPwinDat <- rbind(cumChIPwinDat, ChIPwinDat)
      cuminputwinDat <- rbind(cuminputwinDat, inputwinDat)
      #write.table(ChIPwinDat, file = paste0(outDir, ChIPnames[k], "_chr", chrIndex[i], "_norm_coverage_", winNames[s], ".txt"))
      #write.table(inputwinDat, file = paste0(outDir, inputnames[k], "_chr", chrIndex[i], "_norm_coverage_", winNames[s], ".txt"))
    }
    write.table(cumChIPwinDat, file = paste0(outDir, ChIPnames[k], "_genome_norm_coverage_", winNames[s], ".txt"))
    write.table(cuminputwinDat, file = paste0(outDir, inputnames[k], "_genome_norm_coverage_", winNames[s], ".txt")) 
  }
}

REC8_ChIP_file <- paste0(outDir, "REC8_ChIP_genome_norm_coverage_", winNames[1], ".txt")
REC8_input_file <- paste0(outDir, "REC8_input_genome_norm_coverage_", winNames[1], ".txt")
H2AZ_ChIP_file <- paste0(outDir, "H2AZ_ChIP_genome_norm_coverage_", winNames[1], ".txt")
H2AW_ChIP_file <- paste0(outDir, "H2AW_ChIP_genome_norm_coverage_", winNames[1], ".txt")
H2AX_ChIP_file <- paste0(outDir, "H2AX_ChIP_genome_norm_coverage_", winNames[1], ".txt")
H2A_ChIP_file <- paste0(outDir, "H2A_ChIP_genome_norm_coverage_", winNames[1], ".txt")
H2A_input_file <- paste0(outDir, "H2A_input_genome_norm_coverage_", winNames[1], ".txt")

# plot log2-normalised REC8 and H2AZ, H2AW, H2AX and H2A coverage profiles (genome-scale), adjusted so that the genome-wide average and SD of log2 ratio scores are 0 and 1, respectively
  REC8_ChIP <- read.table(REC8_ChIP_file, header = TRUE)
  REC8_input <- read.table(REC8_input_file, header = TRUE)
  H2AZ_ChIP <- read.table(H2AZ_ChIP_file, header = TRUE)
  H2AW_ChIP <- read.table(H2AW_ChIP_file, header = TRUE)
  H2AX_ChIP <- read.table(H2AX_ChIP_file, header = TRUE)
  H2A_ChIP <- read.table(H2A_ChIP_file, header = TRUE)
  H2A_input <- read.table(H2A_input_file, header = TRUE)
  REC8_norm <- log2(REC8_ChIP[,2]/REC8_input[,2])
  REC8_norm <- (REC8_norm-mean(REC8_norm, na.rm = T))/sd(REC8_norm, na.rm = T)
  #REC8_norm_noNA <- REC8_norm[!is.na(REC8_norm)]
  H2AZ_norm <- log2(H2AZ_ChIP[,2]/H2A_input[,2])
  H2AZ_norm <- (H2AZ_norm-mean(H2AZ_norm, na.rm = T))/sd(H2AZ_norm, na.rm = T)
  #H2AZ_norm_noNA <- H2AZ_norm[!is.na(H2AZ_norm)]
  H2AW_norm <- log2(H2AW_ChIP[,2]/H2A_input[,2])
  H2AW_norm <- (H2AW_norm-mean(H2AW_norm, na.rm = T))/sd(H2AW_norm, na.rm = T)
  #H2AW_norm_noNA <- H2AW_norm[!is.na(H2AW_norm)]
  H2AX_norm <- log2(H2AX_ChIP[,2]/H2A_input[,2])
  H2AX_norm <- (H2AX_norm-mean(H2AX_norm, na.rm = T))/sd(H2AX_norm, na.rm = T)
  #H2AX_norm_noNA <- H2AX_norm[!is.na(H2AX_norm)]
  H2A_norm <- log2(H2A_ChIP[,2]/H2A_input[,2])
  H2A_norm <- (H2A_norm-mean(H2A_norm, na.rm = T))/sd(H2A_norm, na.rm = T)
  #H2A_norm_noNA <- H2A_norm[!is.na(H2A_norm)]

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

filt_H2AZ_norm <- stats::filter(H2AZ_norm, ma)
which_na <- which(is.na(filt_H2AZ_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AZ_norm[left_na[length(left_na)]+1]
filt_H2AZ_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AZ_norm[right_na[1]-1]
filt_H2AZ_norm[right_na] <- right_val
filt_H2AZ_norm_noNA <- filt_H2AZ_norm[!is.na(filt_H2AZ_norm)]

filt_H2AW_norm <- stats::filter(H2AW_norm, ma)
which_na <- which(is.na(filt_H2AW_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AW_norm[left_na[length(left_na)]+1]
filt_H2AW_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AW_norm[right_na[1]-1]
filt_H2AW_norm[right_na] <- right_val
filt_H2AW_norm_noNA <- filt_H2AW_norm[!is.na(filt_H2AW_norm)]

filt_H2AX_norm <- stats::filter(H2AX_norm, ma)
which_na <- which(is.na(filt_H2AX_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AX_norm[left_na[length(left_na)]+1]
filt_H2AX_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AX_norm[right_na[1]-1]
filt_H2AX_norm[right_na] <- right_val
filt_H2AX_norm_noNA <- filt_H2AX_norm[!is.na(filt_H2AX_norm)]

filt_H2A_norm <- stats::filter(H2A_norm, ma)
which_na <- which(is.na(filt_H2A_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2A_norm[left_na[length(left_na)]+1]
filt_H2A_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2A_norm[right_na[1]-1]
filt_H2A_norm[right_na] <- right_val
filt_H2A_norm_noNA <- filt_H2A_norm[!is.na(filt_H2A_norm)]

ymin_REC8 <- min(filt_REC8_norm_noNA)
ymax_REC8 <- max(filt_REC8_norm_noNA)
ymin_H2AZ <- min(filt_H2AZ_norm_noNA)
ymax_H2AZ <- max(filt_H2AZ_norm_noNA)
ymin_H2AW <- min(filt_H2AW_norm_noNA)
ymax_H2AW <- max(filt_H2AW_norm_noNA)
ymin_H2AX <- min(filt_H2AX_norm_noNA)
ymax_H2AX <- max(filt_H2AX_norm_noNA)
ymin_H2A <- min(filt_H2A_norm_noNA)
ymax_H2A <- max(filt_H2A_norm_noNA)

pdf(file = paste0(outDir, "REC8_H2AZ_log2norm_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(H2AZ_ChIP[,1], filt_H2AZ_norm, type = "l", lwd = 1.5, col = "black",
     ylim = c(ymin_H2AZ, ymax_H2AZ),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_H2AZ_norm_noNA), lty = 2, lwd = 1, col = "black")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(H2AZ ChIP/input)")), col = "black")
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

pdf(file = paste0(outDir, "REC8_H2AW_log2norm_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(H2AW_ChIP[,1], filt_H2AW_norm, type = "l", lwd = 1.5, col = "black",
     ylim = c(ymin_H2AW, ymax_H2AW),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_H2AW_norm_noNA), lty = 2, lwd = 1, col = "black")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(H2AW ChIP/input)")), col = "black")
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

pdf(file = paste0(outDir, "REC8_H2AX_log2norm_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(H2AX_ChIP[,1], filt_H2AX_norm, type = "l", lwd = 1.5, col = "black",
     ylim = c(ymin_H2AX, ymax_H2AX),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_H2AX_norm_noNA), lty = 2, lwd = 1, col = "black")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(H2AX ChIP/input)")), col = "black")
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

pdf(file = paste0(outDir, "REC8_H2A_log2norm_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(H2A_ChIP[,1], filt_H2A_norm, type = "l", lwd = 1.5, col = "black",
     ylim = c(ymin_H2A, ymax_H2A),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_H2A_norm_noNA), lty = 2, lwd = 1, col = "black")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(H2A ChIP/input)")), col = "black")
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


print(cor.test(filt_REC8_norm, filt_H2AZ_norm, method = "pearson"))
#        Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_H2AZ_norm
#t = -236.83, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.9112459 -0.9049517
#sample estimates:
#       cor
#-0.9081501
print(cor.test(filt_REC8_norm, filt_H2AZ_norm, method = "spearman"))
#        Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_H2AZ_norm
#S = 4.0029e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho
#-0.4177272

print(cor.test(filt_REC8_norm, filt_H2AW_norm, method = "pearson"))
#        Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_H2AW_norm
#t = 205.2, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8788059 0.8867285
#sample estimates:
#    cor
#0.88283
print(cor.test(filt_REC8_norm, filt_H2AW_norm, method = "spearman"))
#        Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_H2AW_norm
#S = 1.3909e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.5073728

print(cor.test(filt_REC8_norm, filt_H2AX_norm, method = "pearson"))
#        Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_H2AX_norm
#t = 37.927, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.3120484 0.3440860
#sample estimates:
#      cor
#0.3281615
print(cor.test(filt_REC8_norm, filt_H2AX_norm, method = "spearman"))
#        Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_H2AX_norm
#S = 1.325e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.5307339

print(cor.test(filt_REC8_norm, filt_H2A_norm, method = "pearson"))
#        Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_H2A_norm
#t = -97.355, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.6754321 -0.6554298
#sample estimates:
#       cor
#-0.6655505
print(cor.test(filt_REC8_norm, filt_H2A_norm, method = "spearman"))
#        Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_H2A_norm
#S = 3.417e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#-0.210208


sessionInfo()

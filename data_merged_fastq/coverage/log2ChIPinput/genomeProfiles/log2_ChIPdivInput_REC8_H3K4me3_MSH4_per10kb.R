#############################################################################################
# Calculate log2 ratio of ChIP divided by input coverage values in 100-kb windows           #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
inDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/"
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
inNames <- c("WT_H3K4me3_ChIP_input", "MSH4_ChIP_input")
ChIPnames <- c("WT_H3K4me3_ChIP", "MSH4_ChIP")
inputnames <- c("WT_H3K4me3_input", "MSH4_input")
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
WT_H3K4me3_ChIP_file <- paste0(outDir, "WT_H3K4me3_ChIP_genome_norm_coverage_", winNames[1], ".txt")
WT_H3K4me3_input_file <- paste0(outDir, "WT_H3K4me3_input_genome_norm_coverage_", winNames[1], ".txt")
MSH4_ChIP_file <- paste0(outDir, "MSH4_ChIP_genome_norm_coverage_", winNames[1], ".txt")
MSH4_input_file <- paste0(outDir, "MSH4_input_genome_norm_coverage_", winNames[1], ".txt")

# plot log2-normalised REC8 and H3K4me3 coverage profiles (genome-scale), adjusted so that the genome-wide average and SD of log2 ratio scores are 0 and 1, respectively
  REC8_ChIP <- read.table(REC8_ChIP_file, header = TRUE)
  REC8_input <- read.table(REC8_input_file, header = TRUE)
  WT_H3K4me3_ChIP <- read.table(WT_H3K4me3_ChIP_file, header = TRUE)
  WT_H3K4me3_input <- read.table(WT_H3K4me3_input_file, header = TRUE)
  MSH4_ChIP <- read.table(MSH4_ChIP_file, header = TRUE)
  MSH4_input <- read.table(MSH4_input_file, header = TRUE)  
  REC8_norm <- log2(REC8_ChIP[,2]/REC8_input[,2])
  REC8_norm <- (REC8_norm-mean(REC8_norm, na.rm = T))/sd(REC8_norm, na.rm = T)
  #REC8_norm_noNA <- REC8_norm[!is.na(REC8_norm)]
  H3K4me3_norm <- log2(WT_H3K4me3_ChIP[,2]/WT_H3K4me3_input[,2])
  H3K4me3_norm <- (H3K4me3_norm-mean(H3K4me3_norm, na.rm = T))/sd(H3K4me3_norm, na.rm = T)
  #H3K4me3_norm_noNA <- H3K4me3_norm[!is.na(H3K4me3_norm)]
  MSH4_norm <- log2(MSH4_ChIP[,2]/MSH4_input[,2])
  MSH4_norm <- (MSH4_norm-mean(MSH4_norm, na.rm = T))/sd(MSH4_norm, na.rm = T)
  #MSH4_norm_noNA <- MSH4_norm[!is.na(MSH4_norm)]

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

filt_H3K4me3_norm <- stats::filter(H3K4me3_norm, ma)
which_na <- which(is.na(filt_H3K4me3_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H3K4me3_norm[left_na[length(left_na)]+1]
filt_H3K4me3_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H3K4me3_norm[right_na[1]-1]
filt_H3K4me3_norm[right_na] <- right_val
filt_H3K4me3_norm_noNA <- filt_H3K4me3_norm[!is.na(filt_H3K4me3_norm)]

filt_MSH4_norm <- stats::filter(MSH4_norm, ma)
which_na <- which(is.na(filt_MSH4_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_MSH4_norm[left_na[length(left_na)]+1]
filt_MSH4_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_MSH4_norm[right_na[1]-1]
filt_MSH4_norm[right_na] <- right_val
filt_MSH4_norm_noNA <- filt_MSH4_norm[!is.na(filt_MSH4_norm)]

ymin_REC8 <- min(filt_REC8_norm_noNA)
ymax_REC8 <- max(filt_REC8_norm_noNA)
ymin_H3K4me3 <- min(filt_H3K4me3_norm_noNA)
ymax_H3K4me3 <- max(filt_H3K4me3_norm_noNA)
ymin_MSH4 <- min(filt_MSH4_norm_noNA)
ymax_MSH4 <- max(filt_MSH4_norm_noNA)

pdf(file = paste0(outDir, "REC8_H3K4me3_log2norm_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(WT_H3K4me3_ChIP[,1], filt_H3K4me3_norm, type = "l", lwd = 1.5, col = "navy",
     ylim = c(ymin_H3K4me3, ymax_H3K4me3),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_H3K4me3_norm_noNA), lty = 2, lwd = 1, col = "navy")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(H3K4me3 ChIP/H3K9me2 input)")), col = "navy")
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

pdf(file = paste0(outDir, "REC8_MSH4_log2norm_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(MSH4_ChIP[,1], filt_MSH4_norm, type = "l", lwd = 1.5, col = "green",
     ylim = c(ymin_MSH4, ymax_MSH4),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_MSH4_norm_noNA), lty = 2, lwd = 1, col = "green")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(MSH4 ChIP/MSH4 input)")), col = "green")
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

print(cor.test(filt_REC8_norm, filt_H3K4me3_norm, method = "pearson"))
#
#	Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_H3K4me3_norm
#t = -174.28, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.8524359 -0.8423152
#sample estimates:
#       cor 
#-0.8474525 
print(cor.test(filt_REC8_norm, filt_H3K4me3_norm, method = "spearman"))
#
#	Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_H3K4me3_norm
#S = 3.8587e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho 
#-0.3666387
#
#Warning message:
#In cor.test.default(filt_REC8_norm, filt_H3K4me3_norm, method = "spearman") :
#  Cannot compute exact p-value with ties

print(cor.test(filt_REC8_norm, filt_MSH4_norm, method = "pearson"))
#
#	Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_MSH4_norm
#t = 365.37, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9566439 0.9595876
#sample estimates:
#      cor 
#0.9581411 
print(cor.test(filt_REC8_norm, filt_MSH4_norm, method = "spearman"))
#
#	Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_MSH4_norm
#S = 7.6753e+10, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.7281646 
#
#Warning message:
#In cor.test.default(filt_REC8_norm, filt_MSH4_norm, method = "spearman") :
#  Cannot compute exact p-value with ties


sessionInfo()

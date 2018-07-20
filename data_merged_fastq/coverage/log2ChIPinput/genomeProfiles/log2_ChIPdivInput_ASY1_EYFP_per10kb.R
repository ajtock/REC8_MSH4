#############################################################################################
# Calculate log2 ratio of ChIP divided by input coverage values in 100-kb windows           #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
inDir <- "/projects/ajt200/BAM_masters/ASY1_AllanWest/coverage/"
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
inNames <- c("EYFP_control")
ChIPnames <- c("EYFP_control")
#inputnames <- c("ASY1_input")
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

ASY1_ChIP_file <- paste0(outDir, "ASY1_ChIP_genome_norm_coverage_", winNames[1], ".txt")
ASY1_input_file <- paste0(outDir, "ASY1_input_genome_norm_coverage_", winNames[1], ".txt")
EYFP_control_file <- paste0(outDir, "EYFP_control_genome_norm_coverage_", winNames[1], ".txt")

# plot log2-normalised ASY1 and EYFP coverage profiles (genome-scale), adjusted so that the genome-wide average and SD of log2 ratio scores are 0 and 1, respectively
  ASY1_ChIP <- read.table(ASY1_ChIP_file, header = TRUE)
  ASY1_input <- read.table(ASY1_input_file, header = TRUE)
  EYFP_control <- read.table(EYFP_control_file, header = TRUE)
  ASY1_norm <- log2(ASY1_ChIP[,2]/ASY1_input[,2])
  ASY1_norm <- (ASY1_norm-mean(ASY1_norm, na.rm = T))/sd(ASY1_norm, na.rm = T)
  #ASY1_norm_noNA <- ASY1_norm[!is.na(ASY1_norm)]
  EYFP_norm <- EYFP_control[,2]
  EYFP_norm <- (EYFP_norm-mean(EYFP_norm, na.rm = T))/sd(EYFP_norm, na.rm = T)
  #EYFP_norm_noNA <- EYFP_norm[!is.na(EYFP_norm)]

test <- seq(1, 1000, by = 1)
j = 100
ma <- rep(1, test[j])/test[j]

filt_ASY1_norm <- stats::filter(ASY1_norm, ma)
which_na <- which(is.na(filt_ASY1_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_ASY1_norm[left_na[length(left_na)]+1]
filt_ASY1_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_ASY1_norm[right_na[1]-1]
filt_ASY1_norm[right_na] <- right_val
filt_ASY1_norm_noNA <- filt_ASY1_norm[!is.na(filt_ASY1_norm)]

filt_EYFP_norm <- stats::filter(EYFP_norm, ma)
which_na <- which(is.na(filt_EYFP_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_EYFP_norm[left_na[length(left_na)]+1]
filt_EYFP_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_EYFP_norm[right_na[1]-1]
filt_EYFP_norm[right_na] <- right_val
filt_EYFP_norm_noNA <- filt_EYFP_norm[!is.na(filt_EYFP_norm)]

ymin_ASY1 <- min(filt_ASY1_norm_noNA)
ymax_ASY1 <- max(filt_ASY1_norm_noNA)
ymin_EYFP <- min(filt_EYFP_norm_noNA)
ymax_EYFP <- max(filt_EYFP_norm_noNA)

pdf(file = paste0(outDir, "ASY1_log2norm_EYFP_mean0_sd1_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(ASY1_ChIP[,1], filt_ASY1_norm, type = "l", lwd = 1.5, col = "navy",
     ylim = c(ymin_ASY1, ymax_ASY1),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_ASY1_norm_noNA), lty = 2, lwd = 1, col = "navy")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("log2(ASY1 ChIP/ASY1 input)")), col = "navy")
par(new = T)
plot(EYFP_control[,1], filt_EYFP_norm, type = "l", lwd = 1.5, col = "goldenrod1",
     ylim = c(ymin_EYFP, ymax_EYFP),
     main = paste0("Spearman's rho = ", round(cor(filt_EYFP_norm, filt_ASY1_norm, method = "spearman"), digits = 2)),
     xlab = "",
     ylab = "EYFP control", col.lab = "goldenrod1")
abline(v = sumchr, lty = 1, lwd = 1)
abline(v = centromeres, lty = 2, lwd = 1, col = "black")
abline(v = pericenStart, lty = 3, lwd = 2, col = "black")
abline(v = pericenEnd, lty = 3, lwd = 2, col = "black")
abline(h = mean(filt_EYFP_norm_noNA), lty = 2, lwd = 1, col = "goldenrod1")
dev.off()

print(cor.test(filt_EYFP_norm, filt_ASY1_norm, method = "pearson"))
#        Pearson's product-moment correlation
#
#data:  filt_EYFP_norm and filt_ASY1_norm
#t = 262.29, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9205255 0.9258285
#sample estimates:
#     cor
#0.923221
print(cor.test(filt_EYFP_norm, filt_ASY1_norm, method = "spearman"))
#        Spearman's rank correlation rho
#
#data:  filt_EYFP_norm and filt_ASY1_norm
#S = 1.5028e+10, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.9467749


sessionInfo()

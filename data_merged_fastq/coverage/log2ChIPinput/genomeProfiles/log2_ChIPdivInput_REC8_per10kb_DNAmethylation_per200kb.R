#############################################################################################
# Calculate log2 ratio of ChIP divided by input coverage values in 100-kb windows           #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
bedDir <- "/home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/"
inDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/"
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
#inNames <- c("REC8_ChIP_input", "WT_nuc_nakedDNA")
#ChIPnames <- c("REC8_ChIP", "WT_nuc")
#inputnames <- c("REC8_input", "WT_nakedDNA")
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

REC8_ChIP_file <- paste0(outDir, "REC8_ChIP_genome_norm_coverage_", winNames[1], ".txt")
REC8_input_file <- paste0(outDir, "REC8_input_genome_norm_coverage_", winNames[1], ".txt")

# plot log2-normalised REC8 and nucleosomes coverage profiles (genome-scale), adjusted so that the genome-wide average and SD of log2 ratio scores are 0 and 1, respectively
  REC8_ChIP <- read.table(REC8_ChIP_file, header = TRUE)
  REC8_input <- read.table(REC8_input_file, header = TRUE)
  REC8_norm <- log2(REC8_ChIP[,2]/REC8_input[,2])
  REC8_norm <- (REC8_norm-mean(REC8_norm, na.rm = T))/sd(REC8_norm, na.rm = T)
  #REC8_norm_noNA <- REC8_norm[!is.na(REC8_norm)]

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

ymin_REC8 <- min(filt_REC8_norm_noNA)
ymax_REC8 <- max(filt_REC8_norm_noNA)

####################################################################
# Make cumulative genes using GRanges object containing genes      #
####################################################################

GSM980986_WT_rep2_CG <- read.table(file = paste0(bedDir, "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"))
GSM980986_WT_rep2_CHG <- read.table(file = paste0(bedDir, "GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed"))
GSM980986_WT_rep2_CHH <- read.table(file = paste0(bedDir, "GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed"))

cum_meth_dat <- NULL
for(i in 1:5) {
  # define windows2 as GRanges object
  windows2 <- seq(1, chrLens[i], by = 200000)
  windows2 <- c(windows2, chrLens[i])
  cum_windows2 <- windows2 + sumchr[i]
  windows2_iranges <- IRanges(start = windows2, width = 200000)
  windows2_granges <- GRanges(seqnames = chrs[i], strand = "+", ranges = windows2_iranges)

  # calculate mean methylation levels (all contexts) within windows2
  wtCG <- GSM980986_WT_rep2_CG[which(GSM980986_WT_rep2_CG$V1 == paste0("chr", i)),]
  wtCG_ir_coords <- IRanges(start = wtCG$V2, width = 1)
  wtCG_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = wtCG_ir_coords)
  overlapswtCG <- getOverlaps(windows2_granges, wtCG_gr_coords, whichOverlaps = TRUE)
  wtCG_win_vals <- sapply(overlapswtCG, function(x) mean(as.numeric(wtCG$V4[x])))

  wtCHG <- GSM980986_WT_rep2_CHG[which(GSM980986_WT_rep2_CHG$V1 == paste0("chr", i)),]
  wtCHG_ir_coords <- IRanges(start = wtCHG$V2, width = 1)
  wtCHG_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = wtCHG_ir_coords)
  overlapswtCHG <- getOverlaps(windows2_granges, wtCHG_gr_coords, whichOverlaps = TRUE)
  wtCHG_win_vals <- sapply(overlapswtCHG, function(x) mean(as.numeric(wtCHG$V4[x])))

  wtCHH <- GSM980986_WT_rep2_CHH[which(GSM980986_WT_rep2_CHH$V1 == paste0("chr", i)),]
  wtCHH_ir_coords <- IRanges(start = wtCHH$V2, width = 1)
  wtCHH_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = wtCHH_ir_coords)
  overlapswtCHH <- getOverlaps(windows2_granges, wtCHH_gr_coords, whichOverlaps = TRUE)
  wtCHH_win_vals <- sapply(overlapswtCHH, function(x) mean(as.numeric(wtCHH$V4[x])))

  dat2 <- cbind(cum_windows2, wtCG_win_vals, wtCHG_win_vals, wtCHH_win_vals)
  dat2 <- cbind(dat2, rowMeans(dat2[,2:4]))
  colnames(dat2) <- c("cum_windows2", "wtCG_win_vals", "wtCHG_win_vals", "wtCHH_win_vals", "allcntxts_win_vals")
  cum_meth_dat <- rbind(cum_meth_dat, dat2)
}
write.table(cum_meth_dat, file = paste0(outDir, "wtmeth_genome_200kb.txt"))

test <- seq(1, 1000, by = 1)
j = 5
ma <- rep(1, test[j])/test[j]
filt_wtCG <- stats::filter(cum_meth_dat[,2], ma); filt_wtCG_n <- filt_wtCG[!is.na(filt_wtCG)]
filt_wtCHG <- stats::filter(cum_meth_dat[,3], ma); filt_wtCHG_n <- filt_wtCHG[!is.na(filt_wtCHG)]
filt_wtCHH <- stats::filter(cum_meth_dat[,4], ma); filt_wtCHH_n <- filt_wtCHH[!is.na(filt_wtCHH)]
filt_allcntxts <- stats::filter(cum_meth_dat[,5], ma); filt_allcntxts_n <- filt_allcntxts[!is.na(filt_allcntxts)]
filt_wtMeth <- cbind(cum_meth_dat[,1], filt_wtCG, filt_wtCHG, filt_wtCHH, filt_allcntxts)
write.table(filt_wtMeth, file = paste0(outDir, "filt_wtMeth_genome_200kb.txt"))

ymin <- min(filt_wtCG_n, filt_wtCHG_n, filt_wtCHH_n)
ymax <- max(filt_wtCG_n, filt_wtCHG_n, filt_wtCHH_n)
ymin2 <- min(filt_allcntxts_n)
ymax2 <- max(filt_allcntxts_n)

pdf(file = paste0(outDir, "REC8_log2norm_mean0_sd1_10kb_methylation_200kb_genomeplot.pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(filt_wtMeth[,1], filt_wtMeth[,2], type = "l", lwd = 1.5, col = "navy",
     ylim = c(ymin, ymax), ylab = "",
     xlab = "Coordinates (bp)", xaxt = "n", yaxt = "n")
abline(h = mean(filt_wtCG_n), lty = 2, lwd = 1, col = "navy")
par(new = TRUE)
plot(filt_wtMeth[,1], filt_wtMeth[,3], type = "l", lwd = 1.5, col = "blue",
     ylim = c(ymin, ymax), ylab = "",
     ann = F, xaxt = "n", yaxt = "n")
abline(h = mean(filt_wtCHG_n), lty = 2, lwd = 1, col = "blue")
par(new = TRUE)
plot(filt_wtMeth[,1], filt_wtMeth[,4], type = "l", lwd = 1.5, col = "steelblue1",
     ylim = c(ymin, ymax), ylab = "",
     ann = F, xaxt = "n", yaxt = "n")
abline(h = mean(filt_wtCHH_n), lty = 2, lwd = 1, col = "steelblue1")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = "DNA methylation proportion")
par(new = TRUE)
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
legend("topleft",
       legend = c("wt CG", "wt CHG", "wt CHH"),
       col = c("navy", "blue", "steelblue1"),
       text.col = c("navy", "blue", "steelblue1"),
       ncol = 1, cex = 0.6, lwd = 1.2, bty = "n")
dev.off()

pdf(file = paste0(outDir, "REC8_log2norm_mean0_sd1_10kb_meanMethylation_200kb_genomeplot.pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(filt_wtMeth[,1], filt_wtMeth[,5], type = "l", lwd = 1.5, col = "deepskyblue4",
     ylim = c(ymin2, ymax2), ylab = "",
     xlab = "Coordinates (bp)", xaxt = "n", yaxt = "n")
abline(h = mean(filt_allcntxts_n), lty = 2, lwd = 1, col = "deepskyblue4")
axis(side = 4)
mtext(side = 4, line = 3.3, cex = 1, text = expression(atop(paste("DNA methylation proportion"), "(all contexts)")), col = "deepskyblue4")
par(new = TRUE)
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


sessionInfo()


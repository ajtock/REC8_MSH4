#############################################################################################
# Calculate log2 ratio of ChIP divided by input coverage values in 100-kb windows           #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
library(genomation)
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

representative_genes_uniq <- system("ls /projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", intern=T)
genes <- readGeneric(representative_genes_uniq, header=TRUE, strand=4, meta.col=list(gene_model=5))
seqlevels(genes) <- sub("", "Chr", seqlevels(genes))
print("******genes******")
print(genes)

for(s in 1:length(windows)) {
  print(s)
  wins <- windows[s]
  print(winNames[s])
  cumGeneDat <- NULL
  for(i in 1:5) {
    # define windows as GRanges object
    seqWindows <- seq(1, chrLens[i], by = wins)
    seqWindows <- c(seqWindows, chrLens[i])
    cumWindows <- seqWindows + sumchr[i]
    windowsIRanges <- IRanges(start = seqWindows, width = wins)
    windowsGRanges <- GRanges(seqnames = chrs[i], strand = "+", ranges = windowsIRanges)

    # define and count genes within windows
    chrGenes <- genes[seqnames(genes) == chrs[i]]    
    winGenes <- countOverlaps(windowsGRanges, chrGenes)
    dat <- cbind(cumWindows, winGenes)
    cumGeneDat <- rbind(cumGeneDat, dat)
  }
write.table(cumGeneDat, file = paste0(outDir, "gene_density_genome_", winNames[s], ".txt"))
}

filt_cumGeneDat <- stats::filter(cumGeneDat[,2], ma)
which_na <- which(is.na(filt_cumGeneDat) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumGeneDat[left_na[length(left_na)]+1]
filt_cumGeneDat[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumGeneDat[right_na[1]-1]
filt_cumGeneDat[right_na] <- right_val
filt_cumGeneDat_noNA <- filt_cumGeneDat[!is.na(filt_cumGeneDat)]

ymin_genes <- min(filt_cumGeneDat_noNA)
ymax_genes <- max(filt_cumGeneDat_noNA)

pdf(file = paste0(outDir, "REC8_log2norm_mean0_sd1_gene_density_genomeplot_", winNames, ".pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(cumGeneDat[,1], filt_cumGeneDat, type = "l", lwd = 1.5, col = "black",
     ylim = c(ymin_genes, ymax_genes),
     xlab = "Coordinates (bp)",
     ylab = "", xaxt = "n", yaxt = "n")
abline(h = mean(filt_cumGeneDat_noNA), lty = 2, lwd = 1, col = "black")
axis(side = 4)
mtext(side = 4, line = 3, cex = 1, text = expression(paste("Gene density 10 kb"^{-1})), col = "black")
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

print(cor.test(filt_REC8_norm, filt_cumGeneDat, method = "pearson"))
#
#	Pearson's product-moment correlation
#
#data:  filt_REC8_norm and filt_cumGeneDat
#t = -169.9, df = 11919, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.8464497 -0.8359549
#sample estimates:
#       cor 
#-0.8412815 

print(cor.test(filt_REC8_norm, filt_cumGeneDat, method = "spearman"))
#
#	Spearman's rank correlation rho
#
#data:  filt_REC8_norm and filt_cumGeneDat
#S = 4.0793e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho 
#-0.4447738 
#
#Warning message:
#In cor.test.default(filt_REC8_norm, filt_cumGeneDat, method = "spearman") :
#  Cannot compute exact p-value with ties


sessionInfo()

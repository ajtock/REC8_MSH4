#############################################################################################
# Calculate log2 ratio of ChIP divided by input coverage values in 100-kb windows           #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
library(corrplot)
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
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
WT_nuc_file <- paste0(outDir, "WT_nuc_genome_norm_coverage_", winNames[1], ".txt")
WT_nakedDNA_file <- paste0(outDir, "WT_nakedDNA_genome_norm_coverage_", winNames[1], ".txt")
H2AZ_ChIP_file <- paste0(outDir, "H2AZ_ChIP_genome_norm_coverage_", winNames[1], ".txt")
H2A_input_file <- paste0(outDir, "H2A_input_genome_norm_coverage_", winNames[1], ".txt")
WT_H3K4me3_ChIP_file <- paste0(outDir, "WT_H3K4me3_ChIP_genome_norm_coverage_", winNames[1], ".txt")
WT_H3K4me3_input_file <- paste0(outDir, "WT_H3K4me3_input_genome_norm_coverage_", winNames[1], ".txt")
gene_density_file <- paste0(outDir, "gene_density_genome_", winNames[1], ".txt")
DNA_meth_file200 <- paste0(outDir, "wtmeth_genome_200kb.txt")
DNA_meth_file10 <- paste0(outDir, "wtmeth_genome_10kb.txt")
H2AW_ChIP_file <- paste0(outDir, "H2AW_ChIP_genome_norm_coverage_", winNames[1], ".txt")
WT_H3K9me2_ChIP_file <- paste0(outDir, "WT_H3K9me2_ChIP_genome_norm_coverage_", winNames[1], ".txt")
WT_H3K9me2_input_file <- paste0(outDir, "WT_H3K9me2_input_genome_norm_coverage_", winNames[1], ".txt")
TE_density_file <- paste0(outDir, "TE_density_genome_", winNames[1], ".txt")

WT_SPO11_ChIP4_file <- paste0(outDir, "WT_SPO11_ChIP4_genome_norm_coverage_", winNames[1], ".txt")
WT_SPO11_oligos_RPI1_file <- paste0(outDir, "WT_SPO11_oligos_RPI1_genome_norm_coverage_", winNames[1], ".txt")
WT_nakedDNA_T50R1_file <- paste0(outDir, "WT_nakedDNA_T50R1_genome_norm_coverage_", winNames[1], ".txt")
CO_density_file <- paste0(outDir, "WTCO_density_genome_", winNames[1], ".txt")
COReduced_density_file <- paste0(outDir, "WTCOReduced_density_genome_", winNames[1], ".txt")
MSH4_ChIP_file <- paste0(outDir, "MSH4_ChIP_genome_norm_coverage_", winNames[1], ".txt")
MSH4_input_file <- paste0(outDir, "MSH4_input_genome_norm_coverage_", winNames[1], ".txt")

# plot log2-normalised REC8 and other coverage profiles (genome-scale)
  REC8_ChIP <- read.table(REC8_ChIP_file, header = TRUE)
  REC8_input <- read.table(REC8_input_file, header = TRUE)
  REC8_norm <- log2(REC8_ChIP[,2]/REC8_input[,2])
  #REC8_norm <- (REC8_norm-mean(REC8_norm, na.rm = T))/sd(REC8_norm, na.rm = T)

  WT_nuc <- read.table(WT_nuc_file, header = TRUE)
  WT_nakedDNA <- read.table(WT_nakedDNA_file, header = TRUE)
  nuc_norm <- log2(WT_nuc[,2]/WT_nakedDNA[,2])
  #nuc_norm <- (nuc_norm-mean(nuc_norm, na.rm = T))/sd(nuc_norm, na.rm = T)

  H2AZ_ChIP <- read.table(H2AZ_ChIP_file, header = TRUE)
  H2A_input <- read.table(H2A_input_file, header = TRUE)
  H2AZ_norm <- log2(H2AZ_ChIP[,2]/H2A_input[,2])
  #H2AZ_norm <- (H2AZ_norm-mean(H2AZ_norm, na.rm = T))/sd(H2AZ_norm, na.rm = T)

  WT_H3K4me3_ChIP <- read.table(WT_H3K4me3_ChIP_file, header = TRUE)
  WT_H3K4me3_input <- read.table(WT_H3K4me3_input_file, header = TRUE)
  H3K4me3_norm <- log2(WT_H3K4me3_ChIP[,2]/WT_H3K4me3_input[,2])
  #H3K4me3_norm <- (H3K4me3_norm-mean(H3K4me3_norm, na.rm = T))/sd(H3K4me3_norm, na.rm = T)

  cumGeneDat <- read.table(gene_density_file, header = TRUE)
  
  cumMethDat200 <- read.table(DNA_meth_file200, header = TRUE)
  cumMethDat10 <- read.table(DNA_meth_file10, header = TRUE)

  H2AW_ChIP <- read.table(H2AW_ChIP_file, header = TRUE)
  H2AW_norm <- log2(H2AW_ChIP[,2]/H2A_input[,2])
  #H2AW_norm <- (H2AW_norm-mean(H2AW_norm, na.rm = T))/sd(H2AW_norm, na.rm = T)

  WT_H3K9me2_ChIP <- read.table(WT_H3K9me2_ChIP_file, header = TRUE)
  WT_H3K9me2_input <- read.table(WT_H3K9me2_input_file, header = TRUE)
  H3K9me2_norm <- log2(WT_H3K9me2_ChIP[,2]/WT_H3K9me2_input[,2])
  #H3K9me2_norm <- (H3K9me2_norm-mean(H3K9me2_norm, na.rm = T))/sd(H3K9me2_norm, na.rm = T)

  cumTEDat <- read.table(TE_density_file, header = TRUE)

  WT_SPO11_ChIP4 <- read.table(WT_SPO11_ChIP4_file, header = TRUE)
  SPO11_ChIP4_norm <- log2(WT_SPO11_ChIP4[,2]/REC8_input[,2])
  #SPO11_ChIP4_norm <- (SPO11_ChIP4_norm-mean(SPO11_ChIP4_norm, na.rm = T))/sd(SPO11_ChIP4_norm, na.rm = T)

  WT_SPO11_oligos_RPI1 <- read.table(WT_SPO11_oligos_RPI1_file, header = TRUE)
  WT_nakedDNA_T50R1 <- read.table(WT_nakedDNA_T50R1_file, header = TRUE)
  SPO11_oligos_RPI1_norm <- log2(WT_SPO11_oligos_RPI1[,2]/WT_nakedDNA_T50R1[,2])
  #SPO11_oligos_RPI1_norm <- (SPO11_oligos_RPI1_norm-mean(SPO11_oligos_RPI1_norm, na.rm = T))/sd(SPO11_oligos_RPI1_norm, na.rm = T)

  cumCODat <- read.table(CO_density_file, header = TRUE)
  cumCODatRed <- read.table(COReduced_density_file, header = TRUE)

  MSH4_ChIP <- read.table(MSH4_ChIP_file, header = TRUE)
  MSH4_input <- read.table(MSH4_input_file, header = TRUE)
  MSH4_norm <- log2(MSH4_ChIP[,2]/MSH4_input[,2])
  #MSH4_norm <- (MSH4_norm-mean(MSH4_norm, na.rm = T))/sd(MSH4_norm, na.rm = T)

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

filt_nuc_norm <- stats::filter(nuc_norm, ma)
which_na <- which(is.na(filt_nuc_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_nuc_norm[left_na[length(left_na)]+1]
filt_nuc_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_nuc_norm[right_na[1]-1]
filt_nuc_norm[right_na] <- right_val
filt_nuc_norm_noNA <- filt_nuc_norm[!is.na(filt_nuc_norm)]

filt_H2AZ_norm <- stats::filter(H2AZ_norm, ma)
which_na <- which(is.na(filt_H2AZ_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AZ_norm[left_na[length(left_na)]+1]
filt_H2AZ_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AZ_norm[right_na[1]-1]
filt_H2AZ_norm[right_na] <- right_val
filt_H2AZ_norm_noNA <- filt_H2AZ_norm[!is.na(filt_H2AZ_norm)]

filt_H3K4me3_norm <- stats::filter(H3K4me3_norm, ma)
which_na <- which(is.na(filt_H3K4me3_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H3K4me3_norm[left_na[length(left_na)]+1]
filt_H3K4me3_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H3K4me3_norm[right_na[1]-1]
filt_H3K4me3_norm[right_na] <- right_val
filt_H3K4me3_norm_noNA <- filt_H3K4me3_norm[!is.na(filt_H3K4me3_norm)]

filt_cumGeneDat <- stats::filter(cumGeneDat[,2], ma)
which_na <- which(is.na(filt_cumGeneDat) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumGeneDat[left_na[length(left_na)]+1]
filt_cumGeneDat[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumGeneDat[right_na[1]-1]
filt_cumGeneDat[right_na] <- right_val
filt_cumGeneDat_noNA <- filt_cumGeneDat[!is.na(filt_cumGeneDat)]

j = 5
ma <- rep(1, test[j])/test[j]
filt_wtCG200 <- stats::filter(cumMethDat200[,2], ma); filt_wtCG200_n <- filt_wtCG200[!is.na(filt_wtCG200)]
filt_wtCHG200 <- stats::filter(cumMethDat200[,3], ma); filt_wtCHG200_n <- filt_wtCHG200[!is.na(filt_wtCHG200)]
filt_wtCHH200 <- stats::filter(cumMethDat200[,4], ma); filt_wtCHH200_n <- filt_wtCHH200[!is.na(filt_wtCHH200)]
filt_allcntxts200 <- stats::filter(cumMethDat200[,5], ma); filt_allcntxts200_n <- filt_allcntxts200[!is.na(filt_allcntxts200)]
filt_wtMeth200 <- cbind(cumMethDat200[,1], filt_wtCG200, filt_wtCHG200, filt_wtCHH200, filt_allcntxts200)
#write.table(filt_wtMeth200, file = paste0(outDir, "filt_wtMeth_genome_200kb.txt"))
#filt_wtMeth200 <- read.table(file = paste0(outDir, "filt_wtMeth_genome_200kb.txt"))

j = 100
ma <- rep(1, test[j])/test[j]
filt_wtCG10 <- stats::filter(cumMethDat10[,2], ma); filt_wtCG10_n <- filt_wtCG10[!is.na(filt_wtCG10)]
filt_wtCHG10 <- stats::filter(cumMethDat10[,3], ma); filt_wtCHG10_n <- filt_wtCHG10[!is.na(filt_wtCHG10)]
filt_wtCHH10 <- stats::filter(cumMethDat10[,4], ma); filt_wtCHH10_n <- filt_wtCHH10[!is.na(filt_wtCHH10)]
filt_allcntxts10 <- stats::filter(cumMethDat10[,5], ma); filt_allcntxts10_n <- filt_allcntxts10[!is.na(filt_allcntxts10)]
filt_wtMeth10 <- cbind(cumMethDat10[,1], filt_wtCG10, filt_wtCHG10, filt_wtCHH10, filt_allcntxts10)
#write.table(filt_wtMeth10, file = paste0(outDir, "filt_wtMeth_genome_10kb.txt"))
#filt_wtMeth10 <- read.table(file = paste0(outDir, "filt_wtMeth_genome_10kb.txt"))

filt_H2AW_norm <- stats::filter(H2AW_norm, ma)
which_na <- which(is.na(filt_H2AW_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AW_norm[left_na[length(left_na)]+1]
filt_H2AW_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AW_norm[right_na[1]-1]
filt_H2AW_norm[right_na] <- right_val
filt_H2AW_norm_noNA <- filt_H2AW_norm[!is.na(filt_H2AW_norm)]

filt_H3K9me2_norm <- stats::filter(H3K9me2_norm, ma)
which_na <- which(is.na(filt_H3K9me2_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H3K9me2_norm[left_na[length(left_na)]+1]
filt_H3K9me2_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H3K9me2_norm[right_na[1]-1]
filt_H3K9me2_norm[right_na] <- right_val
filt_H3K9me2_norm_noNA <- filt_H3K9me2_norm[!is.na(filt_H3K9me2_norm)]

filt_cumTEDat <- stats::filter(cumTEDat[,2], ma)
which_na <- which(is.na(filt_cumTEDat) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumTEDat[left_na[length(left_na)]+1]
filt_cumTEDat[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumTEDat[right_na[1]-1]
filt_cumTEDat[right_na] <- right_val
filt_cumTEDat_noNA <- filt_cumTEDat[!is.na(filt_cumTEDat)]

filt_SPO11_ChIP4_norm <- stats::filter(SPO11_ChIP4_norm, ma)
which_na <- which(is.na(filt_SPO11_ChIP4_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_SPO11_ChIP4_norm[left_na[length(left_na)]+1]
filt_SPO11_ChIP4_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_SPO11_ChIP4_norm[right_na[1]-1]
filt_SPO11_ChIP4_norm[right_na] <- right_val
filt_SPO11_ChIP4_norm_noNA <- filt_SPO11_ChIP4_norm[!is.na(filt_SPO11_ChIP4_norm)]

filt_SPO11_oligos_RPI1_norm <- stats::filter(SPO11_oligos_RPI1_norm, ma)
which_na <- which(is.na(filt_SPO11_oligos_RPI1_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_SPO11_oligos_RPI1_norm[left_na[length(left_na)]+1]
filt_SPO11_oligos_RPI1_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_SPO11_oligos_RPI1_norm[right_na[1]-1]
filt_SPO11_oligos_RPI1_norm[right_na] <- right_val
filt_SPO11_oligos_RPI1_norm_noNA <- filt_SPO11_oligos_RPI1_norm[!is.na(filt_SPO11_oligos_RPI1_norm)]

filt_MSH4_norm <- stats::filter(MSH4_norm, ma)
which_na <- which(is.na(filt_MSH4_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_MSH4_norm[left_na[length(left_na)]+1]
filt_MSH4_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_MSH4_norm[right_na[1]-1]
filt_MSH4_norm[right_na] <- right_val
filt_MSH4_norm_noNA <- filt_MSH4_norm[!is.na(filt_MSH4_norm)]

j = 200
ma <- rep(1, test[j])/test[j]
filt_cumCODat <- stats::filter(cumCODat[,2], ma)
which_na <- which(is.na(filt_cumCODat) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumCODat[left_na[length(left_na)]+1]
filt_cumCODat[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumCODat[right_na[1]-1]
filt_cumCODat[right_na] <- right_val
filt_cumCODat_noNA <- filt_cumCODat[!is.na(filt_cumCODat)]

filt_cumCODatRed <- stats::filter(cumCODatRed[,2], ma)
which_na <- which(is.na(filt_cumCODatRed) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumCODatRed[left_na[length(left_na)]+1]
filt_cumCODatRed[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumCODatRed[right_na[1]-1]
filt_cumCODatRed[right_na] <- right_val
filt_cumCODatRed_noNA <- filt_cumCODatRed[!is.na(filt_cumCODatRed)]


REC8vsOthersGenomePlot <- function(xplot, otherFiltNorm, otherFiltNormNoNA, otherYlabel) {
  plot(xplot, otherFiltNorm, type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(otherFiltNormNoNA), max(otherFiltNormNoNA)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(otherFiltNorm), lwd.tick = 1.5)
  #p <- par('usr')
  #text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.5), labels = otherYlabel, xpd = NA, srt = -90, col = "blue")
  mtext(side = 4, line = 2.5, cex = 1, text = otherYlabel, col = "blue")
  par(new = T)
  plot(xplot, filt_REC8_norm, type = "l", lwd = 1.5, col = "red",
       ylim = c(min(filt_REC8_norm_noNA), max(filt_REC8_norm_noNA)),
       xlab = "",
       ylab = "")
  axis(side = 2, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = "REC8", col = "red")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
}

REC8vsDNAmethGenomePlot <- function(xplot, xplotMeth, otherFiltNorm, otherFiltNormNoNA, otherYlabel) {
  plot(xplotMeth, filt_wtMeth200[,2], type = "l", lwd = 1.5, col = "navy",
       ylim = c(min(filt_wtCG200_n, filt_wtCHG200_n, filt_wtCHH200_n), max(filt_wtCG200_n, filt_wtCHG200_n, filt_wtCHH200_n)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(filt_wtMeth200[,2]), lwd.tick = 1.5)
  mtext(side = 4, line = 2.5, cex = 1, text = otherYlabel, col = "blue")
  par(new = T)
  plot(xplotMeth, filt_wtMeth200[,3], type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(filt_wtCG200_n, filt_wtCHG200_n, filt_wtCHH200_n), max(filt_wtCG200_n, filt_wtCHG200_n, filt_wtCHH200_n)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n")
  par(new = T)
  plot(xplotMeth, filt_wtMeth200[,4], type = "l", lwd = 1.5, col = "steelblue1",
       ylim = c(min(filt_wtCG200_n, filt_wtCHG200_n, filt_wtCHH200_n), max(filt_wtCG200_n, filt_wtCHG200_n, filt_wtCHH200_n)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n") 
  par(new = T)
  plot(xplot, filt_REC8_norm, type = "l", lwd = 1.5, col = "red",
       ylim = c(min(filt_REC8_norm_noNA), max(filt_REC8_norm_noNA)),
       xlab = "",
       ylab = "")
  axis(side = 2, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = "REC8", col = "red")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
  legend("topleft",
         legend = c("CG", "CHG", "CHH"),
         col = c("navy", "blue", "steelblue1"),
         text.col = c("navy", "blue", "steelblue1"),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}


pdf(file = paste0(outDir, "REC8_vs_Others_log2norm_genomeplot_", winNames[1], "_h10w15.pdf"), height = 10, width = 15)
par(mfcol = c(4, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_nuc_norm, otherFiltNormNoNA = filt_nuc_norm_noNA, otherYlabel = "MNase")
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_H2AZ_norm, otherFiltNormNoNA = filt_H2AZ_norm_noNA, otherYlabel = "H2A.Z")
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_H3K4me3_norm, otherFiltNormNoNA = filt_H3K4me3_norm_noNA, otherYlabel = "H3K4me3")
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_cumGeneDat, otherFiltNormNoNA = filt_cumGeneDat_noNA, otherYlabel = "Genes")
REC8vsDNAmethGenomePlot(xplot = REC8_ChIP[,1], xplotMeth = filt_wtMeth200[,1], otherYlabel = "DNA methylation proportion") 
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_H2AW_norm, otherFiltNormNoNA = filt_H2AW_norm_noNA, otherYlabel = "H2A.W")
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_H3K9me2_norm, otherFiltNormNoNA = filt_H3K9me2_norm_noNA, otherYlabel = "H3K9me2")
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_cumTEDat, otherFiltNormNoNA = filt_cumTEDat_noNA, otherYlabel = "TEs")

dev.off()

pdf(file = paste0(outDir, "REC8_vs_SPO11-1_COs_MSH4_log2_genomeplot_", winNames[1], "_h10w15.pdf"), height = 10, width = 15)
par(mfcol = c(4, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_SPO11_oligos_RPI1_norm, otherFiltNormNoNA = filt_SPO11_oligos_RPI1_norm_noNA, otherYlabel = "SPO11-1-oligos")
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_cumCODatRed, otherFiltNormNoNA = filt_cumCODatRed_noNA, otherYlabel = "Crossovers")
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_SPO11_ChIP4_norm, otherFiltNormNoNA = filt_SPO11_ChIP4_norm_noNA, otherYlabel = "SPO11-1 ChIP")
REC8vsOthersGenomePlot(xplot = REC8_ChIP[,1], otherFiltNorm = filt_MSH4_norm, otherFiltNormNoNA = filt_MSH4_norm_noNA, otherYlabel = "MSH4")

dev.off()

# Create genome-wide Spearman's rho correlation matrix
allDF <- data.frame(filt_REC8_norm, filt_nuc_norm, filt_H2AZ_norm, filt_H3K4me3_norm, filt_cumGeneDat, filt_H2AW_norm, filt_H3K9me2_norm, filt_cumTEDat, filt_wtCG10, filt_wtCHG10, filt_wtCHH10, filt_allcntxts10)
colnames(allDF) <- c("REC8", "MNase", "H2A.Z", "H3K4me3", "Genes", "H2A.W", "H3K9me2", "TEs", "CG", "CHG", "CHH", "DNA meth")
allDF_corMat <- cor(allDF, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
#pdf(file = paste0(outDir, "genome-wide_correlation_matrix_", winNames[1], "_circle_noZscore.pdf"))
#corrplot(allDF_corMat, method = "circle", type = "upper", col = col1(20), tl.col = "black",
#         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0),
#         title = paste0("Genome-wide Spearman correlation matrix (10-kb windows)"))
#dev.off()
#
#pdf(file = paste0(outDir, "genome-wide_correlation_matrix_", winNames[1], "_square_noZscore.pdf"))
#corrplot(allDF_corMat, method = "square", type = "upper", col = col1(20), tl.col = "black",
#         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0),
#         title = paste0("Genome-wide Spearman correlation matrix (10-kb windows)"))
#dev.off()

pdf(file = paste0(outDir, "genome-wide_correlation_matrix_", winNames[1], "_colouronly_noZscore.pdf"))
corrplot(allDF_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0), tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.8,
         title = paste0("Genome-wide Spearman correlation matrix (10-kb windows)"))
dev.off()

#pdf(file = paste0(outDir, "genome-wide_correlation_matrix_", winNames[1], "_number_noZscore.pdf"))
#corrplot(allDF_corMat, method = "number", type = "upper", col = col1(20), tl.col = "black",
#         addgrid.col = "white", mar = c(0,0,1,0),
#         title = paste0("Genome-wide Spearman correlation matrix (10-kb windows)"))
#dev.off()
#
#pdf(file = paste0(outDir, "genome-wide_correlation_matrix_", winNames[1], "_pie_noZscore.pdf"))
#corrplot(allDF_corMat, method = "pie", type = "upper", col = col1(20), tl.col = "black",
#         addgrid.col = "white", addCoef.col = "black", mar = c(0,0,1,0),
#         title = paste0("Genome-wide Spearman correlation matrix (10-kb windows)"))
#dev.off()

# Create genome-wide Spearman's rho correlation matrix including SPO11-1-oligos, SPO11-1 ChIP, crossovers, and MSH4
allDF <- data.frame(filt_REC8_norm, filt_nuc_norm, filt_SPO11_ChIP4_norm, filt_MSH4_norm, filt_H2AW_norm, filt_allcntxts10, filt_H3K9me2_norm, filt_cumGeneDat, filt_H2AZ_norm, filt_H3K4me3_norm, filt_SPO11_oligos_RPI1_norm, filt_cumCODatRed)
colnames(allDF) <- c("REC8", "MNase", "SPO11-1 ChIP", "MSH4", "H2A.W", "DNA meth", "H3K9me2", "Genes", "H2A.Z", "H3K4me3", "SPO11-1-oligos", "Crossovers")
allDF_corMat <- cor(allDF, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "genome-wide_correlation_matrix_", winNames[1], "_colouronly_v02_noZscore.pdf"))
corrplot(allDF_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0), tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.8,
         title = paste0("Genome-wide Spearman correlation matrix (10-kb windows)"))
dev.off()


# Separate arms and pericentromeric regions to create separate Spearman's rho correlation matrices
REC8_ChIP_armL <- NULL
REC8_ChIP_armR <- NULL
REC8_ChIP_arm <- NULL
REC8_ChIP_peri <- NULL
REC8_input_armL <- NULL
REC8_input_armR <- NULL
REC8_input_arm <- NULL
REC8_input_peri <- NULL

WT_nuc_armL <- NULL
WT_nuc_armR <- NULL
WT_nuc_arm <- NULL
WT_nuc_peri <- NULL
WT_nakedDNA_armL <- NULL
WT_nakedDNA_armR <- NULL
WT_nakedDNA_arm <- NULL
WT_nakedDNA_peri <- NULL

H2AZ_ChIP_armL <- NULL
H2AZ_ChIP_armR <- NULL
H2AZ_ChIP_arm <- NULL
H2AZ_ChIP_peri <- NULL
H2A_input_armL <- NULL
H2A_input_armR <- NULL
H2A_input_arm <- NULL
H2A_input_peri <- NULL

WT_H3K4me3_ChIP_armL <- NULL
WT_H3K4me3_ChIP_armR <- NULL
WT_H3K4me3_ChIP_arm <- NULL
WT_H3K4me3_ChIP_peri <- NULL
WT_H3K4me3_input_armL <- NULL
WT_H3K4me3_input_armR <- NULL
WT_H3K4me3_input_arm <- NULL
WT_H3K4me3_input_peri <- NULL

cumGeneDat_armL <- NULL
cumGeneDat_armR <- NULL
cumGeneDat_arm <- NULL
cumGeneDat_peri <- NULL

cumMethDat10_armL <- NULL
cumMethDat10_armR <- NULL
cumMethDat10_arm <- NULL
cumMethDat10_peri <- NULL

H2AW_ChIP_armL <- NULL
H2AW_ChIP_armR <- NULL
H2AW_ChIP_arm <- NULL
H2AW_ChIP_peri <- NULL

WT_H3K9me2_ChIP_armL <- NULL
WT_H3K9me2_ChIP_armR <- NULL
WT_H3K9me2_ChIP_arm <- NULL
WT_H3K9me2_ChIP_peri <- NULL
WT_H3K9me2_input_armL <- NULL
WT_H3K9me2_input_armR <- NULL
WT_H3K9me2_input_arm <- NULL
WT_H3K9me2_input_peri <- NULL

cumTEDat_armL <- NULL
cumTEDat_armR <- NULL
cumTEDat_arm <- NULL
cumTEDat_peri <- NULL

WT_SPO11_ChIP4_armL <- NULL
WT_SPO11_ChIP4_armR <- NULL
WT_SPO11_ChIP4_arm <- NULL
WT_SPO11_ChIP4_peri <- NULL

WT_SPO11_oligos_RPI1_armL <- NULL
WT_SPO11_oligos_RPI1_armR <- NULL
WT_SPO11_oligos_RPI1_arm <- NULL
WT_SPO11_oligos_RPI1_peri <- NULL
WT_nakedDNA_T50R1_armL <- NULL
WT_nakedDNA_T50R1_armR <- NULL
WT_nakedDNA_T50R1_arm <- NULL
WT_nakedDNA_T50R1_peri <- NULL

cumCODat_armL <- NULL
cumCODat_armR <- NULL
cumCODat_arm <- NULL
cumCODat_peri <- NULL

cumCODatRed_armL <- NULL
cumCODatRed_armR <- NULL
cumCODatRed_arm <- NULL
cumCODatRed_peri <- NULL

MSH4_ChIP_armL <- NULL
MSH4_ChIP_armR <- NULL
MSH4_ChIP_arm <- NULL
MSH4_ChIP_peri <- NULL
MSH4_input_armL <- NULL
MSH4_input_armR <- NULL
MSH4_input_arm <- NULL
MSH4_input_peri <- NULL
for(i in 1:5) {
  print(i)
  REC8_ChIP_armL_tmp <- REC8_ChIP[REC8_ChIP[,1] > sumchr[i] & REC8_ChIP[,1] < pericenStart[i] & REC8_ChIP[,1] < sumchr[i+1],]
  REC8_ChIP_armR_tmp <- REC8_ChIP[REC8_ChIP[,1] > pericenEnd[i] & REC8_ChIP[,1] <= sumchr[i+1],]
  REC8_ChIP_arm_tmp <- REC8_ChIP[(REC8_ChIP[,1] > sumchr[i] & REC8_ChIP[,1] < pericenStart[i] & REC8_ChIP[,1] < sumchr[i+1]) |
                                 (REC8_ChIP[,1] > pericenEnd[i] & REC8_ChIP[,1] <= sumchr[i+1]),]
  REC8_ChIP_peri_tmp <- REC8_ChIP[REC8_ChIP[,1] >= pericenStart[i] & REC8_ChIP[,1] <= pericenEnd[i],]
  REC8_ChIP_peri <- rbind(REC8_ChIP_peri, REC8_ChIP_peri_tmp)
  REC8_ChIP_armL <- rbind(REC8_ChIP_armL, REC8_ChIP_armL_tmp)
  REC8_ChIP_armR <- rbind(REC8_ChIP_armR, REC8_ChIP_armR_tmp)
  REC8_ChIP_arm <- rbind(REC8_ChIP_arm, REC8_ChIP_arm_tmp)
  print(i)
  REC8_input_armL_tmp <- REC8_input[REC8_input[,1] > sumchr[i] & REC8_input[,1] < pericenStart[i] & REC8_input[,1] < sumchr[i+1],]
  REC8_input_armR_tmp <- REC8_input[REC8_input[,1] > pericenEnd[i] & REC8_input[,1] <= sumchr[i+1],]
  REC8_input_arm_tmp <- REC8_input[(REC8_input[,1] > sumchr[i] & REC8_input[,1] < pericenStart[i] & REC8_input[,1] < sumchr[i+1]) |
                                 (REC8_input[,1] > pericenEnd[i] & REC8_input[,1] <= sumchr[i+1]),]
  REC8_input_peri_tmp <- REC8_input[REC8_input[,1] >= pericenStart[i] & REC8_input[,1] <= pericenEnd[i],]
  REC8_input_peri <- rbind(REC8_input_peri, REC8_input_peri_tmp)
  REC8_input_armL <- rbind(REC8_input_armL, REC8_input_armL_tmp)
  REC8_input_armR <- rbind(REC8_input_armR, REC8_input_armR_tmp)
  REC8_input_arm <- rbind(REC8_input_arm, REC8_input_arm_tmp)

  print(i)
  WT_nuc_armL_tmp <- WT_nuc[WT_nuc[,1] > sumchr[i] & WT_nuc[,1] < pericenStart[i] & WT_nuc[,1] < sumchr[i+1],]
  WT_nuc_armR_tmp <- WT_nuc[WT_nuc[,1] > pericenEnd[i] & WT_nuc[,1] <= sumchr[i+1],]
  WT_nuc_arm_tmp <- WT_nuc[(WT_nuc[,1] > sumchr[i] & WT_nuc[,1] < pericenStart[i] & WT_nuc[,1] < sumchr[i+1]) |
                                 (WT_nuc[,1] > pericenEnd[i] & WT_nuc[,1] <= sumchr[i+1]),]
  WT_nuc_peri_tmp <- WT_nuc[WT_nuc[,1] >= pericenStart[i] & WT_nuc[,1] <= pericenEnd[i],]
  WT_nuc_peri <- rbind(WT_nuc_peri, WT_nuc_peri_tmp)
  WT_nuc_armL <- rbind(WT_nuc_armL, WT_nuc_armL_tmp)
  WT_nuc_armR <- rbind(WT_nuc_armR, WT_nuc_armR_tmp)
  WT_nuc_arm <- rbind(WT_nuc_arm, WT_nuc_arm_tmp)
  print(i)
  WT_nakedDNA_armL_tmp <- WT_nakedDNA[WT_nakedDNA[,1] > sumchr[i] & WT_nakedDNA[,1] < pericenStart[i] & WT_nakedDNA[,1] < sumchr[i+1],]
  WT_nakedDNA_armR_tmp <- WT_nakedDNA[WT_nakedDNA[,1] > pericenEnd[i] & WT_nakedDNA[,1] <= sumchr[i+1],]
  WT_nakedDNA_arm_tmp <- WT_nakedDNA[(WT_nakedDNA[,1] > sumchr[i] & WT_nakedDNA[,1] < pericenStart[i] & WT_nakedDNA[,1] < sumchr[i+1]) | 
                                 (WT_nakedDNA[,1] > pericenEnd[i] & WT_nakedDNA[,1] <= sumchr[i+1]),]
  WT_nakedDNA_peri_tmp <- WT_nakedDNA[WT_nakedDNA[,1] >= pericenStart[i] & WT_nakedDNA[,1] <= pericenEnd[i],]
  WT_nakedDNA_peri <- rbind(WT_nakedDNA_peri, WT_nakedDNA_peri_tmp)
  WT_nakedDNA_armL <- rbind(WT_nakedDNA_armL, WT_nakedDNA_armL_tmp)
  WT_nakedDNA_armR <- rbind(WT_nakedDNA_armR, WT_nakedDNA_armR_tmp)
  WT_nakedDNA_arm <- rbind(WT_nakedDNA_arm, WT_nakedDNA_arm_tmp)

  print(i)
  H2AZ_ChIP_armL_tmp <- H2AZ_ChIP[H2AZ_ChIP[,1] > sumchr[i] & H2AZ_ChIP[,1] < pericenStart[i] & H2AZ_ChIP[,1] < sumchr[i+1],]
  H2AZ_ChIP_armR_tmp <- H2AZ_ChIP[H2AZ_ChIP[,1] > pericenEnd[i] & H2AZ_ChIP[,1] <= sumchr[i+1],]
  H2AZ_ChIP_arm_tmp <- H2AZ_ChIP[(H2AZ_ChIP[,1] > sumchr[i] & H2AZ_ChIP[,1] < pericenStart[i] & H2AZ_ChIP[,1] < sumchr[i+1]) |
                                 (H2AZ_ChIP[,1] > pericenEnd[i] & H2AZ_ChIP[,1] <= sumchr[i+1]),]
  H2AZ_ChIP_peri_tmp <- H2AZ_ChIP[H2AZ_ChIP[,1] >= pericenStart[i] & H2AZ_ChIP[,1] <= pericenEnd[i],]
  H2AZ_ChIP_peri <- rbind(H2AZ_ChIP_peri, H2AZ_ChIP_peri_tmp)
  H2AZ_ChIP_armL <- rbind(H2AZ_ChIP_armL, H2AZ_ChIP_armL_tmp)
  H2AZ_ChIP_armR <- rbind(H2AZ_ChIP_armR, H2AZ_ChIP_armR_tmp)
  H2AZ_ChIP_arm <- rbind(H2AZ_ChIP_arm, H2AZ_ChIP_arm_tmp)
  print(i)
  H2A_input_armL_tmp <- H2A_input[H2A_input[,1] > sumchr[i] & H2A_input[,1] < pericenStart[i] & H2A_input[,1] < sumchr[i+1],]
  H2A_input_armR_tmp <- H2A_input[H2A_input[,1] > pericenEnd[i] & H2A_input[,1] <= sumchr[i+1],]
  H2A_input_arm_tmp <- H2A_input[(H2A_input[,1] > sumchr[i] & H2A_input[,1] < pericenStart[i] & H2A_input[,1] < sumchr[i+1]) | 
                                 (H2A_input[,1] > pericenEnd[i] & H2A_input[,1] <= sumchr[i+1]),]
  H2A_input_peri_tmp <- H2A_input[H2A_input[,1] >= pericenStart[i] & H2A_input[,1] <= pericenEnd[i],]
  H2A_input_peri <- rbind(H2A_input_peri, H2A_input_peri_tmp)
  H2A_input_armL <- rbind(H2A_input_armL, H2A_input_armL_tmp)
  H2A_input_armR <- rbind(H2A_input_armR, H2A_input_armR_tmp)
  H2A_input_arm <- rbind(H2A_input_arm, H2A_input_arm_tmp)

  print(i)
  WT_H3K4me3_ChIP_armL_tmp <- WT_H3K4me3_ChIP[WT_H3K4me3_ChIP[,1] > sumchr[i] & WT_H3K4me3_ChIP[,1] < pericenStart[i] & WT_H3K4me3_ChIP[,1] < sumchr[i+1],]
  WT_H3K4me3_ChIP_armR_tmp <- WT_H3K4me3_ChIP[WT_H3K4me3_ChIP[,1] > pericenEnd[i] & WT_H3K4me3_ChIP[,1] <= sumchr[i+1],]
  WT_H3K4me3_ChIP_arm_tmp <- WT_H3K4me3_ChIP[(WT_H3K4me3_ChIP[,1] > sumchr[i] & WT_H3K4me3_ChIP[,1] < pericenStart[i] & WT_H3K4me3_ChIP[,1] < sumchr[i+1]) |
                                 (WT_H3K4me3_ChIP[,1] > pericenEnd[i] & WT_H3K4me3_ChIP[,1] <= sumchr[i+1]),]
  WT_H3K4me3_ChIP_peri_tmp <- WT_H3K4me3_ChIP[WT_H3K4me3_ChIP[,1] >= pericenStart[i] & WT_H3K4me3_ChIP[,1] <= pericenEnd[i],]
  WT_H3K4me3_ChIP_peri <- rbind(WT_H3K4me3_ChIP_peri, WT_H3K4me3_ChIP_peri_tmp)
  WT_H3K4me3_ChIP_armL <- rbind(WT_H3K4me3_ChIP_armL, WT_H3K4me3_ChIP_armL_tmp)
  WT_H3K4me3_ChIP_armR <- rbind(WT_H3K4me3_ChIP_armR, WT_H3K4me3_ChIP_armR_tmp)
  WT_H3K4me3_ChIP_arm <- rbind(WT_H3K4me3_ChIP_arm, WT_H3K4me3_ChIP_arm_tmp)
  print(i)
  WT_H3K4me3_input_armL_tmp <- WT_H3K4me3_input[WT_H3K4me3_input[,1] > sumchr[i] & WT_H3K4me3_input[,1] < pericenStart[i] & WT_H3K4me3_input[,1] < sumchr[i+1],]
  WT_H3K4me3_input_armR_tmp <- WT_H3K4me3_input[WT_H3K4me3_input[,1] > pericenEnd[i] & WT_H3K4me3_input[,1] <= sumchr[i+1],]
  WT_H3K4me3_input_arm_tmp <- WT_H3K4me3_input[(WT_H3K4me3_input[,1] > sumchr[i] & WT_H3K4me3_input[,1] < pericenStart[i] & WT_H3K4me3_input[,1] < sumchr[i+1]) | 
                                 (WT_H3K4me3_input[,1] > pericenEnd[i] & WT_H3K4me3_input[,1] <= sumchr[i+1]),]
  WT_H3K4me3_input_peri_tmp <- WT_H3K4me3_input[WT_H3K4me3_input[,1] >= pericenStart[i] & WT_H3K4me3_input[,1] <= pericenEnd[i],]
  WT_H3K4me3_input_peri <- rbind(WT_H3K4me3_input_peri, WT_H3K4me3_input_peri_tmp)
  WT_H3K4me3_input_armL <- rbind(WT_H3K4me3_input_armL, WT_H3K4me3_input_armL_tmp)
  WT_H3K4me3_input_armR <- rbind(WT_H3K4me3_input_armR, WT_H3K4me3_input_armR_tmp)
  WT_H3K4me3_input_arm <- rbind(WT_H3K4me3_input_arm, WT_H3K4me3_input_arm_tmp)

  print(i)
  cumGeneDat_armL_tmp <- cumGeneDat[cumGeneDat[,1] > sumchr[i] & cumGeneDat[,1] < pericenStart[i] & cumGeneDat[,1] < sumchr[i+1],]
  cumGeneDat_armR_tmp <- cumGeneDat[cumGeneDat[,1] > pericenEnd[i] & cumGeneDat[,1] <= sumchr[i+1],]
  cumGeneDat_arm_tmp <- cumGeneDat[(cumGeneDat[,1] > sumchr[i] & cumGeneDat[,1] < pericenStart[i] & cumGeneDat[,1] < sumchr[i+1]) |
                                 (cumGeneDat[,1] > pericenEnd[i] & cumGeneDat[,1] <= sumchr[i+1]),]
  cumGeneDat_peri_tmp <- cumGeneDat[cumGeneDat[,1] >= pericenStart[i] & cumGeneDat[,1] <= pericenEnd[i],]
  cumGeneDat_peri <- rbind(cumGeneDat_peri, cumGeneDat_peri_tmp)
  cumGeneDat_armL <- rbind(cumGeneDat_armL, cumGeneDat_armL_tmp)
  cumGeneDat_armR <- rbind(cumGeneDat_armR, cumGeneDat_armR_tmp)
  cumGeneDat_arm <- rbind(cumGeneDat_arm, cumGeneDat_arm_tmp)

  print(i)
  cumMethDat10_armL_tmp <- cumMethDat10[cumMethDat10[,1] > sumchr[i] & cumMethDat10[,1] < pericenStart[i] & cumMethDat10[,1] < sumchr[i+1],]
  cumMethDat10_armR_tmp <- cumMethDat10[cumMethDat10[,1] > pericenEnd[i] & cumMethDat10[,1] <= sumchr[i+1],]
  cumMethDat10_arm_tmp <- cumMethDat10[(cumMethDat10[,1] > sumchr[i] & cumMethDat10[,1] < pericenStart[i] & cumMethDat10[,1] < sumchr[i+1]) | 
                                 (cumMethDat10[,1] > pericenEnd[i] & cumMethDat10[,1] <= sumchr[i+1]),]
  cumMethDat10_peri_tmp <- cumMethDat10[cumMethDat10[,1] >= pericenStart[i] & cumMethDat10[,1] <= pericenEnd[i],]
  cumMethDat10_peri <- rbind(cumMethDat10_peri, cumMethDat10_peri_tmp)
  cumMethDat10_armL <- rbind(cumMethDat10_armL, cumMethDat10_armL_tmp)
  cumMethDat10_armR <- rbind(cumMethDat10_armR, cumMethDat10_armR_tmp)
  cumMethDat10_arm <- rbind(cumMethDat10_arm, cumMethDat10_arm_tmp)

  print(i)
  H2AW_ChIP_armL_tmp <- H2AW_ChIP[H2AW_ChIP[,1] > sumchr[i] & H2AW_ChIP[,1] < pericenStart[i] & H2AW_ChIP[,1] < sumchr[i+1],]
  H2AW_ChIP_armR_tmp <- H2AW_ChIP[H2AW_ChIP[,1] > pericenEnd[i] & H2AW_ChIP[,1] <= sumchr[i+1],]
  H2AW_ChIP_arm_tmp <- H2AW_ChIP[(H2AW_ChIP[,1] > sumchr[i] & H2AW_ChIP[,1] < pericenStart[i] & H2AW_ChIP[,1] < sumchr[i+1]) |
                                 (H2AW_ChIP[,1] > pericenEnd[i] & H2AW_ChIP[,1] <= sumchr[i+1]),]
  H2AW_ChIP_peri_tmp <- H2AW_ChIP[H2AW_ChIP[,1] >= pericenStart[i] & H2AW_ChIP[,1] <= pericenEnd[i],]
  H2AW_ChIP_peri <- rbind(H2AW_ChIP_peri, H2AW_ChIP_peri_tmp)
  H2AW_ChIP_armL <- rbind(H2AW_ChIP_armL, H2AW_ChIP_armL_tmp)
  H2AW_ChIP_armR <- rbind(H2AW_ChIP_armR, H2AW_ChIP_armR_tmp)
  H2AW_ChIP_arm <- rbind(H2AW_ChIP_arm, H2AW_ChIP_arm_tmp)

  print(i)
  WT_H3K9me2_ChIP_armL_tmp <- WT_H3K9me2_ChIP[WT_H3K9me2_ChIP[,1] > sumchr[i] & WT_H3K9me2_ChIP[,1] < pericenStart[i] & WT_H3K9me2_ChIP[,1] < sumchr[i+1],]
  WT_H3K9me2_ChIP_armR_tmp <- WT_H3K9me2_ChIP[WT_H3K9me2_ChIP[,1] > pericenEnd[i] & WT_H3K9me2_ChIP[,1] <= sumchr[i+1],]
  WT_H3K9me2_ChIP_arm_tmp <- WT_H3K9me2_ChIP[(WT_H3K9me2_ChIP[,1] > sumchr[i] & WT_H3K9me2_ChIP[,1] < pericenStart[i] & WT_H3K9me2_ChIP[,1] < sumchr[i+1]) |
                                 (WT_H3K9me2_ChIP[,1] > pericenEnd[i] & WT_H3K9me2_ChIP[,1] <= sumchr[i+1]),]
  WT_H3K9me2_ChIP_peri_tmp <- WT_H3K9me2_ChIP[WT_H3K9me2_ChIP[,1] >= pericenStart[i] & WT_H3K9me2_ChIP[,1] <= pericenEnd[i],]
  WT_H3K9me2_ChIP_peri <- rbind(WT_H3K9me2_ChIP_peri, WT_H3K9me2_ChIP_peri_tmp)
  WT_H3K9me2_ChIP_armL <- rbind(WT_H3K9me2_ChIP_armL, WT_H3K9me2_ChIP_armL_tmp)
  WT_H3K9me2_ChIP_armR <- rbind(WT_H3K9me2_ChIP_armR, WT_H3K9me2_ChIP_armR_tmp)
  WT_H3K9me2_ChIP_arm <- rbind(WT_H3K9me2_ChIP_arm, WT_H3K9me2_ChIP_arm_tmp)
  print(i)
  WT_H3K9me2_input_armL_tmp <- WT_H3K9me2_input[WT_H3K9me2_input[,1] > sumchr[i] & WT_H3K9me2_input[,1] < pericenStart[i] & WT_H3K9me2_input[,1] < sumchr[i+1],]
  WT_H3K9me2_input_armR_tmp <- WT_H3K9me2_input[WT_H3K9me2_input[,1] > pericenEnd[i] & WT_H3K9me2_input[,1] <= sumchr[i+1],]
  WT_H3K9me2_input_arm_tmp <- WT_H3K9me2_input[(WT_H3K9me2_input[,1] > sumchr[i] & WT_H3K9me2_input[,1] < pericenStart[i] & WT_H3K9me2_input[,1] < sumchr[i+1]) |
                                 (WT_H3K9me2_input[,1] > pericenEnd[i] & WT_H3K9me2_input[,1] <= sumchr[i+1]),]
  WT_H3K9me2_input_peri_tmp <- WT_H3K9me2_input[WT_H3K9me2_input[,1] >= pericenStart[i] & WT_H3K9me2_input[,1] <= pericenEnd[i],]
  WT_H3K9me2_input_peri <- rbind(WT_H3K9me2_input_peri, WT_H3K9me2_input_peri_tmp)
  WT_H3K9me2_input_armL <- rbind(WT_H3K9me2_input_armL, WT_H3K9me2_input_armL_tmp)
  WT_H3K9me2_input_armR <- rbind(WT_H3K9me2_input_armR, WT_H3K9me2_input_armR_tmp)
  WT_H3K9me2_input_arm <- rbind(WT_H3K9me2_input_arm, WT_H3K9me2_input_arm_tmp)

  print(i)
  cumTEDat_armL_tmp <- cumTEDat[cumTEDat[,1] > sumchr[i] & cumTEDat[,1] < pericenStart[i] & cumTEDat[,1] < sumchr[i+1],]
  cumTEDat_armR_tmp <- cumTEDat[cumTEDat[,1] > pericenEnd[i] & cumTEDat[,1] <= sumchr[i+1],]
  cumTEDat_arm_tmp <- cumTEDat[(cumTEDat[,1] > sumchr[i] & cumTEDat[,1] < pericenStart[i] & cumTEDat[,1] < sumchr[i+1]) |
                                 (cumTEDat[,1] > pericenEnd[i] & cumTEDat[,1] <= sumchr[i+1]),]
  cumTEDat_peri_tmp <- cumTEDat[cumTEDat[,1] >= pericenStart[i] & cumTEDat[,1] <= pericenEnd[i],]
  cumTEDat_peri <- rbind(cumTEDat_peri, cumTEDat_peri_tmp)
  cumTEDat_armL <- rbind(cumTEDat_armL, cumTEDat_armL_tmp)
  cumTEDat_armR <- rbind(cumTEDat_armR, cumTEDat_armR_tmp)
  cumTEDat_arm <- rbind(cumTEDat_arm, cumTEDat_arm_tmp)

  print(i)
  WT_SPO11_ChIP4_armL_tmp <- WT_SPO11_ChIP4[WT_SPO11_ChIP4[,1] > sumchr[i] & WT_SPO11_ChIP4[,1] < pericenStart[i] & WT_SPO11_ChIP4[,1] < sumchr[i+1],]
  WT_SPO11_ChIP4_armR_tmp <- WT_SPO11_ChIP4[WT_SPO11_ChIP4[,1] > pericenEnd[i] & WT_SPO11_ChIP4[,1] <= sumchr[i+1],]
  WT_SPO11_ChIP4_arm_tmp <- WT_SPO11_ChIP4[(WT_SPO11_ChIP4[,1] > sumchr[i] & WT_SPO11_ChIP4[,1] < pericenStart[i] & WT_SPO11_ChIP4[,1] < sumchr[i+1]) |
                                 (WT_SPO11_ChIP4[,1] > pericenEnd[i] & WT_SPO11_ChIP4[,1] <= sumchr[i+1]),]
  WT_SPO11_ChIP4_peri_tmp <- WT_SPO11_ChIP4[WT_SPO11_ChIP4[,1] >= pericenStart[i] & WT_SPO11_ChIP4[,1] <= pericenEnd[i],]
  WT_SPO11_ChIP4_peri <- rbind(WT_SPO11_ChIP4_peri, WT_SPO11_ChIP4_peri_tmp)
  WT_SPO11_ChIP4_armL <- rbind(WT_SPO11_ChIP4_armL, WT_SPO11_ChIP4_armL_tmp)
  WT_SPO11_ChIP4_armR <- rbind(WT_SPO11_ChIP4_armR, WT_SPO11_ChIP4_armR_tmp)
  WT_SPO11_ChIP4_arm <- rbind(WT_SPO11_ChIP4_arm, WT_SPO11_ChIP4_arm_tmp)

  print(i)
  WT_SPO11_oligos_RPI1_armL_tmp <- WT_SPO11_oligos_RPI1[WT_SPO11_oligos_RPI1[,1] > sumchr[i] & WT_SPO11_oligos_RPI1[,1] < pericenStart[i] & WT_SPO11_oligos_RPI1[,1] < sumchr[i+1],]
  WT_SPO11_oligos_RPI1_armR_tmp <- WT_SPO11_oligos_RPI1[WT_SPO11_oligos_RPI1[,1] > pericenEnd[i] & WT_SPO11_oligos_RPI1[,1] <= sumchr[i+1],]
  WT_SPO11_oligos_RPI1_arm_tmp <- WT_SPO11_oligos_RPI1[(WT_SPO11_oligos_RPI1[,1] > sumchr[i] & WT_SPO11_oligos_RPI1[,1] < pericenStart[i] & WT_SPO11_oligos_RPI1[,1] < sumchr[i+1]) |
                                 (WT_SPO11_oligos_RPI1[,1] > pericenEnd[i] & WT_SPO11_oligos_RPI1[,1] <= sumchr[i+1]),]
  WT_SPO11_oligos_RPI1_peri_tmp <- WT_SPO11_oligos_RPI1[WT_SPO11_oligos_RPI1[,1] >= pericenStart[i] & WT_SPO11_oligos_RPI1[,1] <= pericenEnd[i],]
  WT_SPO11_oligos_RPI1_peri <- rbind(WT_SPO11_oligos_RPI1_peri, WT_SPO11_oligos_RPI1_peri_tmp)
  WT_SPO11_oligos_RPI1_armL <- rbind(WT_SPO11_oligos_RPI1_armL, WT_SPO11_oligos_RPI1_armL_tmp)
  WT_SPO11_oligos_RPI1_armR <- rbind(WT_SPO11_oligos_RPI1_armR, WT_SPO11_oligos_RPI1_armR_tmp)
  WT_SPO11_oligos_RPI1_arm <- rbind(WT_SPO11_oligos_RPI1_arm, WT_SPO11_oligos_RPI1_arm_tmp)
  print(i)
  WT_nakedDNA_T50R1_armL_tmp <- WT_nakedDNA_T50R1[WT_nakedDNA_T50R1[,1] > sumchr[i] & WT_nakedDNA_T50R1[,1] < pericenStart[i] & WT_nakedDNA_T50R1[,1] < sumchr[i+1],]
  WT_nakedDNA_T50R1_armR_tmp <- WT_nakedDNA_T50R1[WT_nakedDNA_T50R1[,1] > pericenEnd[i] & WT_nakedDNA_T50R1[,1] <= sumchr[i+1],]
  WT_nakedDNA_T50R1_arm_tmp <- WT_nakedDNA_T50R1[(WT_nakedDNA_T50R1[,1] > sumchr[i] & WT_nakedDNA_T50R1[,1] < pericenStart[i] & WT_nakedDNA_T50R1[,1] < sumchr[i+1]) |
                                 (WT_nakedDNA_T50R1[,1] > pericenEnd[i] & WT_nakedDNA_T50R1[,1] <= sumchr[i+1]),]
  WT_nakedDNA_T50R1_peri_tmp <- WT_nakedDNA_T50R1[WT_nakedDNA_T50R1[,1] >= pericenStart[i] & WT_nakedDNA_T50R1[,1] <= pericenEnd[i],]
  WT_nakedDNA_T50R1_peri <- rbind(WT_nakedDNA_T50R1_peri, WT_nakedDNA_T50R1_peri_tmp)
  WT_nakedDNA_T50R1_armL <- rbind(WT_nakedDNA_T50R1_armL, WT_nakedDNA_T50R1_armL_tmp)
  WT_nakedDNA_T50R1_armR <- rbind(WT_nakedDNA_T50R1_armR, WT_nakedDNA_T50R1_armR_tmp)
  WT_nakedDNA_T50R1_arm <- rbind(WT_nakedDNA_T50R1_arm, WT_nakedDNA_T50R1_arm_tmp)

  print(i)
  cumCODat_armL_tmp <- cumCODat[cumCODat[,1] > sumchr[i] & cumCODat[,1] < pericenStart[i] & cumCODat[,1] < sumchr[i+1],]
  cumCODat_armR_tmp <- cumCODat[cumCODat[,1] > pericenEnd[i] & cumCODat[,1] <= sumchr[i+1],]
  cumCODat_arm_tmp <- cumCODat[(cumCODat[,1] > sumchr[i] & cumCODat[,1] < pericenStart[i] & cumCODat[,1] < sumchr[i+1]) |
                                 (cumCODat[,1] > pericenEnd[i] & cumCODat[,1] <= sumchr[i+1]),]
  cumCODat_peri_tmp <- cumCODat[cumCODat[,1] >= pericenStart[i] & cumCODat[,1] <= pericenEnd[i],]
  cumCODat_peri <- rbind(cumCODat_peri, cumCODat_peri_tmp)
  cumCODat_armL <- rbind(cumCODat_armL, cumCODat_armL_tmp)
  cumCODat_armR <- rbind(cumCODat_armR, cumCODat_armR_tmp)
  cumCODat_arm <- rbind(cumCODat_arm, cumCODat_arm_tmp)

  print(i)
  cumCODatRed_armL_tmp <- cumCODatRed[cumCODatRed[,1] > sumchr[i] & cumCODatRed[,1] < pericenStart[i] & cumCODatRed[,1] < sumchr[i+1],]
  cumCODatRed_armR_tmp <- cumCODatRed[cumCODatRed[,1] > pericenEnd[i] & cumCODatRed[,1] <= sumchr[i+1],]
  cumCODatRed_arm_tmp <- cumCODatRed[(cumCODatRed[,1] > sumchr[i] & cumCODatRed[,1] < pericenStart[i] & cumCODatRed[,1] < sumchr[i+1]) |
                                 (cumCODatRed[,1] > pericenEnd[i] & cumCODatRed[,1] <= sumchr[i+1]),]
  cumCODatRed_peri_tmp <- cumCODatRed[cumCODatRed[,1] >= pericenStart[i] & cumCODatRed[,1] <= pericenEnd[i],]
  cumCODatRed_peri <- rbind(cumCODatRed_peri, cumCODatRed_peri_tmp)
  cumCODatRed_armL <- rbind(cumCODatRed_armL, cumCODatRed_armL_tmp)
  cumCODatRed_armR <- rbind(cumCODatRed_armR, cumCODatRed_armR_tmp)
  cumCODatRed_arm <- rbind(cumCODatRed_arm, cumCODatRed_arm_tmp)

  print(i)
  MSH4_ChIP_armL_tmp <- MSH4_ChIP[MSH4_ChIP[,1] > sumchr[i] & MSH4_ChIP[,1] < pericenStart[i] & MSH4_ChIP[,1] < sumchr[i+1],]
  MSH4_ChIP_armR_tmp <- MSH4_ChIP[MSH4_ChIP[,1] > pericenEnd[i] & MSH4_ChIP[,1] <= sumchr[i+1],]
  MSH4_ChIP_arm_tmp <- MSH4_ChIP[(MSH4_ChIP[,1] > sumchr[i] & MSH4_ChIP[,1] < pericenStart[i] & MSH4_ChIP[,1] < sumchr[i+1]) |
                                 (MSH4_ChIP[,1] > pericenEnd[i] & MSH4_ChIP[,1] <= sumchr[i+1]),]
  MSH4_ChIP_peri_tmp <- MSH4_ChIP[MSH4_ChIP[,1] >= pericenStart[i] & MSH4_ChIP[,1] <= pericenEnd[i],]
  MSH4_ChIP_peri <- rbind(MSH4_ChIP_peri, MSH4_ChIP_peri_tmp)
  MSH4_ChIP_armL <- rbind(MSH4_ChIP_armL, MSH4_ChIP_armL_tmp)
  MSH4_ChIP_armR <- rbind(MSH4_ChIP_armR, MSH4_ChIP_armR_tmp)
  MSH4_ChIP_arm <- rbind(MSH4_ChIP_arm, MSH4_ChIP_arm_tmp)
  print(i)
  MSH4_input_armL_tmp <- MSH4_input[MSH4_input[,1] > sumchr[i] & MSH4_input[,1] < pericenStart[i] & MSH4_input[,1] < sumchr[i+1],]
  MSH4_input_armR_tmp <- MSH4_input[MSH4_input[,1] > pericenEnd[i] & MSH4_input[,1] <= sumchr[i+1],]
  MSH4_input_arm_tmp <- MSH4_input[(MSH4_input[,1] > sumchr[i] & MSH4_input[,1] < pericenStart[i] & MSH4_input[,1] < sumchr[i+1]) |
                                 (MSH4_input[,1] > pericenEnd[i] & MSH4_input[,1] <= sumchr[i+1]),]
  MSH4_input_peri_tmp <- MSH4_input[MSH4_input[,1] >= pericenStart[i] & MSH4_input[,1] <= pericenEnd[i],]
  MSH4_input_peri <- rbind(MSH4_input_peri, MSH4_input_peri_tmp)
  MSH4_input_armL <- rbind(MSH4_input_armL, MSH4_input_armL_tmp)
  MSH4_input_armR <- rbind(MSH4_input_armR, MSH4_input_armR_tmp)
  MSH4_input_arm <- rbind(MSH4_input_arm, MSH4_input_arm_tmp)
}

# log2-normalised REC8 and other coverage profiles (genome-scale)
  REC8_arm_norm <- log2(REC8_ChIP_arm[,2]/REC8_input_arm[,2])
  #REC8_arm_norm <- (REC8_arm_norm-mean(REC8_arm_norm, na.rm = T))/sd(REC8_arm_norm, na.rm = T)
  REC8_peri_norm <- log2(REC8_ChIP_peri[,2]/REC8_input_peri[,2])
  #REC8_peri_norm <- (REC8_peri_norm-mean(REC8_peri_norm, na.rm = T))/sd(REC8_peri_norm, na.rm = T)

  nuc_arm_norm <- log2(WT_nuc_arm[,2]/WT_nakedDNA_arm[,2])
  #nuc_arm_norm <- (nuc_arm_norm-mean(nuc_arm_norm, na.rm = T))/sd(nuc_arm_norm, na.rm = T)
  nuc_peri_norm <- log2(WT_nuc_peri[,2]/WT_nakedDNA_peri[,2])
  #nuc_peri_norm <- (nuc_peri_norm-mean(nuc_peri_norm, na.rm = T))/sd(nuc_peri_norm, na.rm = T)

  H2AZ_arm_norm <- log2(H2AZ_ChIP_arm[,2]/H2A_input_arm[,2])
  #H2AZ_arm_norm <- (H2AZ_arm_norm-mean(H2AZ_arm_norm, na.rm = T))/sd(H2AZ_arm_norm, na.rm = T)
  H2AZ_peri_norm <- log2(H2AZ_ChIP_peri[,2]/H2A_input_peri[,2])
  #H2AZ_peri_norm <- (H2AZ_peri_norm-mean(H2AZ_peri_norm, na.rm = T))/sd(H2AZ_peri_norm, na.rm = T)

  H3K4me3_arm_norm <- log2(WT_H3K4me3_ChIP_arm[,2]/WT_H3K4me3_input_arm[,2])
  #H3K4me3_arm_norm <- (H3K4me3_arm_norm-mean(H3K4me3_arm_norm, na.rm = T))/sd(H3K4me3_arm_norm, na.rm = T)
  H3K4me3_peri_norm <- log2(WT_H3K4me3_ChIP_peri[,2]/WT_H3K4me3_input_peri[,2])
  #H3K4me3_peri_norm <- (H3K4me3_peri_norm-mean(H3K4me3_peri_norm, na.rm = T))/sd(H3K4me3_peri_norm, na.rm = T)

  H2AW_arm_norm <- log2(H2AW_ChIP_arm[,2]/H2A_input_arm[,2])
  #H2AW_arm_norm <- (H2AW_arm_norm-mean(H2AW_arm_norm, na.rm = T))/sd(H2AW_arm_norm, na.rm = T)
  H2AW_peri_norm <- log2(H2AW_ChIP_peri[,2]/H2A_input_peri[,2])
  #H2AW_peri_norm <- (H2AW_peri_norm-mean(H2AW_peri_norm, na.rm = T))/sd(H2AW_peri_norm, na.rm = T)
 
  H3K9me2_arm_norm <- log2(WT_H3K9me2_ChIP_arm[,2]/WT_H3K9me2_input_arm[,2])
  #H3K9me2_arm_norm <- (H3K9me2_arm_norm-mean(H3K9me2_arm_norm, na.rm = T))/sd(H3K9me2_arm_norm, na.rm = T)
  H3K9me2_peri_norm <- log2(WT_H3K9me2_ChIP_peri[,2]/WT_H3K9me2_input_peri[,2])
  #H3K9me2_peri_norm <- (H3K9me2_peri_norm-mean(H3K9me2_peri_norm, na.rm = T))/sd(H3K9me2_peri_norm, na.rm = T)

  SPO11_ChIP4_arm_norm <- log2(WT_SPO11_ChIP4_arm[,2]/REC8_input_arm[,2])
  #SPO11_ChIP4_arm_norm <- (SPO11_ChIP4_arm_norm-mean(SPO11_ChIP4_arm_norm, na.rm = T))/sd(SPO11_ChIP4_arm_norm, na.rm = T)
  SPO11_ChIP4_peri_norm <- log2(WT_SPO11_ChIP4_peri[,2]/REC8_input_peri[,2])
  #SPO11_ChIP4_peri_norm <- (SPO11_ChIP4_peri_norm-mean(SPO11_ChIP4_peri_norm, na.rm = T))/sd(SPO11_ChIP4_peri_norm, na.rm = T)

  SPO11_oligos_RPI1_arm_norm <- log2(WT_SPO11_oligos_RPI1_arm[,2]/WT_nakedDNA_T50R1_arm[,2])
  #SPO11_oligos_RPI1_arm_norm <- (SPO11_oligos_RPI1_arm_norm-mean(SPO11_oligos_RPI1_arm_norm, na.rm = T))/sd(SPO11_oligos_RPI1_arm_norm, na.rm = T)
  SPO11_oligos_RPI1_peri_norm <- log2(WT_SPO11_oligos_RPI1_peri[,2]/WT_nakedDNA_T50R1_peri[,2])
  #SPO11_oligos_RPI1_peri_norm <- (SPO11_oligos_RPI1_peri_norm-mean(SPO11_oligos_RPI1_peri_norm, na.rm = T))/sd(SPO11_oligos_RPI1_peri_norm, na.rm = T)

  MSH4_arm_norm <- log2(MSH4_ChIP_arm[,2]/MSH4_input_arm[,2])
  #MSH4_arm_norm <- (MSH4_arm_norm-mean(MSH4_arm_norm, na.rm = T))/sd(MSH4_arm_norm, na.rm = T)
  MSH4_peri_norm <- log2(MSH4_ChIP_peri[,2]/MSH4_input_peri[,2])
  #MSH4_peri_norm <- (MSH4_peri_norm-mean(MSH4_peri_norm, na.rm = T))/sd(MSH4_peri_norm, na.rm = T)

# NOT SMOOTHED with MA filter
allDF_arm <- data.frame(REC8_arm_norm, nuc_arm_norm, H2AZ_arm_norm, H3K4me3_arm_norm, cumGeneDat_arm$winGenes, H2AW_arm_norm, H3K9me2_arm_norm, cumTEDat_arm$winTEs, cumMethDat10_arm$wtCG_win_vals, cumMethDat10_arm$wtCHG_win_vals, cumMethDat10_arm$wtCHH_win_vals, cumMethDat10_arm$allcntxts_win_vals)
colnames(allDF_arm) <- c("REC8", "MNase", "H2A.Z", "H3K4me3", "Genes", "H2A.W", "H3K9me2", "TEs", "CG", "CHG", "CHH", "DNA meth")
allDF_arm_corMat <- cor(allDF_arm, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "arm_correlation_matrix_", winNames[1], "_colouronly_unsmoothed_noZscore.pdf"))
corrplot(allDF_arm_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0),
         title = paste0("Arms Spearman correlation matrix (10-kb windows)"))
dev.off()

allDF_peri <- data.frame(REC8_peri_norm, nuc_peri_norm, H2AZ_peri_norm, H3K4me3_peri_norm, cumGeneDat_peri$winGenes, H2AW_peri_norm, H3K9me2_peri_norm, cumTEDat_peri$winTEs, cumMethDat10_peri$wtCG_win_vals, cumMethDat10_peri$wtCHG_win_vals, cumMethDat10_peri$wtCHH_win_vals, cumMethDat10_peri$allcntxts_win_vals)
colnames(allDF_peri) <- c("REC8", "MNase", "H2A.Z", "H3K4me3", "Genes", "H2A.W", "H3K9me2", "TEs", "CG", "CHG", "CHH", "DNA meth")
allDF_peri_corMat <- cor(allDF_peri, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "peri_correlation_matrix_", winNames[1], "_colouronly_unsmoothed_noZscore.pdf"))
corrplot(allDF_peri_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0),
         title = paste0("Pericentromeres Spearman correlation matrix (10-kb windows)"))
dev.off()


# SMOOTHED with MA filter
# ARMS
test <- seq(1, 1000, by = 1)
j = 100
ma <- rep(1, test[j])/test[j]
filt_REC8_arm_norm <- stats::filter(REC8_arm_norm, ma)
which_na <- which(is.na(filt_REC8_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_REC8_arm_norm[left_na[length(left_na)]+1]
filt_REC8_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_REC8_arm_norm[right_na[1]-1]
filt_REC8_arm_norm[right_na] <- right_val
filt_REC8_arm_norm_noNA <- filt_REC8_arm_norm[!is.na(filt_REC8_arm_norm)]

filt_nuc_arm_norm <- stats::filter(nuc_arm_norm, ma)
which_na <- which(is.na(filt_nuc_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_nuc_arm_norm[left_na[length(left_na)]+1]
filt_nuc_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_nuc_arm_norm[right_na[1]-1]
filt_nuc_arm_norm[right_na] <- right_val
filt_nuc_arm_norm_noNA <- filt_nuc_arm_norm[!is.na(filt_nuc_arm_norm)]

filt_H2AZ_arm_norm <- stats::filter(H2AZ_arm_norm, ma)
which_na <- which(is.na(filt_H2AZ_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AZ_arm_norm[left_na[length(left_na)]+1]
filt_H2AZ_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AZ_arm_norm[right_na[1]-1]
filt_H2AZ_arm_norm[right_na] <- right_val
filt_H2AZ_arm_norm_noNA <- filt_H2AZ_arm_norm[!is.na(filt_H2AZ_arm_norm)]

filt_H3K4me3_arm_norm <- stats::filter(H3K4me3_arm_norm, ma)
which_na <- which(is.na(filt_H3K4me3_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H3K4me3_arm_norm[left_na[length(left_na)]+1]
filt_H3K4me3_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H3K4me3_arm_norm[right_na[1]-1]
filt_H3K4me3_arm_norm[right_na] <- right_val
filt_H3K4me3_arm_norm_noNA <- filt_H3K4me3_arm_norm[!is.na(filt_H3K4me3_arm_norm)]

filt_cumGeneDat_arm <- stats::filter(cumGeneDat_arm[,2], ma)
which_na <- which(is.na(filt_cumGeneDat_arm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumGeneDat_arm[left_na[length(left_na)]+1]
filt_cumGeneDat_arm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumGeneDat_arm[right_na[1]-1]
filt_cumGeneDat_arm[right_na] <- right_val
filt_cumGeneDat_arm_noNA <- filt_cumGeneDat_arm[!is.na(filt_cumGeneDat_arm)]

filt_wtCG10_arm <- stats::filter(cumMethDat10_arm[,2], ma); filt_wtCG10_arm_n <- filt_wtCG10_arm[!is.na(filt_wtCG10_arm)]
filt_wtCHG10_arm <- stats::filter(cumMethDat10_arm[,3], ma); filt_wtCHG10_arm_n <- filt_wtCHG10_arm[!is.na(filt_wtCHG10_arm)]
filt_wtCHH10_arm <- stats::filter(cumMethDat10_arm[,4], ma); filt_wtCHH10_arm_n <- filt_wtCHH10_arm[!is.na(filt_wtCHH10_arm)]
filt_allcntxts10_arm <- stats::filter(cumMethDat10_arm[,5], ma); filt_allcntxts10_arm_n <- filt_allcntxts10_arm[!is.na(filt_allcntxts10_arm)]
filt_wtMeth10_arm <- cbind(cumMethDat10_arm[,1], filt_wtCG10_arm, filt_wtCHG10_arm, filt_wtCHH10_arm, filt_allcntxts10_arm)
write.table(filt_wtMeth10_arm, file = paste0(outDir, "filt_wtMeth_arm_genome_10kb.txt"))
#filt_wtMeth10_arm <- read.table(file = paste0(outDir, "filt_wtMeth_arm_genome_10kb.txt"))

filt_H2AW_arm_norm <- stats::filter(H2AW_arm_norm, ma)
which_na <- which(is.na(filt_H2AW_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AW_arm_norm[left_na[length(left_na)]+1]
filt_H2AW_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AW_arm_norm[right_na[1]-1]
filt_H2AW_arm_norm[right_na] <- right_val
filt_H2AW_arm_norm_noNA <- filt_H2AW_arm_norm[!is.na(filt_H2AW_arm_norm)]

filt_H3K9me2_arm_norm <- stats::filter(H3K9me2_arm_norm, ma)
which_na <- which(is.na(filt_H3K9me2_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H3K9me2_arm_norm[left_na[length(left_na)]+1]
filt_H3K9me2_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H3K9me2_arm_norm[right_na[1]-1]
filt_H3K9me2_arm_norm[right_na] <- right_val
filt_H3K9me2_arm_norm_noNA <- filt_H3K9me2_arm_norm[!is.na(filt_H3K9me2_arm_norm)]

filt_cumTEDat_arm <- stats::filter(cumTEDat_arm[,2], ma)
which_na <- which(is.na(filt_cumTEDat_arm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumTEDat_arm[left_na[length(left_na)]+1]
filt_cumTEDat_arm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumTEDat_arm[right_na[1]-1]
filt_cumTEDat_arm[right_na] <- right_val
filt_cumTEDat_arm_noNA <- filt_cumTEDat_arm[!is.na(filt_cumTEDat_arm)]

filt_SPO11_ChIP4_arm_norm <- stats::filter(SPO11_ChIP4_arm_norm, ma)
which_na <- which(is.na(filt_SPO11_ChIP4_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_SPO11_ChIP4_arm_norm[left_na[length(left_na)]+1]
filt_SPO11_ChIP4_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_SPO11_ChIP4_arm_norm[right_na[1]-1]
filt_SPO11_ChIP4_arm_norm[right_na] <- right_val
filt_SPO11_ChIP4_arm_norm_noNA <- filt_SPO11_ChIP4_arm_norm[!is.na(filt_SPO11_ChIP4_arm_norm)]

filt_SPO11_oligos_RPI1_arm_norm <- stats::filter(SPO11_oligos_RPI1_arm_norm, ma)
which_na <- which(is.na(filt_SPO11_oligos_RPI1_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_SPO11_oligos_RPI1_arm_norm[left_na[length(left_na)]+1]
filt_SPO11_oligos_RPI1_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_SPO11_oligos_RPI1_arm_norm[right_na[1]-1]
filt_SPO11_oligos_RPI1_arm_norm[right_na] <- right_val
filt_SPO11_oligos_RPI1_arm_norm_noNA <- filt_SPO11_oligos_RPI1_arm_norm[!is.na(filt_SPO11_oligos_RPI1_arm_norm)]

filt_MSH4_arm_norm <- stats::filter(MSH4_arm_norm, ma)
which_na <- which(is.na(filt_MSH4_arm_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_MSH4_arm_norm[left_na[length(left_na)]+1]
filt_MSH4_arm_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_MSH4_arm_norm[right_na[1]-1]
filt_MSH4_arm_norm[right_na] <- right_val
filt_MSH4_arm_norm_noNA <- filt_MSH4_arm_norm[!is.na(filt_MSH4_arm_norm)]

j = 200
ma <- rep(1, test[j])/test[j]
filt_cumCODatRed_arm <- stats::filter(cumCODatRed_arm[,2], ma)
which_na <- which(is.na(filt_cumCODatRed_arm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumCODatRed_arm[left_na[length(left_na)]+1]
filt_cumCODatRed_arm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumCODatRed_arm[right_na[1]-1]
filt_cumCODatRed_arm[right_na] <- right_val
filt_cumCODatRed_arm_noNA <- filt_cumCODatRed_arm[!is.na(filt_cumCODatRed_arm)]

filt_cumCODat_arm <- stats::filter(cumCODat_arm[,2], ma)
which_na <- which(is.na(filt_cumCODat_arm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumCODat_arm[left_na[length(left_na)]+1]
filt_cumCODat_arm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumCODat_arm[right_na[1]-1]
filt_cumCODat_arm[right_na] <- right_val
filt_cumCODat_arm_noNA <- filt_cumCODat_arm[!is.na(filt_cumCODat_arm)]

# PERICENTROMERES
test <- seq(1, 1000, by = 1)
j = 100
ma <- rep(1, test[j])/test[j]
filt_REC8_peri_norm <- stats::filter(REC8_peri_norm, ma)
which_na <- which(is.na(filt_REC8_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_REC8_peri_norm[left_na[length(left_na)]+1]
filt_REC8_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_REC8_peri_norm[right_na[1]-1]
filt_REC8_peri_norm[right_na] <- right_val
filt_REC8_peri_norm_noNA <- filt_REC8_peri_norm[!is.na(filt_REC8_peri_norm)]

filt_nuc_peri_norm <- stats::filter(nuc_peri_norm, ma)
which_na <- which(is.na(filt_nuc_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_nuc_peri_norm[left_na[length(left_na)]+1]
filt_nuc_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_nuc_peri_norm[right_na[1]-1]
filt_nuc_peri_norm[right_na] <- right_val
filt_nuc_peri_norm_noNA <- filt_nuc_peri_norm[!is.na(filt_nuc_peri_norm)]

filt_H2AZ_peri_norm <- stats::filter(H2AZ_peri_norm, ma)
which_na <- which(is.na(filt_H2AZ_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AZ_peri_norm[left_na[length(left_na)]+1]
filt_H2AZ_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AZ_peri_norm[right_na[1]-1]
filt_H2AZ_peri_norm[right_na] <- right_val
filt_H2AZ_peri_norm_noNA <- filt_H2AZ_peri_norm[!is.na(filt_H2AZ_peri_norm)]

filt_H3K4me3_peri_norm <- stats::filter(H3K4me3_peri_norm, ma)
which_na <- which(is.na(filt_H3K4me3_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H3K4me3_peri_norm[left_na[length(left_na)]+1]
filt_H3K4me3_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H3K4me3_peri_norm[right_na[1]-1]
filt_H3K4me3_peri_norm[right_na] <- right_val
filt_H3K4me3_peri_norm_noNA <- filt_H3K4me3_peri_norm[!is.na(filt_H3K4me3_peri_norm)]

filt_cumGeneDat_peri <- stats::filter(cumGeneDat_peri[,2], ma)
which_na <- which(is.na(filt_cumGeneDat_peri) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumGeneDat_peri[left_na[length(left_na)]+1]
filt_cumGeneDat_peri[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumGeneDat_peri[right_na[1]-1]
filt_cumGeneDat_peri[right_na] <- right_val
filt_cumGeneDat_peri_noNA <- filt_cumGeneDat_peri[!is.na(filt_cumGeneDat_peri)]

filt_wtCG10_peri <- stats::filter(cumMethDat10_peri[,2], ma); filt_wtCG10_peri_n <- filt_wtCG10_peri[!is.na(filt_wtCG10_peri)]
filt_wtCHG10_peri <- stats::filter(cumMethDat10_peri[,3], ma); filt_wtCHG10_peri_n <- filt_wtCHG10_peri[!is.na(filt_wtCHG10_peri)]
filt_wtCHH10_peri <- stats::filter(cumMethDat10_peri[,4], ma); filt_wtCHH10_peri_n <- filt_wtCHH10_peri[!is.na(filt_wtCHH10_peri)]
filt_allcntxts10_peri <- stats::filter(cumMethDat10_peri[,5], ma); filt_allcntxts10_peri_n <- filt_allcntxts10_peri[!is.na(filt_allcntxts10_peri)]
filt_wtMeth10_peri <- cbind(cumMethDat10_peri[,1], filt_wtCG10_peri, filt_wtCHG10_peri, filt_wtCHH10_peri, filt_allcntxts10_peri)
write.table(filt_wtMeth10_peri, file = paste0(outDir, "filt_wtMeth_peri_genome_10kb.txt"))
#filt_wtMeth10_peri <- read.table(file = paste0(outDir, "filt_wtMeth_peri_genome_10kb.txt"))

filt_H2AW_peri_norm <- stats::filter(H2AW_peri_norm, ma)
which_na <- which(is.na(filt_H2AW_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H2AW_peri_norm[left_na[length(left_na)]+1]
filt_H2AW_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H2AW_peri_norm[right_na[1]-1]
filt_H2AW_peri_norm[right_na] <- right_val
filt_H2AW_peri_norm_noNA <- filt_H2AW_peri_norm[!is.na(filt_H2AW_peri_norm)]

filt_H3K9me2_peri_norm <- stats::filter(H3K9me2_peri_norm, ma)
which_na <- which(is.na(filt_H3K9me2_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_H3K9me2_peri_norm[left_na[length(left_na)]+1]
filt_H3K9me2_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_H3K9me2_peri_norm[right_na[1]-1]
filt_H3K9me2_peri_norm[right_na] <- right_val
filt_H3K9me2_peri_norm_noNA <- filt_H3K9me2_peri_norm[!is.na(filt_H3K9me2_peri_norm)]

filt_cumTEDat_peri <- stats::filter(cumTEDat_peri[,2], ma)
which_na <- which(is.na(filt_cumTEDat_peri) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumTEDat_peri[left_na[length(left_na)]+1]
filt_cumTEDat_peri[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumTEDat_peri[right_na[1]-1]
filt_cumTEDat_peri[right_na] <- right_val
filt_cumTEDat_peri_noNA <- filt_cumTEDat_peri[!is.na(filt_cumTEDat_peri)]

filt_SPO11_ChIP4_peri_norm <- stats::filter(SPO11_ChIP4_peri_norm, ma)
which_na <- which(is.na(filt_SPO11_ChIP4_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_SPO11_ChIP4_peri_norm[left_na[length(left_na)]+1]
filt_SPO11_ChIP4_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_SPO11_ChIP4_peri_norm[right_na[1]-1]
filt_SPO11_ChIP4_peri_norm[right_na] <- right_val
filt_SPO11_ChIP4_peri_norm_noNA <- filt_SPO11_ChIP4_peri_norm[!is.na(filt_SPO11_ChIP4_peri_norm)]

filt_SPO11_oligos_RPI1_peri_norm <- stats::filter(SPO11_oligos_RPI1_peri_norm, ma)
which_na <- which(is.na(filt_SPO11_oligos_RPI1_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_SPO11_oligos_RPI1_peri_norm[left_na[length(left_na)]+1]
filt_SPO11_oligos_RPI1_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_SPO11_oligos_RPI1_peri_norm[right_na[1]-1]
filt_SPO11_oligos_RPI1_peri_norm[right_na] <- right_val
filt_SPO11_oligos_RPI1_peri_norm_noNA <- filt_SPO11_oligos_RPI1_peri_norm[!is.na(filt_SPO11_oligos_RPI1_peri_norm)]

filt_MSH4_peri_norm <- stats::filter(MSH4_peri_norm, ma)
which_na <- which(is.na(filt_MSH4_peri_norm) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_MSH4_peri_norm[left_na[length(left_na)]+1]
filt_MSH4_peri_norm[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_MSH4_peri_norm[right_na[1]-1]
filt_MSH4_peri_norm[right_na] <- right_val
filt_MSH4_peri_norm_noNA <- filt_MSH4_peri_norm[!is.na(filt_MSH4_peri_norm)]

j = 200
ma <- rep(1, test[j])/test[j]
filt_cumCODatRed_peri <- stats::filter(cumCODatRed_peri[,2], ma)
which_na <- which(is.na(filt_cumCODatRed_peri) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumCODatRed_peri[left_na[length(left_na)]+1]
filt_cumCODatRed_peri[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumCODatRed_peri[right_na[1]-1]
filt_cumCODatRed_peri[right_na] <- right_val
filt_cumCODatRed_peri_noNA <- filt_cumCODatRed_peri[!is.na(filt_cumCODatRed_peri)]

filt_cumCODat_peri <- stats::filter(cumCODat_peri[,2], ma)
which_na <- which(is.na(filt_cumCODat_peri) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_cumCODat_peri[left_na[length(left_na)]+1]
filt_cumCODat_peri[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_cumCODat_peri[right_na[1]-1]
filt_cumCODat_peri[right_na] <- right_val
filt_cumCODat_peri_noNA <- filt_cumCODat_peri[!is.na(filt_cumCODat_peri)]

# Create arm and peri Spearman's rho correlation matrices
allDF_arm <- data.frame(filt_REC8_arm_norm, filt_nuc_arm_norm, filt_H2AZ_arm_norm, filt_H3K4me3_arm_norm, filt_cumGeneDat_arm, filt_H2AW_arm_norm, filt_H3K9me2_arm_norm, filt_cumTEDat_arm, filt_wtCG10_arm, filt_wtCHG10_arm, filt_wtCHH10_arm, filt_allcntxts10_arm)
colnames(allDF_arm) <- c("REC8", "MNase", "H2A.Z", "H3K4me3", "Genes", "H2A.W", "H3K9me2", "TEs", "CG", "CHG", "CHH", "DNA meth")
allDF_arm_corMat <- cor(allDF_arm, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "arm_correlation_matrix_", winNames[1], "_colouronly_noZscore.pdf"))
corrplot(allDF_arm_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0), tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.8,
         title = paste0("Arms Spearman correlation matrix (10-kb windows)"))
dev.off()

allDF_peri <- data.frame(filt_REC8_peri_norm, filt_nuc_peri_norm, filt_H2AZ_peri_norm, filt_H3K4me3_peri_norm, filt_cumGeneDat_peri, filt_H2AW_peri_norm, filt_H3K9me2_peri_norm, filt_cumTEDat_peri, filt_wtCG10_peri, filt_wtCHG10_peri, filt_wtCHH10_peri, filt_allcntxts10_peri)
colnames(allDF_peri) <- c("REC8", "MNase", "H2A.Z", "H3K4me3", "Genes", "H2A.W", "H3K9me2", "TEs", "CG", "CHG", "CHH", "DNA meth")
allDF_peri_corMat <- cor(allDF_peri, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "peri_correlation_matrix_", winNames[1], "_colouronly_noZscore.pdf"))
corrplot(allDF_peri_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0), tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.8,
         title = paste0("Pericentromeres Spearman correlation matrix (10-kb windows)"))
dev.off()

# Create arm and peri Spearman's rho correlation matrices including SPO11-1-oligos, SPO11-1 ChIP, crossovers, and MSH4
allDF_arm <- data.frame(filt_REC8_arm_norm, filt_nuc_arm_norm, filt_SPO11_ChIP4_arm_norm, filt_MSH4_arm_norm, filt_H2AW_arm_norm, filt_allcntxts10_arm, filt_H3K9me2_arm_norm, filt_cumGeneDat_arm, filt_H2AZ_arm_norm, filt_H3K4me3_arm_norm, filt_SPO11_oligos_RPI1_arm_norm, filt_cumCODatRed_arm)
colnames(allDF_arm) <- c("REC8", "MNase", "SPO11-1 ChIP", "MSH4", "H2A.W", "DNA meth", "H3K9me2", "Genes", "H2A.Z", "H3K4me3", "SPO11-1-oligos", "Crossovers")
allDF_arm_corMat <- cor(allDF_arm, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "arm_correlation_matrix_", winNames[1], "_colouronly_v02_noZscore.pdf"))
corrplot(allDF_arm_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0), tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.8,
         title = paste0("Arms Spearman correlation matrix (10-kb windows)"))
dev.off()

allDF_peri <- data.frame(filt_REC8_peri_norm, filt_nuc_peri_norm, filt_SPO11_ChIP4_peri_norm, filt_MSH4_peri_norm, filt_H2AW_peri_norm, filt_allcntxts10_peri, filt_H3K9me2_peri_norm, filt_cumGeneDat_peri, filt_H2AZ_peri_norm, filt_H3K4me3_peri_norm, filt_SPO11_oligos_RPI1_peri_norm, filt_cumCODatRed_peri)
colnames(allDF_peri) <- c("REC8", "MNase", "SPO11-1 ChIP", "MSH4", "H2A.W", "DNA meth", "H3K9me2", "Genes", "H2A.Z", "H3K4me3", "SPO11-1-oligos", "Crossovers")
allDF_peri_corMat <- cor(allDF_peri, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "peri_correlation_matrix_", winNames[1], "_colouronly_v02_noZscore.pdf"))
corrplot(allDF_peri_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,1,0), tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.8,
         title = paste0("Pericentromeres Spearman correlation matrix (10-kb windows)"))
dev.off()


# For Chris's poster, create genome-wide, arm and peri correlation matrices including REC8, SPO11-ChIP, SPO11-oligo and crossovers only

smallDF <- data.frame(filt_REC8_norm, filt_SPO11_ChIP4_norm, filt_SPO11_oligos_RPI1_norm, filt_cumCODatRed)
colnames(smallDF) <- c("REC8", "SPO11-1 ChIP", "SPO11-1-oligos", "Crossovers")
smallDF_corMat <- cor(smallDF, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "genome-wide_correlation_matrix_REC8_SPO11_crossovers_", winNames[1], "_colouronly_noZscore.pdf"))
corrplot(smallDF_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,0,0), tl.cex = 0.8)
dev.off()

smallDF_arm <- data.frame(filt_REC8_arm_norm, filt_SPO11_ChIP4_arm_norm, filt_SPO11_oligos_RPI1_arm_norm, filt_cumCODatRed_arm)
colnames(smallDF_arm) <- c("REC8", "SPO11-1 ChIP", "SPO11-1-oligos", "Crossovers")
smallDF_arm_corMat <- cor(smallDF_arm, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "arm_correlation_matrix_REC8_SPO11_crossovers_", winNames[1], "_colouronly_noZscore.pdf"))
corrplot(smallDF_arm_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,0,0), tl.cex = 0.8)
dev.off()

smallDF_peri <- data.frame(filt_REC8_peri_norm, filt_SPO11_ChIP4_peri_norm, filt_SPO11_oligos_RPI1_peri_norm, filt_cumCODatRed_peri)
colnames(smallDF_peri) <- c("REC8", "SPO11-1 ChIP", "SPO11-1-oligos", "Crossovers")
smallDF_peri_corMat <- cor(smallDF_peri, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(outDir, "peri_correlation_matrix_REC8_SPO11_crossovers_", winNames[1], "_colouronly_noZscore.pdf"))
corrplot(smallDF_peri_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "white", mar = c(0,0,0,0), tl.cex = 0.8)
dev.off()


sessionInfo()


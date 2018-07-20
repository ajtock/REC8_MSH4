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

####################################################################
# DNA methylation proportion in 10-kb                              #
####################################################################

GSM980986_WT_rep2_CG <- read.table(file = paste0(bedDir, "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"))
GSM980986_WT_rep2_CHG <- read.table(file = paste0(bedDir, "GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed"))
GSM980986_WT_rep2_CHH <- read.table(file = paste0(bedDir, "GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed"))

cum_meth_dat <- NULL
for(i in 1:5) {
  # define windows2 as GRanges object
  windows2 <- seq(1, chrLens[i], by = 10000)
  windows2 <- c(windows2, chrLens[i])
  cum_windows2 <- windows2 + sumchr[i]
  windows2_iranges <- IRanges(start = windows2, width = 10000)
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
  colnames(dat2) <- c("cum_windows", "wtCG_win_vals", "wtCHG_win_vals", "wtCHH_win_vals", "allcntxts_win_vals")
  cum_meth_dat <- rbind(cum_meth_dat, dat2)
}
write.table(cum_meth_dat, file = paste0(outDir, "wtmeth_genome_10kb.txt"))


sessionInfo()


#############################################################################################
# Calculate log2 ratio of ChIP divided by input coverage values in 100-kb windows           #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/"
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
inNames <- c("SPO11-oligo_RPI1_nakedDNA_R1")
ChIPnames <- c("WT_SPO11_oligos_RPI1")
inputnames <- c("WT_nakedDNA_T50R1")
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


sessionInfo()

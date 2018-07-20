## Convert peak coordinate .txt files to .gff files for visualisation in IGV

library(doParallel)
registerDoParallel(cores=4)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

in.dir <- "/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/no_trim/analysis_01/"

peakTabList <- list("all.fft.peaks.arms.S54.txt", "all.fft.peaks.arms.S60.txt", "all.fft.peaks.arms.S70.txt", "all.fft.peaks.arms.SH975.txt")
peakGFFList <- list("all.fft.peaks.arms.S54.gff", "all.fft.peaks.arms.S60.gff", "all.fft.peaks.arms.S70.gff", "all.fft.peaks.arms.SH975.gff")

foreach(i = seq_along(peakTabList), .combine = 'c') %dopar% {
  print(i)
  peakTab <- read.table(file = paste0(in.dir, peakTabList[[i]]))
  peakGFF <- cbind(peakTab[,1], rep(".", length(peakTab[,1])), rep("nuc_peak", length(peakTab[,1])), peakTab[,2], peakTab[,3], peakTab[,8], rep("+", length(peakTab[,1])), rep(".", length(peakTab[,1])), rep(".", length(peakTab[,1])))
  #colnames(peakGFF) <- c("chr", "source", "feature", "start", "end", "REC8.fft", "strand", "frame", "rank")
  write.table(peakGFF, file = paste0(in.dir, peakGFFList[[i]]), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}


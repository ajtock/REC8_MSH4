library(GenomicRanges)

inDir <- "/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/"
outDir <- "/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/motif_analysis/"

all.fft.peaks.arms.SH99 <- read.table(file = paste0(inDir, "all.fft.peaks.arms.SH99.txt"))
all.fft.peaks.peri.SH99 <- read.table(file = paste0(inDir, "all.fft.peaks.peri.SH99.txt")) 

###Create GRanges objects to merge overlapping peaks
#score_h > 0.99
armPeaksSH99GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.arms.SH99$space), ranges = IRanges(start = all.fft.peaks.arms.SH99$start, end = all.fft.peaks.arms.SH99$end), strand = "+")
periPeaksSH99GR <- GRanges(seqnames = paste0("Chr", all.fft.peaks.peri.SH99$space), ranges = IRanges(start = all.fft.peaks.peri.SH99$start, end = all.fft.peaks.peri.SH99$end), strand = "+")
print("Unmerged arm peaks with score_h > 0.99")
print(length(armPeaksSH99GR))
#[1] 39011
print("Unmerged peri peaks with score_h > 0.99")
print(length(periPeaksSH99GR))
#[1] 22428
print("Unmerged peaks with score_h > 0.99")
print(length(armPeaksSH99GR)+length(periPeaksSH99GR))
#[1] 61439

armPeaksSH99GRmerge <- reduce(armPeaksSH99GR)
#save(armPeaksSH99GRmerge, file = paste0(inDir, "armPeaksSH99GRmerge.RData"))

periPeaksSH99GRmerge <- reduce(periPeaksSH99GR)
#save(periPeaksSH99GRmerge, file = paste0(inDir, "periPeaksSH99GRmerge.RData"))

print("Merged arm peaks with score_h > 0.99")
print(length(armPeaksSH99GRmerge))
#[1] 37529
print("Merged peri peaks with score_h > 0.99")
print(length(periPeaksSH99GRmerge))
#[1] 20205
print("Merged peaks with score_h > 0.99")
print(length(armPeaksSH99GRmerge)+length(periPeaksSH99GRmerge))
#[1] 57734

armdf <- data.frame(chr = sub("Chr", "", seqnames(armPeaksSH99GR)), start = start(armPeaksSH99GR)-1, end = end(armPeaksSH99GR),
                    midpoint = round((start(armPeaksSH99GR)-1)+((end(armPeaksSH99GR)-(start(armPeaksSH99GR)-1))/2))
                   )
peridf <- data.frame(chr = sub("Chr", "", seqnames(periPeaksSH99GR)), start = start(periPeaksSH99GR)-1, end = end(periPeaksSH99GR),
                     midpoint = round((start(periPeaksSH99GR)-1)+((end(periPeaksSH99GR)-(start(periPeaksSH99GR)-1))/2))
                    )
write.table(armdf, file = paste0(outDir, "armPeaksSH99_0based.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(peridf, file = paste0(outDir, "periPeaksSH99_0based.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE)

armdf500bp <- cbind(armdf[,1], armdf[,4]-250, armdf[,4]+250)
peridf500bp <- cbind(peridf[,1], peridf[,4]-250, peridf[,4]+250)
write.table(armdf500bp, file = paste0(outDir, "armPeaksSH99_0based_500bp.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(peridf500bp, file = paste0(outDir, "periPeaksSH99_0based_500bp.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE)

armdf <- data.frame(chr = sub("Chr", "", seqnames(armPeaksSH99GRmerge)), start = start(armPeaksSH99GRmerge)-1, end = end(armPeaksSH99GRmerge),
                    midpoint = round((start(armPeaksSH99GRmerge)-1)+((end(armPeaksSH99GRmerge)-(start(armPeaksSH99GRmerge)-1))/2))
                   )
peridf <- data.frame(chr = sub("Chr", "", seqnames(periPeaksSH99GRmerge)), start = start(periPeaksSH99GRmerge)-1, end = end(periPeaksSH99GRmerge),
                     midpoint = round((start(periPeaksSH99GRmerge)-1)+((end(periPeaksSH99GRmerge)-(start(periPeaksSH99GRmerge)-1))/2))
                    )
write.table(armdf, file = paste0(outDir, "armPeaksSH99merge_0based.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(peridf, file = paste0(outDir, "periPeaksSH99merge_0based.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE)

armdf500bp <- cbind(armdf[,1], armdf[,4]-250, armdf[,4]+250)
peridf500bp <- cbind(peridf[,1], peridf[,4]-250, peridf[,4]+250)
write.table(armdf500bp, file = paste0(outDir, "armPeaksSH99merge_0based_500bp.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(peridf500bp, file = paste0(outDir, "periPeaksSH99merge_0based_500bp.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE)



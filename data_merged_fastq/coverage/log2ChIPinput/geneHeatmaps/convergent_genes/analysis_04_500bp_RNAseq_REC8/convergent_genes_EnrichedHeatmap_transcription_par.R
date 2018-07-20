##################################################################################################################################################################
# Plot heatmaps of normalised coverage levels at convergent gene pairs, ordered according to the difference between paired genes in RNA-Seq Rep1 coverage levels #
##################################################################################################################################################################

library(parallel)
library(EnrichedHeatmap)
library(genomation)
library(circlize)
library(RColorBrewer)
sessionInfo()
inDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/geneHeatmaps/convergent_genes/definition/"
plotDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/geneHeatmaps/convergent_genes/analysis_04_500bp_RNAseq_REC8/plots/"
matDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/geneHeatmaps/convergent_genes/analysis_04_500bp_RNAseq_REC8/matrices/"
rowOrderDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/geneHeatmaps/convergent_genes/analysis_04_500bp_RNAseq_REC8/row_order/"

# Read in convergent gene pairs as GRanges objects
load(paste0(inDir, "allConvergentGR_500bp_RNAseq.RData"))
genePairs <- allConvergentGR
seqlevels(genePairs) <- sub("Chr", "", seqlevels(genePairs))
print("******Convergent gene pairs******")
print(genePairs)
# Sort convergent gene pairs by increasing transcript overlap
genePairsSD <- sort(genePairs, by = ~ Δcov, decreasing = T)
genePairsMidpointSD <- GRanges(seqnames = seqnames(genePairsSD), strand = "*",
                               ranges = IRanges(start = genePairsSD$midpoint, end = genePairsSD$midpoint),
                               gene_model_plus = genePairsSD$gene_model_plus,
                               gene_model_minus = genePairsSD$gene_model_minus,
                               TSS_plus = genePairsSD$TSS_plus,
                               TTS_plus = genePairsSD$TTS_plus,
                               TSS_minus = genePairsSD$TSS_minus,
                               TTS_minus = genePairsSD$TTS_minus,
                               TTS_minus_less_TTS_plus = genePairsSD$TTS_minus_less_TTS_plus,
                               TTS_minus_dist_TTS_plus = genePairsSD$TTS_minus_dist_TTS_plus)
rowOrderCovDiff <- list(1:length(genePairsMidpointSD))

REC8 <- system("ls /projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/REC8_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed", intern = T)
WT_RNA_Seq <- system("ls /projects/ajt200/BAM_masters/RNAseq/wt/multi_unique/bam_all/coverage/wt_RNAseq_ATCACG_norm_allchrs_coverage_coord_tab.bed", intern = T)

# Create lists/vectors of path and library names
libPaths <- list(REC8, WT_RNA_Seq)
libNames <- c("REC8", "WT_RNA_Seq")
orderNames <- c("REC8", "WT_RNA_Seq")

# Read in coverage files as GRanges objects and assign to library names
grTmp <- mclapply(seq_along(libPaths), function(x) {
  readGeneric(libPaths[[x]], meta.col = list(coverage = 4))
}, mc.cores = 2, mc.preschedule = F)
for(i in 1:length(grTmp)) {
  assign(paste0(libNames[i]), grTmp[[i]])
}

#####REC8
seqlevels(REC8) <- sub("Chr", "", seqlevels(REC8))
print("******REC8_norm_cov******")
print(REC8)

#####WT_RNA_Seq
seqlevels(WT_RNA_Seq) <- sub("Chr", "", seqlevels(WT_RNA_Seq))
print("******WT_RNA_Seq_norm_cov******")
print(WT_RNA_Seq)

# Create GRangesList object containing per base coverage for each library
grl <- GRangesList("REC8" = REC8, "WT_RNA_Seq" = WT_RNA_Seq)
orderGrl <- GRangesList("REC8" = REC8, "WT_RNA_Seq" = WT_RNA_Seq)

w <- 5
trim <- c(0, 0.01)

# Function to create coverage matrices and heatmaps for convergent gene pairs,
## with rows ordered by the extent of transcript overlap
genePairsHeatmapSD <- function(signal, target, x) {
  set.seed(1739)
  mat1_trim0.01_smoothed <- normalizeToMatrix(signal, target, value_column = "coverage",
                                              extend = 500, mean_mode = "w0", w = w, trim = trim,
                                              empty_value = 0, smooth = TRUE)
  print(mat1_trim0.01_smoothed)
  print(class(mat1_trim0.01_smoothed))
  save(mat1_trim0.01_smoothed, file = outFile[[x]])

  rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")
  #col_fun1 <- colorRamp2(quantile(mat1_trim0.01_smoothed, c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)), rich8to6equal)
  col_fun1 <- colorRamp2(quantile(mat1_trim0.01_smoothed, c(0.1, 0.95)), c("white", "black"))
  pdf(plotFiles[[x]])
  ht1 <- print(EnrichedHeatmap(mat1_trim0.01_smoothed, row_order = rowOrderCovDiff[[1]],
                               col = col_fun1, name = "Coverage",
                               column_title = c(paste0(libNames[x], " at convergent gene pairs ordered by transcriptional difference")),
                               axis_name = c("-500 bp", "0", "+500 bp"), axis_name_rot = 90, border = FALSE))
  dev.off()
  #ht1 <- draw(ht1)
  #row_order_ht1 <- row_order(ht1)
  #save(row_order_ht1, file = rowOrderFilesConstr[[x]])
}

# Run genePairsHeatmap() on all libraries
outFile <- lapply(seq_along(libNames), function(x) {
             paste0(matDir, libNames[x],
                    "_by_convergent_gene_Δtranscription_cov_mat1_trim0.00_0.01_w5_smoothed_flank_0.1to0.95.RData")
           })
rowOrderFilesConstr <- lapply(seq_along(libNames), function(x) {
                         paste0(rowOrderDir, libNames[x],
                                "_convergent_genes_cov_mat1_trim0.00_0.01_w5_smoothed_flank0.1to0.95.RData")
                       })
plotFiles <- lapply(seq_along(libNames), function(x) {
               paste0(plotDir, libNames[x],
                      "_by_convergent_gene_Δtranscription_cov_mat1_trim0.00_0.01_w5_smoothed0.1to0.95.pdf")
             })

mclapply(seq_along(grl), function(x) {
  genePairsHeatmapSD(grl[[x]], genePairsMidpointSD, x)
}, mc.cores = 2, mc.preschedule = F)

                                                
sessionInfo()


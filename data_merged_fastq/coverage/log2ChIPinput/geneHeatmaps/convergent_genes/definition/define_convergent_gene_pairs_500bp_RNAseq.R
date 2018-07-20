###########################################################################################################
# Define convergent gene pairs based on methods described in Busslinger et al. (2017) Nature 544, 503–507 #
###########################################################################################################

library(segmentSeq)
library(parallel)
library(doParallel)
library(EnrichedHeatmap)
library(genomation)
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/geneHeatmaps/convergent_genes/definition/"
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

#Read in genes as GRanges objects
representative_genes_uniq <- system("ls /projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", intern=T)
genes <- readGeneric(representative_genes_uniq, header=TRUE, strand=4, meta.col=list(gene_model=5))
print("******genes******")
print(genes)

#Create lists of path and library names
WT_RNA_Seq <- system("ls /projects/ajt200/BAM_masters/RNAseq/wt/multi_unique/bam_all/coverage/wt_RNAseq_ATCACG_norm_allchrs_coverage_coord_tab.bed", intern = T)
lib.paths <- list(WT_RNA_Seq)
lib.names <- c("WT_RNA_Seq")
ChIP.names <- c("WT_RNA_Seq")

#Read in coverage files as GRanges objects and assign to library names
gr.tmp <- mclapply(seq_along(lib.paths), function(x) {
  readGeneric(lib.paths[[x]], meta.col=list(coverage=4))
}, mc.cores=1, mc.preschedule=F)
for(i in 1:length(gr.tmp)) {
    assign(paste0(lib.names[i]), gr.tmp[[i]])
}

#####WT_RNA_Seq
seqlevels(WT_RNA_Seq) <- sub("Chr", "", seqlevels(WT_RNA_Seq))
print("******WT_RNA_Seq_norm_cov******")
print(WT_RNA_Seq)

#Create list of GRanges objects containing per base coverage for each library
grl <- GRangesList("WT_RNA_Seq"=WT_RNA_Seq)
ChIP.grl <- GRangesList("WT_RNA_Seq"=WT_RNA_Seq)

w <- 20
set.seed(123)
mat1_smoothed <- normalizeToMatrix(ChIP.grl[[1]], genes, value_colum = "coverage",
                                   extend = 0, mean_mode = "w0", w = w,
                                   empty_value = 0, smooth = TRUE,
                                   include_target = TRUE, target_ratio = 1)
geneCovMat <- mat1_smoothed[1:27204,1:20]
geneCovSum <- rowSums(geneCovMat)
values(genes) <- cbind(values(genes), data.frame(geneCovSum))

genes <- as.data.frame(genes)
genes <- genes[,-4]
#genes <- read.table(file = "/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", header = T)
genes <- cbind(chr = paste0("Chr", genes[,1]), genes[,-1])
genesPlus <- genes[genes$strand == "+",]
genesMinus <- genes[genes$strand == "-",]

genesPlusGR <- NULL
genesMinusGR <- NULL
for(i in 1:5) {
  chrGenesPlus <- genesPlus[genesPlus$chr == chrs[i],]
  chrGenesPlusIR <- IRanges(start = chrGenesPlus$start, end = chrGenesPlus$end)
  chrGenesPlusGR <- GRanges(seqnames = chrGenesPlus$chr, strand = chrGenesPlus$strand, ranges = chrGenesPlusIR,
                            gene_model = chrGenesPlus$gene_model, geneCovSum = chrGenesPlus$geneCovSum)
  genesPlusGR <- c(genesPlusGR, chrGenesPlusGR)
  chrGenesMinus <- genesMinus[genesMinus$chr == chrs[i],]
  chrGenesMinusIR <- IRanges(start = chrGenesMinus$start, end = chrGenesMinus$end)
  chrGenesMinusGR <- GRanges(seqnames = chrGenesMinus$chr, strand = chrGenesMinus$strand, ranges = chrGenesMinusIR,
                             gene_model = chrGenesMinus$gene_model, geneCovSum = chrGenesMinus$geneCovSum) 
  genesMinusGR <- c(genesMinusGR, chrGenesMinusGR)
}
  
genesPlusExtGR <- NULL
genesMinusExtGR <- NULL
for(i in 1:5) {
  chrGenesPlus <- genesPlus[genesPlus$chr == chrs[i],]
  # sense gene "convergent region" = TTS +/- 500 bp
  chrGenesPlusExtIR <- IRanges(start = pmax(1, chrGenesPlus$end-500), end = pmin(chrGenesPlus$end+500, chrLens[i]))
  chrGenesPlusExtGR <- GRanges(seqnames = chrGenesPlus$chr, strand = chrGenesPlus$strand, ranges = chrGenesPlusExtIR,
                               gene_model = chrGenesPlus$gene_model, geneCovSum = chrGenesPlus$geneCovSum)
  genesPlusExtGR <- c(genesPlusExtGR, chrGenesPlusExtGR)
  chrGenesMinus <- genesMinus[genesMinus$chr == chrs[i],]
  # antisense gene "convergent region" = TTS +/- 500 bp
  chrGenesMinusExtIR <- IRanges(start = pmax(1, chrGenesMinus$start-500), end = pmin(chrGenesMinus$start+500, chrLens[i]))
  chrGenesMinusExtGR <- GRanges(seqnames = chrGenesMinus$chr, strand = chrGenesMinus$strand, ranges = chrGenesMinusExtIR,
                                gene_model = chrGenesMinus$gene_model, geneCovSum = chrGenesMinus$geneCovSum)
  genesMinusExtGR <- c(genesMinusExtGR, chrGenesMinusExtGR)
}

overlapsPlusExt <- lapply(seq_along(genesPlusExtGR), function(x) {
  findOverlaps(genesPlusExtGR[[x]], genesMinusGR[[x]], ignore.strand = TRUE, select = "all")
})

allConvergentGR <- GRanges()
for(c in 1:length(overlapsPlusExt)) {
  #print(c)
  for(h in 1:length(overlapsPlusExt[[c]])) {
    #print(h)
    if (
      end(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])) < end(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]))
      && start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])])) <= end(ranges(genesPlusExtGR[[c]][queryHits(overlapsPlusExt[[c]][h])]))
      && start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])])) >= start(ranges(genesPlusExtGR[[c]][queryHits(overlapsPlusExt[[c]][h])]))
    ) {
    #print("Convergent")
    convergentIR <- IRanges(start = pmin(start(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])), 
                                         start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]))),
                            end = end(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]))
                           )
    convergentGR <- GRanges(seqnames = seqnames(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])]), strand = "*", ranges = convergentIR,
                            gene_model_plus = genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])]$gene_model,
                            gene_model_minus = genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]$gene_model,
                            geneCovSum_plus = genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])]$geneCovSum,
                            geneCovSum_minus = genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]$geneCovSum,
                            Δcov = genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])]$geneCovSum-genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]$geneCovSum,
                            TSS_plus = start(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])),
                            TTS_plus = end(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])),
                            TTS_minus = start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])])),
                            TSS_minus = end(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])])),
                            TTS_minus_less_TTS_plus  = start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])])) - 
                                                       end(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])),
                            TTS_minus_dist_TTS_plus  = -1 * ( start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])])) -
                                                              end(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])) ),
                            midpoint = round(
                                       pmin(end(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])), 
                                            start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]))) +
                                       ( ( pmax(end(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])), 
                                                start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]))) - 
                                           pmin(end(ranges(genesPlusGR[[c]][queryHits(overlapsPlusExt[[c]][h])])), 
                                                start(ranges(genesMinusGR[[c]][subjectHits(overlapsPlusExt[[c]][h])]))) ) / 2 ) )
                           )
    print(convergentGR)
    allConvergentGR <- append(allConvergentGR, convergentGR)
    }
  }
}
save(allConvergentGR, file = paste0(outDir, "allConvergentGR_500bp_RNAseq.RData"))


print(sessionInfo())


###########################################################################################################
# Define convergent gene pairs based on methods described in Busslinger et al. (2017) Nature 544, 503–507 #
###########################################################################################################

library(segmentSeq)
library(parallel)
library(doParallel)
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/geneProfiles/convergent_genes/"
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

genes <- read.table(file = "/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", header = T)
genes <- cbind(chr = paste0("Chr", genes[,1]), genes[,-1])
genesPlus <- genes[genes$strand == "+",]
genesMinus <- genes[genes$strand == "-",]

genesPlusGR <- NULL
genesMinusGR <- NULL
for(i in 1:5) {
  chrGenesPlus <- genesPlus[genesPlus$chr == chrs[i],]
  chrGenesPlusIR <- IRanges(start = chrGenesPlus$start, end = chrGenesPlus$end)
  chrGenesPlusGR <- GRanges(seqnames = chrGenesPlus$chr, strand = chrGenesPlus$strand, ranges = chrGenesPlusIR, gene_model = chrGenesPlus$gene_model)
  genesPlusGR <- c(genesPlusGR, chrGenesPlusGR)
  chrGenesMinus <- genesMinus[genesMinus$chr == chrs[i],]
  chrGenesMinusIR <- IRanges(start = chrGenesMinus$start, end = chrGenesMinus$end)
  chrGenesMinusGR <- GRanges(seqnames = chrGenesMinus$chr, strand = chrGenesMinus$strand, ranges = chrGenesMinusIR, gene_model = chrGenesMinus$gene_model) 
  genesMinusGR <- c(genesMinusGR, chrGenesMinusGR)
}
  
genesPlusExtGR <- NULL
genesMinusExtGR <- NULL
for(i in 1:5) {
  chrGenesPlus <- genesPlus[genesPlus$chr == chrs[i],]
  # sense gene "convergent region" = TTS +/- 5000 bp
  chrGenesPlusExtIR <- IRanges(start = pmax(1, chrGenesPlus$end-5000), end = pmin(chrGenesPlus$end+5000, chrLens[i]))
  chrGenesPlusExtGR <- GRanges(seqnames = chrGenesPlus$chr, strand = chrGenesPlus$strand, ranges = chrGenesPlusExtIR, gene_model = chrGenesPlus$gene_model)
  genesPlusExtGR <- c(genesPlusExtGR, chrGenesPlusExtGR)
  chrGenesMinus <- genesMinus[genesMinus$chr == chrs[i],]
  # antisense gene "convergent region" = TTS +/- 5000 bp
  chrGenesMinusExtIR <- IRanges(start = pmax(1, chrGenesMinus$start-5000), end = pmin(chrGenesMinus$start+5000, chrLens[i]))
  chrGenesMinusExtGR <- GRanges(seqnames = chrGenesMinus$chr, strand = chrGenesMinus$strand, ranges = chrGenesMinusExtIR, gene_model = chrGenesMinus$gene_model)
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
save(allConvergentGR, file = paste0(outDir, "allConvergentGR_5kb.RData"))


print(sessionInfo())


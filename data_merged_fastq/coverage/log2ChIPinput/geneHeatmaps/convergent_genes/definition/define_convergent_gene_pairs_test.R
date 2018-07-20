###########################################################################################################
# Define convergent gene pairs based on methods described in Busslinger et al. (2017) Nature 544, 503–507 #
###########################################################################################################

library(segmentSeq)
library(parallel)
library(doParallel)
#library(genomation)
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/geneProfiles/convergent_genes/"
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

#representative_genes_uniq <- system("ls /projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", intern=T)
#genes2 <- readGeneric(representative_genes_uniq, header=TRUE, strand=4, meta.col=list(gene_model=5))
#for(i in 1:5) {  
# seqlevels(genes2) <- sub(i, chrs[i], seqlevels(genes2))
#}

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
  # sense gene "convergent region" = TTS +/- 50 bp
  chrGenesPlusExtIR <- IRanges(start = pmax(1, chrGenesPlus$end-50), end = pmin(chrGenesPlus$end+50, chrLens[i]))
  chrGenesPlusExtGR <- GRanges(seqnames = chrGenesPlus$chr, strand = chrGenesPlus$strand, ranges = chrGenesPlusExtIR, gene_model = chrGenesPlus$gene_model)
  genesPlusExtGR <- c(genesPlusExtGR, chrGenesPlusExtGR)
  chrGenesMinus <- genesMinus[genesMinus$chr == chrs[i],]
  # antisense gene "convergent region" = TTS +/- 50 bp
  chrGenesMinusExtIR <- IRanges(start = pmax(1, chrGenesMinus$start-50), end = pmin(chrGenesMinus$start+50, chrLens[i]))
  chrGenesMinusExtGR <- GRanges(seqnames = chrGenesMinus$chr, strand = chrGenesMinus$strand, ranges = chrGenesMinusExtIR, gene_model = chrGenesMinus$gene_model)
  genesMinusExtGR <- c(genesMinusExtGR, chrGenesMinusExtGR)
}

#if(require("parallel"))
#{
#    numCores <- min(8, detectCores())
#    cl <- makeCluster(numCores)
#} else {
#    cl <- NULL
}

#overlapsPlusExt <- lapply(seq_along(genesPlusExtGR), function(x) {
#  getOverlaps(genesPlusExtGR[[x]], genesMinusGR[[x]], overlapType = "overlapping", whichOverlaps = TRUE, ignoreStrand = TRUE, cl = cl)
#})

#is.integer0 <- function(x)
#{
#  is.integer(x) && length(x) == 0L
#}

overlapsPlusExt <- lapply(seq_along(genesPlusExtGR), function(x) {
  findOverlaps(genesPlusExtGR[[x]], genesMinusGR[[x]], ignore.strand = TRUE, select = "all")
})

registerDoParallel(cores=48)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())


#for(c in 1:length(overlapsPlusExt)) {
#  for(gp in 1:length(overlapsPlusExt[[c]])) {
#    for(gm in 1:length(overlapsPlusExt[[c]][[gp]])) {
#      if (
#        end(ranges(genesPlusGR[[c]][gp])) < end(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))
#        && start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])) <= end(ranges(genesPlusExtGR[[c]][gp]))
#        && start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])) >= start(ranges(genesPlusExtGR[[c]][gp]))
#      ) {
#      print("Convergent")
#      }
#    }
#  }
#}

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






#for(c in 1:length(overlapsPlusExt)) {
#  for(gp in 1:length(overlapsPlusExt[[c]])) {
#    for(gm in 1:length(overlapsPlusExt[[c]][[gp]])) {
#      if (
#        end(ranges(genesPlusGR[[c]][gp])) < end(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))
#        && start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])) <= end(ranges(genesPlusExtGR[[c]][gp]))
#        && start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])) >= start(ranges(genesPlusExtGR[[c]][gp]))
#      ) {
#      print("Convergent")
#      }
#    }
#  }
#}

allConvergentGR <- GRanges()
foreach(c = 1:length(overlapsPlusExt)) %:%
  foreach(gp = 1:length(overlapsPlusExt[[c]])) %:%
    foreach(gm = 1:length(overlapsPlusExt[[c]][[gp]])) %dopar% {
      if (
        end(ranges(genesPlusGR[[c]][gp])) < end(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))
        && start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])) <= end(ranges(genesPlusExtGR[[c]][gp]))
        && start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])) >= start(ranges(genesPlusExtGR[[c]][gp])) 
      ) {
      print("Convergent")
      convergentIR <- IRanges(start = pmin(start(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))),
                              end = end(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))
                             )
      convergentGR <- GRanges(seqnames = seqnames(genesPlusGR[[c]][gp]), strand = "*", ranges = convergentIR,
                              gene_model_plus = genesPlusGR[[c]][gp]$gene_model,
                              gene_model_minus = genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]$gene_model,
                              TSS_plus = start(ranges(genesPlusGR[[c]][gp])),
                              TTS_plus = end(ranges(genesPlusGR[[c]][gp])),
                              TTS_minus = start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])),
                              TSS_minus = end(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])),
                              TTS_minus_less_TTS_plus  = start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))-end(ranges(genesPlusGR[[c]][gp])),
                              midpoint = round(
                                         pmin(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))) +
                                         ( ( pmax(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))) -
                                             pmin(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))) ) / 2 ) )
                             )
      print(convergentGR)
      allConvergentGR <- append(allConvergentGR, convergentGR)
      }
    }



allConvergentGR <- GRanges()
for(c in 5:5) {
  for(gp in 1:length(overlapsPlusExt[[c]])) {
    for(gm in 1:length(overlapsPlusExt[[c]][[gp]])) {
      if (
        end(ranges(genesPlusGR[[c]][gp])) < end(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))
        && start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])) <= end(ranges(genesPlusExtGR[[c]][gp]))
        && start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])) >= start(ranges(genesPlusExtGR[[c]][gp]))
      ) {
      print("Convergent")
      convergentIR <- IRanges(start = pmin(start(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))),
                              end = end(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))
                             )
      convergentGR <- GRanges(seqnames = seqnames(genesPlusGR[[c]][gp]), strand = "*", ranges = convergentIR,
                              gene_model_plus = genesPlusGR[[c]][gp]$gene_model,
                              gene_model_minus = genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]$gene_model,
                              TSS_plus = start(ranges(genesPlusGR[[c]][gp])),
                              TTS_plus = end(ranges(genesPlusGR[[c]][gp])),
                              TTS_minus = start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])),
                              TSS_minus = end(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])),
                              TTS_minus_less_TTS_plus  = start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))-end(ranges(genesPlusGR[[c]][gp])),
                              midpoint = round(
                                         pmin(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))) +
                                         ( ( pmax(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))) -
                                             pmin(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]]))) ) / 2 ) )
                             )
      print(convergentGR)
      allConvergentGR <- append(allConvergentGR, convergentGR)
      }
    }
  }
}





      print(genesPlusGR[[c]][gp])
      print("**********AND*********")
      print(genesMinusGR[[c]][overlapsPlusExt[[c]][[gp]][[gm]]])


      
      #if ( end(ranges(genesPlusGR[[c]][gp])) < end(ranges(genesMinusGR[[c]][gm]))
      #  && start(ranges(genesMinusGR[[c]][gp])) <= end(ranges(genesPlusExtGR[[c]][gm]))
      #  && start(ranges(genesMinusGR[[c]][gp])) >= start(ranges(genesPlusExtGR[[c]][gm])) ) {
      #  print("Convergent")
      #}
    }
  }
}

        # TSS of antisense gene must be downstream of TTS of sense gene
        end(ranges(genesPlusGR[[c]][gp])) < end(ranges(genesMinusGR[[c]][gm]))
        # gene on antisense strand must end within the "convergent region" ("start" = TTS of antisense gene)
        && start(ranges(genesMinusGR[[c]][gp])) <= end(ranges(genesPlusExtGR[[c]][gm]))
        && start(ranges(genesMinusGR[[c]][gp])) >= start(ranges(genesPlusExtGR[[c]][gm])) 
      ) {
        print("Convergent")
      }
    }
  }
}
        convergentIR <- IRanges(start = pmin(start(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][gm]))),
                                end = end(ranges(genesMinusGR[[c]][gm]))
                               )
        convergentGR <- GRanges(seqnames = seqnames(genesPlusGR[[c]][gp]), strand = "*", ranges = convergentIR,
                                gene_model_plus = 
      }
    }
}

allConvergentGR <- GRanges()
tmp <- lapply(seq_along(overlapsPlusExt), function(c) {
  mclapply(seq_along(overlapsPlusExt[[c]]), function(gp) {
    lapply(seq_along(overlapsPlusExt[[c]][[gp]]), function(gm) {
      if ( end(ranges(genesPlusGR[[c]][gp])) < end(ranges(genesMinusGR[[c]][gm]))
        && start(ranges(genesMinusGR[[c]][gp])) <= end(ranges(genesPlusExtGR[[c]][gm]))
        && start(ranges(genesMinusGR[[c]][gp])) >= start(ranges(genesPlusExtGR[[c]][gm])) ) {
          #print("Convergent")
            convergentIR <- IRanges(start = pmin(start(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][gm]))),
                                    end = end(ranges(genesMinusGR[[c]][gm]))
                                   )
            convergentGR <- GRanges(seqnames = seqnames(genesPlusGR[[c]][gp]), strand = "*", ranges = convergentIR,
                                    gene_model_plus = genesPlusGR[[c]][gp]$gene_model,
                                    gene_model_minus = genesMinusGR[[c]][gm]$gene_model,
                                    TSS_plus = start(ranges(genesPlusGR[[c]][gp])),
                                    TTS_plus = end(ranges(genesPlusGR[[c]][gp])),
                                    TTS_minus = start(ranges(genesMinusGR[[c]][gm])),
                                    TSS_minus = end(ranges(genesMinusGR[[c]][gm])),
                                    TTS_minus_less_TTS_plus  = start(ranges(genesMinusGR[[c]][gm]))-end(ranges(genesPlusGR[[c]][gp])),
                                    midpoint = round(
                                               pmin(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][gm]))) +
                                               ( ( pmax(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][gm]))) -
                                                   pmin(end(ranges(genesPlusGR[[c]][gp])), start(ranges(genesMinusGR[[c]][gm]))) ) / 2 ) )
                                   )
            print(convergentGR)
            allConvergentGR <- append(allConvergentGR, convergentGR)
      }
    })
  }, mc.cores = 48)
})

allConvergentGR <- GRanges()
tmp2 <- lapply(seq_along(tmp), function(chr) {
          chrConvergentGR <- GRanges()
          lapply(seq_along(tmp[[2]]), function(gp) {
            lapply(seq_along(tmp[[2]][[gp]]), function(gm) {
              #if ( !is.null(tmp[[2]][[gp]][[gm]]) == TRUE) {
              #print(tmp[[2]][[gp]][[gm]])
              #}
            convergentGR <- tmp[[2]][is.null(tmp[[2]][[gp]][[gm]] == FALSE)]
            chrConvergentGR <- append(chrConvergentGR, convergentGR)
            save(chrConvergentGR, file = "chrConvergentGR.RData")
            })
          })
        })


tmp5 <- lapply(seq_along(overlapsPlusExt), function(c) {
  mclapply(seq_along(overlapsPlusExt[[5]]), function(gp) {
    lapply(seq_along(overlapsPlusExt[[5]][[gp]]), function(gm) {
      if ( end(ranges(genesPlusGR[[5]][gp])) < end(ranges(genesMinusGR[[5]][gm]))
        && start(ranges(genesMinusGR[[5]][gp])) <= end(ranges(genesPlusExtGR[[5]][gm]))
        && start(ranges(genesMinusGR[[5]][gp])) >= start(ranges(genesPlusExtGR[[5]][gm])) ) {
          #print("Convergent")
            convergentIR <- IRanges(start = pmin(start(ranges(genesPlusGR[[5]][gp])), start(ranges(genesMinusGR[[5]][gm]))),
                                    end = end(ranges(genesMinusGR[[5]][gm]))
                                   )
            convergentGR <- GRanges(seqnames = seqnames(genesPlusGR[[5]][gp]), strand = "*", ranges = convergentIR,
                                    gene_model_plus = genesPlusGR[[5]][gp]$gene_model,
                                    gene_model_minus = genesMinusGR[[5]][gm]$gene_model,
                                    TSS_plus = start(ranges(genesPlusGR[[5]][gp])),
                                    TTS_plus = end(ranges(genesPlusGR[[5]][gp])),
                                    TTS_minus = start(ranges(genesMinusGR[[5]][gm])),
                                    TSS_minus = end(ranges(genesMinusGR[[5]][gm])),
                                    TTS_minus_less_TTS_plus  = start(ranges(genesMinusGR[[5]][gm]))-end(ranges(genesPlusGR[[5]][gp])),
                                    midpoint = round(
                                               pmin(end(ranges(genesPlusGR[[5]][gp])), start(ranges(genesMinusGR[[5]][gm]))) +
                                               ( ( pmax(end(ranges(genesPlusGR[[5]][gp])), start(ranges(genesMinusGR[[5]][gm]))) -
                                                   pmin(end(ranges(genesPlusGR[[5]][gp])), start(ranges(genesMinusGR[[5]][gm]))) ) / 2 ) )

                                   )
            print(convergentGR)
            allConvergentGR <- append(allConvergentGR, convergentGR)
      }
    })
  }, mc.cores = 48)
})





allConvergentGR <- GRanges()
mclapply(seq_along(genesPlusGR), function(x) {
  lapply(seq_along(genesMinusGR[[x]]), function(y) {
    mclapply(seq_along(genesPlusGR[[x]]), function(z) {
      if ( end(ranges(genesPlusGR[[x]][z])) < end(ranges(genesMinusGR[[x]][y]))
        && start(ranges(genesMinusGR[[x]][y])) <= end(ranges(genesPlusExtGR[[x]][z]))
        && start(ranges(genesMinusGR[[x]][y])) >= start(ranges(genesPlusExtGR[[x]][z])) ) {
          print("Convergent")
            convergentIR <- IRanges(start = pmin(start(ranges(genesPlusGR[[x]][z])), start(ranges(genesMinusGR[[x]][y]))), end = end(ranges(genesMinusGR[[x]][y])))
            convergentGR <- GRanges(seqnames = genesPlusGR$chr, strand = "*", ranges = convergentIR, 
                                    gene_model_plus = genesPlusGR[[x]][z]$gene_model,
                                    TSS_plus = start(ranges(genesPlusGR[[x]][z])),
                                    TTS_plus = end(ranges(genesPlusGR[[x]][z])),
                                    gene_model_minus = genesMinusGR[[x]][y]$gene_model,
                                    TTS_minus = start(ranges(genesMinusGR[[x]][y])),
                                    TSS_minus = end(ranges(genesMinusGR[[x]][y])),
                                    TTS_minus_less_TTS_plus  = start(ranges(genesMinusGR[[x]][y]))-end(ranges(genesPlusGR[[x]][z])),
                                    midpoint = round( 
                                               pmin(end(ranges(genesPlusGR[[x]][z])), start(ranges(genesMinusGR[[x]][y]))) +
                                               ( ( pmax(end(ranges(genesPlusGR[[x]][z])), start(ranges(genesMinusGR[[x]][y]))) -
                                                   pmin(end(ranges(genesPlusGR[[x]][z])), start(ranges(genesMinusGR[[x]][y]))) ) / 2 ) )
                                   )
            allConvergentGR <- append(allConvergentGR, convergentGR)     
      }
    }, mc.cores = 48)
  })
}, mc.cores = 1)

#interactions <- vector(mode="list", length=N)
#for (i in seq_len(N)) {
#    gr <- ... make some GRanges object ...   
#    interactions[[i]] <- gr
#}
#interactions <- do.call(c, interactions)


allConvergentGR <- GRanges()
mclapply(seq_along(genesPlusGR), function(x) {
  lapply(seq_along(genesMinusGR[[x]]), function(y) {
    lapply(seq_along(genesPlusGR[[x]]), function(z) {
      if ( end(ranges(genesPlusGR[[x]][z])) < end(ranges(genesMinusGR[[x]][y]))
        && start(ranges(genesMinusGR[[x]][y])) <= end(ranges(genesPlusExtGR[[x]][z]))
        && start(ranges(genesMinusGR[[x]][y])) >= start(ranges(genesPlusExtGR[[x]][z])) ) {
          print("Convergent")
            convergentIR <- IRanges(start = pmin(start(ranges(genesPlusGR[[x]][z])), start(ranges(genesMinusGR[[x]][y]))), end = end(ranges(genesMinusGR[[x]][y])))
            convergentGR <- GRanges(seqnames = genesPlusGR$chr, strand = "*", ranges = convergentIR,
                                    gene_model_plus = genesPlusGR[[x]][z]$gene_model,
                                    TSS_plus = start(ranges(genesPlusGR[[x]][z])),
                                    TTS_plus = end(ranges(genesPlusGR[[x]][z])),
                                    gene_model_minus = genesMinusGR[[x]][y]$gene_model,
                                    TTS_minus = start(ranges(genesMinusGR[[x]][y])),
                                    TSS_minus = end(ranges(genesMinusGR[[x]][y])),
                                    TTS_minus_less_TTS_plus  = start(ranges(genesMinusGR[[x]][y]))-end(ranges(genesPlusGR[[x]][z])),
                                    midpoint = round(
                                               pmin(end(ranges(genesPlusGR[[x]][z])), start(ranges(genesMinusGR[[x]][y]))) +
                                               ( ( pmax(end(ranges(genesPlusGR[[x]][z])), start(ranges(genesMinusGR[[x]][y]))) -
                                                   pmin(end(ranges(genesPlusGR[[x]][z])), start(ranges(genesMinusGR[[x]][y]))) ) / 2 ) )
                                   )
            allConvergentGR <- append(allConvergentGR, convergentGR)
      }
    })
  })
})


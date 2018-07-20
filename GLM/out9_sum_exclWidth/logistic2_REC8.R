# GLM to evaluate relationship between crossovers, REC8_ChIP, SPO11-1-oligos and and nucleosomes
# Script adapted from Tom Hardcastle's

library(segmentSeq)
inDir <- "/projects/ajt200/REC8_MSH4/GLM/in/"
outDir <- "/projects/ajt200/REC8_MSH4/GLM/out9_sum_exclWidth/"
genesDir <- "/projects/ajt200/TAIR10/representative_genes/"
tesDir <- "/projects/ajt200/TAIR10/"
covDir1 <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/"
covDir2 <- "/projects/ajt200/BAM_masters/nucleosomes/coverage/"
covDir3 <- "/projects/ajt200/BAM_masters/SPO11_oligo_Ian/coverage/"
chrlens <- c(Chr1 = 30427671, Chr2 = 19698289, Chr3 = 23459830, Chr4 = 18585056, Chr5 = 26975502)
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

####################################################################################################################
# 481252 SNPs -> GRanges object 481257 length (+5 = chrs?) -> defines intervals in which a crossover can de detected #
####################################################################################################################

snps <- read.table(file = paste0(inDir, "BC.complete.tiger.txt"))
snpGR <- do.call("c", lapply(1:5, function(chr) {
    chrsnps <- snps[snps[,1] == chr,]
    GRanges(seqnames = names(chrlens)[chr], IRanges(start = c(1, chrsnps[,2]), end = c(chrsnps[,2], chrlens[chr])))
}))

##############################################################################################
# load 2,499 crossovers -> note Tom reduces to 2,386 ranges -> underestimates hot intervals? #
##############################################################################################

load(file = paste0(inDir, "wt.cos.gr.coords.RData"))
cos.gr.coords <- wt.cos.gr.coords
cos.gr.coords <- sort(cos.gr.coords)
strand(cos.gr.coords) <- "*"
cos.gr.coords <- reduce(cos.gr.coords)

##################################################################################################################
# overlap to identify SNP windows which overlap a crossover, check type = "within", overlap = 3295 SNP intervals #
# overlapType = "within" - SNP windows that are completely contained by at least one CO interval                 #
##################################################################################################################

coevents <- getOverlaps(snpGR, cos.gr.coords, overlapType = "within", whichOverlaps = FALSE)
length(which(coevents[1:length(coevents)] == TRUE))
#[1] 3295
length(which(coevents[1:length(coevents)] == FALSE))
#[1] 477962

######################################################################################
# regions? 480348 regions                                                            # 
# (477962 SNP windows that are not completely contained by at least one CO interval #
# + 2386 (all) CO intervals)                                                         #
######################################################################################

regions <- c(snpGR[!coevents], cos.gr.coords) 
regions <- sort(regions)

#################################################
# genes = 27206, tes = 29150 -> GRanges objects #
#################################################
###############################################################
# split TEs into families, define gene promoters, terminators #
###############################################################

genes <- read.table(file = paste0(genesDir, "representative_genes_uniq_fmt_strand.txt"), header = TRUE)
tes <- read.table(file = paste0(tesDir, "TAIR10_Buisine_TEs_strand_tab_ann.txt"), header = TRUE)
genes <- GRanges(seqnames = paste("Chr", genes[,1], sep = ""), IRanges(start = genes[,2], end = genes[,3]), strand = genes[,4])
#strand = c("+", "-")[as.numeric(genes$strand == "REVERSE") + 1])
genes <- genes[order(as.character(seqnames(genes)), start(genes)),]
tes <- GRanges(seqnames = tes[,1], IRanges(start = tes[,2], end = tes[,3]), strand = tes[,4], sf = tes[,7])
tes <- tes[order(as.character(seqnames(tes)), start(tes)),]
tesList <- lapply(levels(tes$sf), function(sf) {
    tes[values(tes)$sf == sf,]
})
names(tesList) <- levels(tes$sf)

makePromoters <- function(genes) {
    promoters <- genes
    end(promoters[strand(promoters) == "+"]) <- start(genes[strand(genes) == "+"])-1
    start(promoters[strand(promoters) == "+"]) <- start(genes[strand(genes) == "+"])-500
    start(promoters[strand(promoters) == "-"]) <- end(genes[strand(genes) == "-"])+1
    end(promoters[strand(promoters) == "-"]) <- end(genes[strand(genes) == "-"])+500
    start(promoters)[(start(promoters)<1)] <- 1
    promoters
}
promoters <- makePromoters(genes)

makeTerminators <- function(genes) {
    terminators <- genes
    start(terminators[strand(terminators) == "+"]) <- end(genes[strand(genes) == "+"])-1
    end(terminators[strand(terminators) == "+"]) <- end(genes[strand(genes) == "+"])+500
    start(terminators[strand(terminators) == "-"]) <- start(genes[strand(genes) == "-"])-500
    end(terminators[strand(terminators) == "-"]) <- start(genes[strand(genes) == "-"])+1
    start(terminators)[(start(terminators)<1)] <- 1
    terminators
}
terminators <- makeTerminators(genes)

########################################################################
# overlap gene, TE, promoter annotations with regions (~SNP intervals) #
########################################################################

geneOverlap <- getOverlaps(regions, genes, whichOverlaps = FALSE)
promoterOverlap <- getOverlaps(regions, promoters, whichOverlaps = FALSE)
terminatorOverlap <- getOverlaps(regions, terminators, whichOverlaps = FALSE)
tesOverlap <- do.call("cbind", lapply(tesList, function(x) {
    getOverlaps(regions, x, whichOverlaps = FALSE)
}))

##############################################################################################
# load REC8_ChIP data, split according to region, calculate mean, max, min, sum in each region #
##############################################################################################

splitREC8_ChIP <- list()
for(REC8_ChIPChr in dir(path = covDir1, pattern = "REC8_ChIP_norm_Chr[1-5]_coverage.txt")) {
    REC8_ChIP <- read.table(file = paste0(covDir1, REC8_ChIPChr))[,1]
    chreg <- regions[seqnames(regions) == seqlevels(regions)[as.numeric(gsub("REC8_ChIP_norm_Chr([1-5])_coverage.txt", "\\1", REC8_ChIPChr))]]
    splitREC8_ChIP <- c(splitREC8_ChIP, split(REC8_ChIP, c(rep(1:length(chreg), width(chreg)-1), length(chreg))))
}
meanREC8_ChIP <- sapply(splitREC8_ChIP, mean)
maxREC8_ChIP <- sapply(splitREC8_ChIP, max)
minREC8_ChIP <- sapply(splitREC8_ChIP, min)
sumREC8_ChIP <- sapply(splitREC8_ChIP, sum)

# Do the same for MSH4_ChIP data
splitMSH4_ChIP <- list()
for(MSH4_ChIPChr in dir(path = covDir1, pattern = "MSH4_ChIP_norm_Chr[1-5]_coverage.txt")) {
    MSH4_ChIP <- read.table(file = paste0(covDir1, MSH4_ChIPChr))[,1]
    chreg <- regions[seqnames(regions) == seqlevels(regions)[as.numeric(gsub("MSH4_ChIP_norm_Chr([1-5])_coverage.txt", "\\1", MSH4_ChIPChr))]]
    splitMSH4_ChIP <- c(splitMSH4_ChIP, split(MSH4_ChIP, c(rep(1:length(chreg), width(chreg)-1), length(chreg))))
}
meanMSH4_ChIP <- sapply(splitMSH4_ChIP, mean)
maxMSH4_ChIP <- sapply(splitMSH4_ChIP, max)
minMSH4_ChIP <- sapply(splitMSH4_ChIP, min)
sumMSH4_ChIP <- sapply(splitMSH4_ChIP, sum)

# Do the same for nucleosomes data
splitnuc <- list()
for(i in 1:5) {
    nuc <- read.table(file = paste0(covDir2, "WT_nuc_chr", i, "_unique_both_coverage.txt"))[,1]
    chreg <- regions[seqnames(regions) == chrs[i]]
    splitnuc <- c(splitnuc, split(nuc, c(rep(1:length(chreg), width(chreg)-1), length(chreg))))
}
meannuc <- sapply(splitnuc, mean)
maxnuc <- sapply(splitnuc, max)
minnuc <- sapply(splitnuc, min)
sumnuc <- sapply(splitnuc, sum)

# Do the same for SPO11-oligos data
splitspo <- list()
for(spochr in dir(path = covDir3, pattern = "SPO11_oligo_chr[1-5]_cov.norm.both.txt")) {
    spo11 <- read.table(file = paste0(covDir3, spochr))[,1]
    chreg <- regions[seqnames(regions) == seqlevels(regions)[as.numeric(gsub("SPO11_oligo_chr([1-5])_cov.norm.both.txt", "\\1", spochr))]]
    splitspo <- c(splitspo, split(spo11, c(rep(1:length(chreg), width(chreg)-1), length(chreg))))
}
meanspo <- sapply(splitspo,mean)
maxspo <- sapply(splitspo,max)
minspo <- sapply(splitspo,min)
sumspo <- sapply(splitspo,sum)

################################################################################################################################
# define centromeres and banding regions - defines centromeres as 2Mb around CEN, bands are [ajt: 2 Mb] intervals in arms #
################################################################################################################################

centromeres <- read.delim(file = paste0(inDir, "centromere.txt"), header=FALSE)
centGR <- GRanges(seqnames = paste("Chr", centromeres[,1], sep = ""), IRanges(start = centromeres[,2]-2e6, end = centromeres[,3]+2e6))
centromere <- getOverlaps(regions, centGR, whichOverlaps = FALSE)
pericentGR <- sort(c(GRanges(seqnames(centGR), IRanges(start = pmax(1, start(centGR)-1-2e6), end = start(centGR)-1)), GRanges(seqnames(centGR), IRanges(start = end(centGR)+1, end = end(centGR)+1+2e6))))
pericent <- getOverlaps(regions, pericentGR, whichOverlaps = FALSE)
banding <- do.call("c", lapply(seqlevels(regions), function(chr){
    chrlen <- max(end(regions[seqnames(regions) == chr]))
    cent <- centGR[which(seqnames(centGR) == chr)]
    divL <- round((start(cent)-1)/2e6)
    # ajt: added brackets around end(cent)+1
    divR <- round((chrlen-(end(cent)+1))/2e6)
    # ajt: added brackets around end(cent)+1
    ends <- c(seq.int(divL)*(start(cent)-1)/divL, end(cent), seq.int(divR)*((chrlen-(end(cent)+1))/divR)+end(cent))
    GRanges(seqnames = chr, IRanges(start = c(1, ends[-length(ends)]+1), end = ends))
}))
bandRegions <- sapply(getOverlaps(regions, banding, whichOverlaps = TRUE), function(x) x[1])

##########################################################################################
# create data object for model -> contains CO overlap, SPO11-1 data, annotation overlaps #
##########################################################################################

dat <- cbind.data.frame(co = regions %in% cos.gr.coords,
                        meanREC8_ChIP = meanREC8_ChIP, maxREC8_ChIP = maxREC8_ChIP, minREC8_ChIP = minREC8_ChIP, sumREC8_ChIP = sumREC8_ChIP, 
                        meanMSH4_ChIP = meanMSH4_ChIP, maxMSH4_ChIP = maxMSH4_ChIP, minMSH4_ChIP = minMSH4_ChIP, sumMSH4_ChIP = sumMSH4_ChIP, 
                        meanspo = meanspo*1e9, maxspo = maxspo*1e9, minspo = minspo*1e9, sumspo = sumspo*1e9,
                        meannuc = meannuc, maxnuc = maxnuc, minnuc = minnuc, sumnuc = sumnuc,                        
                        gene = geneOverlap, promoter = promoterOverlap, terminator = terminatorOverlap, tesOverlap,
                        band = as.factor(bandRegions),
                        width = width(regions),
                        centromere = centromere,
                        pericent = pericent)
colnames(dat) <- gsub("/", "_", colnames(dat))
colnames(dat) <- gsub("-", ".", colnames(dat))
dat$centromere[dat$centromere] <- as.character(seqnames(regions[centromere]))
dat$centromere[dat$centromere == FALSE] <- "Chr0"
dat$pericent[dat$pericent] <- as.character(seqnames(regions[pericent]))
dat$pericent[dat$pericent == FALSE] <- "Chr0"
save(dat, file = paste0(outDir, "df_for_GLM.RData"))

#############################################################
# run binomial GLM model, link function is logit, plot data #
#############################################################
# predict function produces predicted values 
# Tom uses sumspo - why not meanspo - because width is a variable?

glmCO <- glm(co~band+sumREC8_ChIP*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict1 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_REC8_ChIP_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_REC8_ChIP_coeffs.csv"))
save(predict1, file = paste0(outDir, "GLM_REC8_ChIP_predict.RData"))

glmCO <- glm(co~band+sumnuc*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict2 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_nuc_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_nuc_coeffs.csv"))
save(predict2, file = paste0(outDir, "GLM_nuc_predict.RData"))

glmCO <- glm(co~band+sumspo*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict3 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_SPO11-oligos_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_SPO11-oligos_coeffs.csv"))
save(predict3, file = paste0(outDir, "GLM_SPO11-oligos_predict.RData"))

glmCO <- glm(co~band+sumMSH4_ChIP*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict4 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_MSH4_ChIP_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_MSH4_ChIP_coeffs.csv"))
save(predict4, file = paste0(outDir, "GLM_MSH4_ChIP_predict.RData"))

glmCO <- glm(co~band+(sumREC8_ChIP+sumMSH4_ChIP+sumnuc+sumspo)*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict5 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_REC8_MSH4_nuc_SPO11_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_REC8_MSH4_nuc_SPO11_coeffs.csv"))
save(predict5, file = paste0(outDir, "GLM_REC8_MSH4_nuc_SPO11_predict.RData"))

###################################################
# plot predicted CO overlaps for annotation class #
###################################################

#regid <- c(list(geneOverlap, promoterOverlap, terminatorOverlap), lapply(colnames(tesOverlap), function(tesn){tesOverlap[,tesn]}))
#names(regid) <- c("gene", "promoter", "terminator", colnames(tesOverlap))
#pdf(file = paste0(outDir, "GLM_REC8_ChIP_annotation_CO_boxplot.pdf"))
#boxplot(lapply(regid, function(x) sapply(1:100, function(ii) sum(rbinom(sum(x), 1, prob = predict1[x]))/sum(dat$width[x])*1e6)), main = "predicted CO/Mb by annotation")
#tco <- sapply(regid, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
#points(x = 1:length(regid), y = tco, col = "red", pch = 19)
#dev.off()

#####################################################
# plot predicted CO for regions by REC8_ChIP hexile #
#####################################################

#sort(dat$meanREC8_ChIP)[round(0:7*(nrow(dat)/7))]
#[1]    4.384697    6.054059    7.394007    8.744343   10.205637   12.277201
#[7] 3242.308239
#levels(cut(dat$meanREC8_ChIP, breaks = sort(dat$meanREC8_ChIP)[round(0:7*(nrow(dat)/7))]))
#[1] "(4.38,6.05]"     "(6.05,7.39]"     "(7.39,8.74]"     "(8.74,10.2]"    
#[5] "(10.2,12.3]"     "(12.3,3.24e+03]"
#length(dat$meanREC8_ChIP[dat$meanREC8_ChIP < 4.384697])
#[1] 68613
# omits 68613 of 480348 regions for which meanREC8_ChIP < 4.384697
# so changed below from breaks = sort(dat$meanREC8_ChIP)[round(0:7*(nrow(dat)/7))] to breaks = c(0, sort(dat$meanREC8_ChIP)[round(0:6*(nrow(dat)/6))]) 
# so that...
#levels(cut(dat$meanREC8_ChIP, breaks = c(0, sort(dat$meanREC8_ChIP)[round(0:6*(nrow(dat)/6))])))
#[1] "(0,4.75]"        "(4.75,6.52]"     "(6.52,8.07]"     "(8.07,9.68]"    
#[5] "(9.68,11.8]"     "(11.8,3.24e+03]"
ssID <- split(1:nrow(dat), cut(dat$meanREC8_ChIP, breaks = c(0, sort(dat$meanREC8_ChIP)[round(0:6*(nrow(dat)/6))])))
pdf(file = paste0(outDir, "GLM_REC8_ChIP_hexiles_CO_boxplot.pdf"))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict1[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+REC8 ChIP*(annotation)", xlab = "REC8 ChIP hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

ssID <- split(1:nrow(dat), cut(dat$meanREC8_ChIP, breaks = sort(dat$meanREC8_ChIP)[round(0:7*(nrow(dat)/7))]))
pdf(file = paste0(outDir, "GLM_REC8_ChIP_hexiles_CO_boxplot_lessnoise?.pdf"))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict1[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+REC8 ChIP*(annotation)", xlab = "REC8 ChIP hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

#######################################################
# plot predicted CO for regions by nucleosomes hexile #
#######################################################

ssID <- split(1:nrow(dat), cut(dat$meannuc, breaks = c(0, sort(dat$meannuc)[round(0:6*(nrow(dat)/6))])))
pdf(file = paste0(outDir, "GLM_nuc_hexiles_CO_boxplot.pdf"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict2[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+nucleosomes*(annotation)", xlab = "Nucleosomes hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"), ylim = c(min(tco), 40))
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

ssID <- split(1:nrow(dat), cut(dat$meannuc, breaks = sort(dat$meannuc)[round(0:7*(nrow(dat)/7))]))
pdf(file = paste0(outDir, "GLM_nuc_hexiles_CO_boxplot_lessnoise?.pdf"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict2[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+nucleosomes*(annotation)", xlab = "Nucleosomes hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"), ylim = c(min(tco), max(tco)))
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

##########################################################
# plot predicted CO for regions by SPO11-1-oligos hexile #
##########################################################

ssID <- split(1:nrow(dat), cut(dat$meanspo, breaks = c(0, sort(dat$meanspo)[round(0:6*(nrow(dat)/6))])))
pdf(file = paste0(outDir, "GLM_SPO11-1-oligos_hexiles_CO_boxplot.pdf"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict3[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+SPO11-1-oligos*(annotation)", xlab = "SPO11-1-oligos hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"), ylim = c(min(tco), 53))
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

ssID <- split(1:nrow(dat), cut(dat$meanspo, breaks = sort(dat$meanspo)[round(0:7*(nrow(dat)/7))]))
pdf(file = paste0(outDir, "GLM_SPO11-1-oligos_hexiles_CO_boxplot_lessnoise?.pdf"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict3[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+SPO11-1-oligos*(annotation)", xlab = "SPO11-1-oligos hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

#####################################################
# plot predicted CO for regions by MSH4_ChIP hexile #
#####################################################

ssID <- split(1:nrow(dat), cut(dat$meanMSH4_ChIP, breaks = c(0, sort(dat$meanMSH4_ChIP)[round(0:6*(nrow(dat)/6))])))
pdf(file = paste0(outDir, "GLM_MSH4_ChIP_hexiles_CO_boxplot.pdf"))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict4[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+MSH4 ChIP*(annotation)", xlab = "MSH4 ChIP hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

ssID <- split(1:nrow(dat), cut(dat$meanMSH4_ChIP, breaks = sort(dat$meanMSH4_ChIP)[round(0:7*(nrow(dat)/7))]))
pdf(file = paste0(outDir, "GLM_MSH4_ChIP_hexiles_CO_boxplot_lessnoise?.pdf"))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict4[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+MSH4 ChIP*(annotation)", xlab = "MSH4 ChIP hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

##########################################
# plot predicted CO for chromosome bands #
##########################################

pdf(file = paste0(outDir, "GLM_REC8_ChIP_chrbands_CO_boxplot.pdf"), height = 6, width = 9)
boxplot(lapply(split(1:nrow(dat), dat$band), function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict1[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+REC8 ChIP*(annotation)", xlab = "2-Mb chromosome bands", ylab = "Number of crossovers per Mb")
ssb <- split(1:nrow(dat), dat$band)
points(x = unique(dat$band), y = sapply(ssb, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6), col = "red", pch = 19)
abline(v = c(8, 17, 31, 37, 50), lty = 2, lwd = 2.5)
legend("topleft", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

pdf(file = paste0(outDir, "GLM_REC8_MSH4_nuc_SPO11_chrbands_CO_boxplot.pdf"), height = 6, width = 9)
boxplot(lapply(split(1:nrow(dat), dat$band), function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict5[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+(REC8 ChIP+MSH4 ChIP+nucleosomes+SPO11-1-oligos)*(annotation)", xlab = "2-Mb chromosome bands", ylab = "Number of crossovers per Mb")
ssb <- split(1:nrow(dat), dat$band)
points(x = unique(dat$band), y = sapply(ssb, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6), col = "red", pch = 19)
abline(v = c(8, 17, 31, 37, 50), lty = 2, lwd = 2.5)
legend("topleft", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

#########################################
# plot predicted CO for each SNP region #
#########################################

#pdf(file = paste0(outDir, "GLM_REC8_ChIP_snpreg_CO.pdf"))
#xy <- cbind(x = cumsum(width(regions)), y = (predict1))
#plot(xy[xy[,2]>0.01,], pch = ".", main = "Likelihood of CO at each SNP window");
#abline(v = cumsum(chrlens), col = "red", lty = 3)
#abline(v = c(0, cumsum(chrlens)[-5])+start(centGR), col = "blue", lty = 3)
#abline(v = c(0, cumsum(chrlens)[-5])+end(centGR), col = "blue", lty = 3)
#dev.off()

sessionInfo()


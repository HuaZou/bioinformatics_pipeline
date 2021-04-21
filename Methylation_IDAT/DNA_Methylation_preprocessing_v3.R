#!/usr/bin/Rscript

###Methylation Pre-processing     

###############################################################
###########            loading R packages          ############
###############################################################
library(knitr)
library(dplyr)
library(limma)
library(minfi)
library(missMethyl)
library(minfiData)
library(stringr)
library(ENmix)
library(RColorBrewer)
library(optparse)

# get options
option_list = list(
  make_option(c("-f", "--folder"), type="character", default=".", 
        help="folder with idat files [default= %default]", metavar="character"),
  make_option(c("-t", "--target"), type="character", default=".", 
        help="target file with sample information [default= %default]", metavar="character"),
  make_option(c("-c", "--chip"), type="character", default=".", 
        help="the type of DNA chip array file [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
        help="output directory name [default= %default]", metavar="character")
); 
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
print(opt)

###############################################################
###########      Creating output directory         ############
###############################################################
create_dir <- function(dirpath){
  
  if(!dir.exists(dirpath)){
    dir.create(dirpath, recursive = TRUE)
  }  
}
# current dir
current_dir <- getwd()
# Result & QC dir 
out_dir <- opt$out
qc_dir <- paste0(out_dir, "/QC")
create_dir(qc_dir)
# DataSet dir
set_dir <- paste0(out_dir, "/Set")
create_dir(set_dir)
# Matrix dir
raw_dir <- paste0(out_dir, "/Matrix")
create_dir(raw_dir)

###############################################################
#####   Chip Array annotation and cross reactive files  #######
###############################################################
if(opt$chip == "850k"){
  # For HumanMethylationEPIC bead chip array files (850k) load the following two packages and two files
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  crossreac_probes <- read.csv("/disk/user/zouhua/pipeline/methylseq_IDAT/util/13059_2016_1066_MOESM1_ESM_cross-reactive_probes.csv")
}else if(chip == "450k"){
  # For HumanMethylation450 bead chip array files (450k) load the following packages and two files
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  crossreac_probes <- read.csv("/disk/user/zouhua/pipeline/methylseq_IDAT/util/48639-non-specific-probes-Illumina450k.csv")  
}

###############################################################
###########            Reading idat files          ############
###############################################################
# list the files
list.files(opt$folder, recursive = TRUE)
# read in the sample sheet for the experiment
targets <- read.csv(opt$target)
head(targets)

###############################################################
###########          Creating RGset files          ############
###############################################################
RGSet <- read.metharray.exp(targets = targets)
targets$ID <- paste(targets$Sample_Name)
sampleNames(RGSet) <- targets$ID
print(RGSet)

###############################################################
###########             Quality Checks             ############
###############################################################
# C.1. Plot quality control plots (package ENmix)
setwd(qc_dir)
plotCtrl(RGSet)

# C.2. Make PDF QC report (package minfi)
# To include colouring samples by variable such as sentrix position include argument: sampGroups=targets$sentrix_pos
setwd(current_dir)
qcReport_pdf <- paste0(qc_dir, "/qcReport.pdf")
qcReport(RGSet, sampNames = targets$Sample_Name, pdf = qcReport_pdf)

# C.3. Make pre-normalisation beta density plots
# Here samples are coloured by sentrix position
nb.levels <- length(unique(targets$Array))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.levels)
jpeg(paste0(qc_dir, "/UnormalisedBetaDensityPlot_bySentrixPosition.jpg"), width = 800, height = 800)
densityPlot(RGSet, sampGroups = targets$Array, pal = mycolors, ylim = c(0, 5))
dev.off()

# C.4. Create MethylSet, RatioSet, then a GenomicRatioSet for further QC
MSet <- preprocessRaw(RGSet)
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet, mergeManifest = T)

# C.5. Perform QC on MSet and plot methylated versus unmethylated intensities
qc <- getQC(MSet)
pdf(paste0(qc_dir, "/Meth-unmeth_intensities.pdf"), height = 4, width = 4)
par(mfrow = c(1, 1), family = "Times", las = 1)
plotQC(qc) # If U and/or M intensity log medians are <10.5, sample is by default of bad quality
dev.off()

# C.6. remove any samples with meth or unmeth intensity < 10.5, then remove from all previous files
if(sum(rowMeans(as.data.frame(qc)) < 10.5 ) > 0) {
  print(paste("Warning: sample", rownames(qc)[rowMeans(as.data.frame(qc)) < 10.5], 
              "has mean meth-unmeth signal intensity <10.5. Remove sample."))
}
keep.samples <- apply(as.data.frame(qc), 1, mean) > 10.5
RGSet_remain <- RGSet[, keep.samples]
MSet_remain <- MSet[, keep.samples]
ratioSet_remain <- ratioSet[, keep.samples]
gset_remain <- gset[, keep.samples]
targets_remain <- targets[keep.samples, ]

###############################################################
########    Calculate detection p values and plot      ########
###############################################################
# Here the samples are coloured by sentrix ID to check for poor performing chips
detP <- detectionP(RGSet_remain)
pdf(paste0(qc_dir, "/detection_pvalues.pdf"), height = 4, width = 4)
par(mfrow = c(1, 1), family = "Times", las = 1)
barplot(colMeans(detP), 
        col = as.numeric(factor(targets$Slide)), 
        las = 2, 
        cex.names = 0.8, 
        ylim = c(0, 0.05), 
        ylab = "Mean detection p-values")
abline(h = 0.01 ,col = "red")
dev.off()

if(sum(colMeans(detP) > 0.01) > 0){
  print(paste("Warning: sample", names(colMeans(detP))[colMeans(detP) > 0.01], 
              "has >0.01 mean detection p-value. Remove sample"))
}
# If required: remove any samples with detection p value > 0.01, then remove from all previous files
keep.samples_v2 <- apply(detP, 2, mean) < 0.01 
RGSet_remain_v2 <- RGSet_remain[, keep.samples_v2]
MSet_remain_v2 <- MSet_remain[, keep.samples_v2]
ratioSet_remain_v2 <- ratioSet_remain[, keep.samples_v2]
gset_remain_v2 <- gset_remain[, keep.samples_v2]
targets_remain_v2 <- targets_remain[keep.samples_v2, ]
glimpse(targets_remain_v2)

###############################################################
########             output Raw Set files              ########
###############################################################
save(RGSet, file = paste0(set_dir, "/RGSet_Raw.RData"))
save(MSet, file = paste0(set_dir, "/MSet_Raw.RData"))
save(gset, file = paste0(set_dir, "/gset_Raw.RData"))

save(RGSet_remain, file = paste0(set_dir, "/RGSet_trim_intensity.RData"))
save(MSet_remain, file = paste0(set_dir, "/MSet_trim_intensity.RData"))
save(gset_remain, file = paste0(set_dir, "/gset_trim_intensity.RData"))

save(RGSet_remain_v2, file = paste0(set_dir, "/RGSet_trim_detection.RData"))
save(MSet_remain_v2, file = paste0(set_dir, "/MSet_trim_detection.RData"))
save(gset_remain_v2, file = paste0(set_dir, "/gset_trim_detection.RData"))

size.RGSet_remain_v2 <- as.data.frame(t(dim(RGSet_remain_v2)))
size.MSet_remain_v2 <- as.data.frame(t(dim(MSet_remain_v2)))
size.gset_remain_v2 <- as.data.frame(t(dim(gset_remain_v2)))

#########################################################################
# Create raw unormalised and unfiltered beta table for later comparison #
#########################################################################
betaRaw <- getBeta(gset_remain_v2)
dim(betaRaw)
mRaw <- getM(gset_remain_v2)
dim(mRaw)

save(betaRaw, file = paste0(raw_dir, "/betaRaw.RData"))
save(mRaw, file = paste0(raw_dir, "/mRaw.RData"))

###############################################################
########                  Normalization                ########
###############################################################
RGSet_norm <- preprocessFunnorm(RGSet_remain_v2)
size.RGSet_norm <- as.data.frame(t(dim(RGSet_norm)))
save(RGSet_norm, file = paste0(set_dir, "/RGSet_normalization.RData"))

# visualise what the data looks like before and after normalisation
pdf(paste0(raw_dir, "/normalisation_status.pdf"), height = 4, width = 4)
par(mfrow=c(1,2))
densityPlot(RGSet_remain_v2, sampGroups=targets$Sample_Source, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Source)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(RGSet_norm), sampGroups=targets$Sample_Source,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Source)), 
       text.col=brewer.pal(8, "Dark2"))
dev.off()

###############################################################
########               Data exploration                ########
###############################################################
# MDS plots to look at largest sources of variation
pal <- brewer.pal(8,"Dark2")
pdf(paste0(raw_dir, "/Data_exploration_PCA.pdf"), height = 4, width = 4)
par(mfrow=c(1, 2))
plotMDS(getM(RGSet_norm), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)
plotMDS(getM(RGSet_norm), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

# Examine higher dimensions to look at other sources of variation
pdf(paste0(raw_dir, "/Data_exploration_PCA_v2.pdf"), height = 4, width = 6)
par(mfrow=c(1, 3))
plotMDS(getM(RGSet_norm), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
       cex=0.7, bg="white")

plotMDS(getM(RGSet_norm), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(RGSet_norm), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
dev.off()

###############################################################
########                   filtering                   ########
###############################################################
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(RGSet_norm), rownames(detP)), ] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(RGSet_norm) 
table(keep)
RGSet_normFlt <- RGSet_norm[keep, ]
RGSet_normFlt

# removing duplicated 
pData <- pData(RGSet_normFlt)
pData2 <- pData[order(pData(RGSet_normFlt)$Sample_Name), ]
dim(pData2)
head(pData2)

dups <- pData2$Sample_Name[duplicated(pData2$Sample_Name)]
if( length(dups) > 0 ){
    whdups     <- lapply(dups, function(x){which(pData2$Sample_Name == x)})
    whdups2rem <- sapply(1:length(dups), function(i) rbinom(1, 1, prob = 1/length(whdups[[1]])))+1
    torem <- sapply(1:length(whdups), function(i){whdups[[i]][whdups2rem[i]]})
    pData2 <- pData2[-torem, ]
}
RGSet_norm_duplicate <- RGSet_normFlt[, rownames(pData2)]
save(RGSet_norm_duplicate, file = paste0(set_dir, "/RGSet_norm_duplicate.RData"))
size.detP <- as.data.frame(t(dim(RGSet_norm_duplicate)))

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(RGSet_norm_duplicate) %in% ann$Name[ann$chr %in% 
                                                        c("chrX","chrY")])
table(keep)
RGSet_norm_duplicate <- RGSet_norm_duplicate[keep, ]

# remove probes with SNPs at CpG site
RGSetFlt <- dropLociWithSnps(RGSet_norm_duplicate)
RGSetFlt

# exclude cross reactive probes 
keep <- !(featureNames(RGSetFlt) %in% crossreac_probes$TargetID)
table(keep)
RGSetFlt <- RGSetFlt[keep,] 
RGSetFlt

pdf(paste0(raw_dir, "/filtered_status.pdf"), height = 6, width = 9)
par(mfrow=c(1, 2))
plotMDS(getM(RGSet_norm), top=1000, gene.selection="common", main="Normalized",
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
plotMDS(getM(RGSetFlt), top=1000, gene.selection="common", main="Normalized_filter",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

###############################################################
########                Normalization matrix           ########
###############################################################
betaNorm <- getBeta(RGSetFlt)
dim(betaNorm)
size.betaNorm <- as.data.frame(t(dim(betaNorm)))

mNorm <- getM(RGSetFlt)
dim(mNorm)
size.mNorm <- as.data.frame(t(dim(mNorm)))

# visualization
pdf(paste0(raw_dir, "/BetaM_norm.pdf"), height = 4, width = 4)
par(mfrow=c(1,2))
densityPlot(betaNorm, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mNorm, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

# Check for NA and -Inf values
betaRaw.na <- betaRaw[!complete.cases(betaRaw), ]
dim(betaRaw.na) 
NAbetas <- rownames(betaRaw.na)
betaNAsNorm <- betaNorm[rownames(betaNorm) %in% NAbetas, ] # Should be [1] 0 x ncol(betaNorm)

# Check for infinite values in m table, replace in the no infinite values m table
TestInf <- which(apply(mNorm, 1, function(i)sum(is.infinite(i))) > 0)
save(TestInf, file = paste0(set_dir, "/InfiniteValueProbes.RData"))
mNoInf <- mNorm
mNoInf[!is.finite(mNoInf)] <- min(mNoInf[is.finite(mNoInf)])
TestInf2 <- which(apply(mNoInf, 1, function(i) sum(is.infinite(i))) > 0)
TestInf2 # Should be named integer(0)
size.mNoInf <- as.data.frame(t(dim(mNoInf)))

# Write tables
write.table(targets_remain_v2, file=paste0(raw_dir, "/TargetsFile.csv"), sep=",", col.names=NA)
write.table(betaNorm, file=paste0(raw_dir, "/NormalisedFilteredBetaTable.csv"), sep=",", col.names=NA)
write.table(mNorm, file=paste0(raw_dir, "/NormalisedFilteredMTable.csv"), sep=",", col.names=NA)
write.table(mNoInf, file=paste0(raw_dir, "/NormalisedFilteredMTable_noInf.csv"), sep=",", col.names=NA)

###############################################################
########              Create dimensions table          ########
###############################################################
DimensionsTable <- rbind(size.RGSet_remain_v2, size.MSet_remain_v2, size.gset_remain_v2, 
                         size.RGSet_norm, size.detP,
                         size.betaNorm, size.mNorm, size.mNoInf)
rownames(DimensionsTable) <- c("RGset_size", "MSet_size", "gset_size", 
                               "fun_size", "detP_lost", 
                               "beta_size", "m_size", "mNoInf_size")
colnames(DimensionsTable) <- c("Probe_number", "Sample_number")
write.table(DimensionsTable, file=paste0(out_dir, "/DimensionsTable.txt"), sep="\t", col.names=NA)

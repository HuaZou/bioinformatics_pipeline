#!/usr/bin/Rscript

###############################################################
###########   Methylation Pre-processing script    ############
# https://github.com/IARCbioinfo/Methylation_analysis_scripts #
###############################################################

library(optparse)
# get options
option_list = list(
  make_option(c("-f", "--folder"), type="character", default=".", 
        help="folder with idat files [default= %default]", metavar="character"),
  make_option(c("-t", "--target"), type="character", default=".", 
        help="target file with sample information [default= %default]", metavar="character"),
  make_option(c("-p", "--sep"), type="character", default="\t", 
        help="field separator for target file [default= %default]", metavar="character"),
  make_option(c("-c", "--crossreac"), type="character", default=NULL, 
        help="file with list of cross-reactive probes [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
        help="output directory name [default= %default]", metavar="character"),
  make_option(c("-s", "--snp_filter"), action="store_true", default=TRUE, type="logical", 
        help="filter SNPs-associated probes [default= %default]", metavar="logical"),
  make_option(c("-m", "--multimodal_filter"), action="store_true", default=FALSE, type="logical", 
        help="filter multimodal probes [default= %default]", metavar="logical")
); 
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
print(opt)

create_dir <- function(dirpath){
  
  if(!dir.exists(dirpath)){
    dir.create(dirpath, recursive = TRUE)
  }  
}

###############################################################
###########            loading R packages          ############
###############################################################
library(minfi)
library(ENmix)
library(wateRmelon)
library(MASS)
library(broom)
library(dplyr)
library(RColorBrewer)
# For HumanMethylationEPIC bead chip array files (850k) load the following two packages and two files
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

###############################################################
###########            Reading idat files          ############
###############################################################
targets <- read.table(opt$target, header = T, sep = opt$sep)
colnames(targets) <- tolower(colnames(targets))
if(!any(is.element(colnames(targets), c("sample_id", "sentrix_id", "sentrix_position")))){
  stop("Please checking the input file")
}
targets$barcode <- with(targets, paste(sentrix_id, sentrix_position, sep = "_"))
rownames(targets) <- as.character(targets$barcode)
targets$Basename <- rownames(targets)
glimpse(targets)

###############################################################
###########          Creating RGset files          ############
###############################################################
RGSet <- read.metharray.exp(base = opt$folder, 
                            targets = targets, 
                            recursive = T, 
                            extended = T)
targets$ID <- paste(targets$sample_id)
sampleNames(RGSet) <- targets$ID
print(RGSet)

###############################################################
###########             Quality Checks             ############
###############################################################
# current dir
current_dir <- getwd()

# result dir 
out_dir <- opt$out
qc_dir <- paste0(out_dir, "/QC")
create_dir(qc_dir)

# C.1. Plot quality control plots (package ENmix)
setwd(qc_dir)
plotCtrl(RGSet)

# C.2. Make PDF QC report (package minfi)
# To include colouring samples by variable such as sentrix position include argument: sampGroups=targets$sentrix_pos
setwd(current_dir)
qcReport_pdf <- paste0(qc_dir, "/qcReport.pdf")
qcReport(RGSet, sampNames = targets$sample_id, pdf = qcReport_pdf)

# C.3. Make pre-normalisation beta density plots
# Here samples are coloured by sentrix position
nb.levels <- length(unique(targets$sentrix_position))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.levels)
jpeg(paste0(qc_dir, "/UnormalisedBetaDensityPlot_bySentrixPosition.jpg"), width = 800, height = 800)
densityPlot(RGSet, sampGroups = targets$sentrix_pos, pal = mycolors, ylim = c(0, 5))
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
  print(paste("Warning: sample", rownames(qc)[rowMeans(as.data.frame(qc))<10.5], 
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
        col = as.numeric(factor(targets$sentrix_id)), 
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
set_dir <- paste0(out_dir, "/Set")
create_dir(set_dir)

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

raw_dir <- paste0(out_dir, "/Matrix")
create_dir(raw_dir)
save(betaRaw, file = paste0(raw_dir, "/betaRaw.RData"))
save(mRaw, file = paste0(raw_dir, "/mRaw.RData"))

###############################################################
########                  Normalization                ########
###############################################################
RGSet_norm <- preprocessFunnorm(RGSet_remain_v2)
size.RGSet_norm <- as.data.frame(t(dim(RGSet_norm)))
save(RGSet_norm, file = paste0(set_dir, "/RGSet_normalization.RData"))

###############################################################
########              Remove duplicates                ########
###############################################################
pData <- pData(RGSet_norm)
pData2 <- pData[order(pData(RGSet_norm)$sample_id), ]
dim(pData2)
head(pData2)

dups <- pData2$sample_id[duplicated(pData2$sample_id)]
if( length(dups) > 0 ){
    whdups     <- lapply(dups, function(x){which(pData2$sample_id == x)})
    whdups2rem <- sapply(1:length(dups), function(i) rbinom(1, 1, prob = 1/length(whdups[[1]])))+1
    torem <- sapply(1:length(whdups), function(i){whdups[[i]][whdups2rem[i]]})
    pData2 <- pData2[-torem, ]
}
RGSet_norm_duplicate <- RGSet_norm[, rownames(pData2)]
save(RGSet_norm_duplicate, file = paste0(set_dir, "/RGSet_norm_duplicate.RData"))
size.detP <- as.data.frame(t(dim(RGSet_norm_duplicate)))

###############################################################
########                Normalization matrix           ########
###############################################################
betaNorm <- getBeta(RGSet_norm_duplicate)
dim(betaNorm)
size.betaNorm <- as.data.frame(t(dim(betaNorm)))

mNorm <- getM(RGSet_norm_duplicate)
dim(mNorm)
size.mNorm <- as.data.frame(t(dim(mNorm)))

raw_dir <- paste0(out_dir, "/Matrix")
create_dir(raw_dir)

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

#!/usr/bin/Rscript

###############################################################
###########   Methylation Pre-processing script    ############
# https://github.com/IARCbioinfo/Methylation_analysis_scripts #
###############################################################

library("optparse")
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

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

# load R packages
library(minfi)
library(ENmix)
library(wateRmelon)
library(MASS)
library(broom)
library(RColorBrewer)

# For HumanMethylationEPIC bead chip array files (850k) load the following two packages and two files
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# crossreac file: "13059_2016_1066_MOESM1_ESM_cross-reactive_probes.csv"

# For HumanMethylation450 bead chip array files (450k) load the following packages and two files
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# crossreac file: "48639-non-specific-probes-Illumina450k.csv"

# A) Load targets file
targets <- read.table(opt$target, header=T, sep=opt$sep)
# If targets file is a csv use: targets <- read.csv(opt$targets, header=TRUE)
colnames(targets) <- tolower(colnames(targets))
targets$barcode <- paste(targets$sentrix_id, targets$sentrix_pos, sep="_")
rownames(targets) <- as.character(targets$barcode)
targets$Basename <- rownames(targets)
head(targets)

# B) Create RGSet
RGSet <- read.metharray.exp(base = opt$folder, targets = targets, recursive = T, extended = T)
# Give RGSet meaningful names (here targets$sample_id)
targets$ID = paste(targets$sample_id)
sampleNames(RGSet) = targets$ID
print(RGSet)

# C) Quality Checks
# Create output dir and QC dir within output dir
outdir <- opt$out
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)
dir.create("QC/", showWarnings = FALSE)

# C.1. Plot quality control plots (package ENmix)
plotCtrl(RGSet)

# C.2. Make PDF QC report (package minfi)
# To include colouring samples by variable such as sentrix position include argument: sampGroups=targets$sentrix_pos
qcReport(RGSet, sampNames = targets$sample_id, pdf = "QC/qcReport.pdf")

# C.3. Make pre-normalisation beta density plots
# Here samples are coloured by sentrix position
nb.levels <- length(unique(targets$sentrix_pos))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.levels)
jpeg(paste("QC/UnormalisedBetaDensityPlot_bySentrixPosition.jpg", sep="/"), width=800, height=800)
densityPlot(RGSet, sampGroups = targets$sentrix_pos, pal=mycolors, ylim=c(0,5))
dev.off()

# C.4. Create MethylSet, RatioSet, then a GenomicRatioSet for further QC
MSet <- preprocessRaw(RGSet)
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet, mergeManifest=T)

# C.5. Perform QC on MSet and plot methylated versus unmethylated intensities
qc <- getQC(MSet)
pdf("QC/Meth-unmeth_intensities.pdf",h=4,w=4)
par(mfrow=c(1,1),family="Times",las=1)
plotQC(qc) # If U and/or M intensity log medians are <10.5, sample is by default of bad quality
dev.off()

if(sum(rowMeans(as.data.frame(qc))<10.5)>0) print(paste("Warning: sample", rownames(qc)[rowMeans(as.data.frame(qc))<10.5], "has mean meth-unmeth signal intensity <10.5. Remove sample."))
# If required: remove any samples with meth or unmeth intensity < 10.5, then remove from all previous files
keep.samples <- apply(as.data.frame(qc),1,mean) > 10.5
RGSet <- RGSet[, keep.samples]
MSet <- MSet[, keep.samples]
ratioSet <- ratioSet[, keep.samples]
gset <- gset[, keep.samples]
targets <- targets[keep.samples, ]

# D) Calculate detection p values and plot 
# Here the samples are coloured by sentrix ID to check for poor performing chips
detP <- detectionP(RGSet)
pdf("QC/detection_pvalues.pdf",h=4,w=4)
par(mfrow=c(1,1),family="Times",las=1)
barplot(colMeans(detP), col=as.numeric(factor(targets$sentrix_id)), las=2, cex.names=0.8, ylim=c(0,0.05), ylab="Mean detection p-values")
abline(h=0.01,col="red")
dev.off()

if(sum(colMeans(detP)>0.01)>0 ) print(paste("Warning: sample", names(colMeans(detP))[colMeans(detP)>0.01], "has >0.01 mean detection p-value. Remove sample"))
# If required: remove any samples with detection p value > 0.01, then remove from all previous files
keep.samples <- apply(detP,2,mean) < 0.01 
RGSet <- RGSet[,keep.samples]
MSet <- MSet[,keep.samples]
ratioSet <- ratioSet[,keep.samples]
gset <- gset[,keep.samples]
targets <- targets[keep.samples,]

# E) Write files
save(RGSet, file="RGSet.RData")
save(MSet, file="MSet.RData")
save(gset, file="gset.RData")

size.RGSet <- as.data.frame(t(dim(RGSet)))
size.MSet <- as.data.frame(t(dim(MSet)))
size.gset <- as.data.frame(t(dim(gset)))

# F) Create raw unormalised and unfiltered beta table for later comparison
betaRaw <- getBeta(gset)
dim(betaRaw)

# G) Surrogate variables and Principal components
# Surrogate variables derived from intensity data for non-negative internal control probes.
# These variables can be modeled in association analysis to adjust for experimental batch effects.
sva <- ctrlsva(RGSet, percvar = 0.9, flag = 1)
save(sva, file="sva90.RData")

K = ncol( sva )
# When data comes from different centers:
#lmsvaFull = lapply(1:K, function(i) lm(sva[,i]~Sentrix_id+Sentrix_position+Center,
#                                       data.frame("Sentrix_id"=as.factor(pData(RGSet)$sentrix_id),"Sentrix_position"=as.factor(pData(RGSet)$sentrix_position),"Center"=as.factor(pData(RGSet)$center) ) ) )
# When data comes from the same center or unknown center
lmsvaFull = lapply(1:K, function(i) lm(sva[,i]~Sentrix_id+Sentrix_position,
                                       data.frame("Sentrix_id"=as.factor(pData(RGSet)$sentrix_id),"Sentrix_position"=as.factor(pData(RGSet)$sentrix_position) ) ) )
lmsvaRed = vector("list",K)
for(i in 1:K){
  lmtmp = lmsvaFull[[i]]
  while(1){
    dttmp = dropterm(lmtmp,test = "F")
    if(max(dttmp$`Pr(F)`,na.rm = T)> (0.05) ) ttmp = rownames(dttmp)[which.max(dttmp$`Pr(F)`)]
    else break
    lmtmp = update(lmtmp, paste(".~. - ",ttmp) )
    print(dttmp)
    print(lmtmp)
  }
  lmsvaRed[[i]] = lmtmp
}

for(i in 1:K) write.csv( tidy( anova(lmsvaFull[[i]]) ) ,file = paste("QC/ANOVAfull_sva",i,".csv",sep=""),quote = F,row.names = F)
for(i in 1:K) write.csv( tidy( anova(lmsvaRed[[i]]) ) ,file = paste("QC/ANOVAreduced_sva",i,".csv",sep=""),quote = F,row.names = F)

#require(beeswarm)
pdf("QC/SVA.pdf",h=3*K,w=3*K)
par(mfrow=c(K,K),family="Times",las=1)
for(i in 1:K){
  for(j in 1:K){
    plot(sva[,j], sva[,i], col=rainbow(20)[as.factor(pData(RGSet)$sentrix_id)],pch=as.numeric(as.factor(pData(RGSet)$sentrix_position)) ,   
         xlab=paste("SV",j),ylab=paste("SV",i),main="Effects of Sentrix id (color) & Sentrix position (shape)")
  }
}
#legend("bottomright",legend = c( levels(as.factor(metaM[,22])),levels(as.factor(metaM[,13]))), col=c(colClustLNEN,1,1),pch=c(16,16,16,16,17) )
#beeswarm( sva[,7], pwcol=as.factor(metaM[as.character(pData(fun)$Ms_IDs),2]) , pch=16 ,xlab="",ylab="SV 7",main="Effect of Center of origin",method = "square",spacing = 2)
#legend("bottomright",legend = levels(as.factor(metaM[,2])), col=1:8,pch=16)
dev.off() ## to check

# H) Predict sex, plot predicted sex (and if available, plot clinical sex and identify mismatches)
object = getSex(gset, cutoff = -2) # Default cutoff is -2
pdf("Sex_predicted.pdf",h=10,w=10)
plot(x = object$xMed, y = object$yMed, type = "n", xlab = "X chr, median total intensity (log2)", ylab = "Y chr, median total intensity (log2)")
text(x = object$xMed, y = object$yMed, labels = targets$sample_id, col = ifelse(object$predictedSex == "M", "deepskyblue", "deeppink3"))
legend("bottomleft", c("M", "F"), col = c("deepskyblue", "deeppink3"), pch = 16)
dev.off()

object2 <- as.data.frame(object) # To create clinical sex plot
object2 <- merge(object2, targets, by.x="row.names", by.y="sample_id")
pdf("Sex_clinical.pdf",h=10,w=10)
plot(x = object2$xMed, y = object2$yMed, type = "n", xlab = "X chr, median total intensity (log2)", ylab = "Y chr, median total intensity (log2)")
text(x = object2$xMed, y = object2$yMed, labels = object2$Row.names, col = ifelse(object2$sex == "M", "deepskyblue", "deeppink3"))
legend("bottomleft", c("M", "F"), col = c("deepskyblue", "deeppink3"), pch = 16)
dev.off()

# Bind the predicted sex to the targets file and identify any mismatches 
targets$predSex <- object$predictedSex
targets[targets$sex != targets$predSex, ]

# I) Normalization
# includes NOOB background/dye correction
fun <- preprocessFunnorm(RGSet) # If clinical sex is available can use preprocessFunnorm(RGSet, sex=targets$sex), default of sex=NULL uses function getSex to guess the sex
size.fun <- as.data.frame(t(dim(fun)))
save(fun, file = "Fun.RData")

# Post normalisation beta density plots
nb.levels <- length(unique(targets$sentrix_pos))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.levels)
jpeg(paste("NormalisedBetaDensityPlot_bySentrixPosition.jpg",sep="/"), width=800, height=800)
densityPlot(getBeta(fun), sampGroups = targets$sentrix_pos, pal=mycolors, ylim=c(0,5))
dev.off()

# J) Optional processes
### Optional: Age Prediction
predicted_age <- agep(minfi::getBeta(fun))
pData$predictedAge = predicted_age
write.csv( data.frame("Sample"=as.character(pData$sample_id),
                      "Clinical_data"=as.character(pData$age),
                      "Prediction"=predicted_age,row.names = pData$barcode), file = "Age_pred.csv",quote = F)

pdf("Age_pred.pdf")
par(family="Times",las=1)
plot( pData(fun)$age, predicted_age,pch=16,xlim=c(0,120),ylim=c(0,140), xlab="Age", ylab="Predicted age")
abline(lm(pData(fun)$age~predicted_age[,1]), col="red")
ct = cor.test(as.numeric(pData(fun)$age), as.numeric(predicted_age))
text( 0,10 , paste0("r=",ct$estimate) ,pos = 4 )
text( 0,0 , paste0("p=",ct$p.value) ,pos = 4 )
dev.off()

### Optional: Smoking "prediction"
cg05575921 <- minfi::getBeta(fun)["cg05575921", ] # within AHRR locus
pdf("Smokig_prediction.pdf",h=6,w=6)
par(mfrow=c(1,1), mar=c(8,4,4,4),family="NimbusSan",las=2)
stripchart(cg05575921~pData(fun)$smoking_status,
             method="jitter",#group.names=c("Current","Former","Never","Passive"),
             pch=16,cex=1,col=c(4,2,3,1),ylab= "Normalized Beta values",
           main="cg05575921", vertical=TRUE,cex.axis=1.2,cex.lab=1)
boxplot(cg05575921,ylab= "Normalized Beta values", main="cg05575921 (smoking related)")
dev.off()

### If required: Remove duplicates
pData2 <- pData[order(pData(fun)$sample_id), ]
dim(pData2)
head(pData2)

dups = pData2$sample_id[duplicated(pData2$sample_id)]
if( length(dups)>0 ){
    whdups     = lapply(dups, function(x) which(pData2$sample_id==x) )
    whdups2rem = sapply( 1:length(dups) , function(i) rbinom( 1 , 1 ,  prob = 1/length(whdups[[1]])  ) )+1
    torem = sapply( 1:length(whdups) , function(i) whdups[[i]][whdups2rem[i]] )
    pData2 <- pData2[-torem,]
}
fun <- fun[, rownames(pData2)]

# K) Remove poor performing probes
detP2 = detP[rownames(fun),colnames(fun)]
failed = which(rowSums(detP2 < 0.01) != ncol(fun))
fun1 <- fun[ -failed, ]
print(fun1)
# Check that the correct probes have been removed - so no probe names in failed should be in fun1
failedDF <- as.data.frame(failed)
check <- fun1[featureNames(fun1) %in% rownames(failedDF), ]
if(sum(rownames(check)>0)>0) print(paste("Warning: poor performing probes remain in normalised data"))

size.fun1 <- as.data.frame(t(dim(fun1)))

save(detP, detP2, file="PdetectionTables.RData")

# L) Remove cross-reactive probes
if(!is.null(opt$crossreac)){
    print("Remove cross-reactive probes")
  Cross_reactive <- read.csv(opt$crossreac,header=F)$V1
  fun1 <- fun1[ ! featureNames(fun1) %in% Cross_reactive, ] 
  print(fun1)
}else{
    print("No file with cross-reactive probes supplied; to remove cross-reactive probes, use the -c option")
}

size.fun1a <- as.data.frame(t(dim(fun1)))

save(fun1, file="Fun1.RData")

# M) Remove XY probes
autosomes = !(rownames(fun1) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
fun2 = fun1[autosomes,]
print(fun2)

size.fun2 <- as.data.frame(t(dim(fun2)))

save(fun2, file="Fun2.RData")

# N) Remove SNP-containing probes (mandatory)
print("Remove SNP-associated probes")
fun3 <- dropLociWithSnps(fun2, snps=c("SBE", "CpG"), maf = 0.05) 
print(fun3)

size.fun3 <- as.data.frame(t(dim(fun3)))

save(fun3, file="Fun3.RData")

### Optional: Remove multimodal CpGs
if(as.logical(opt$multimodal_filter)){
    print("Remove multimodal beta-values probes")
    nmode<- nmode.mc(minfi::getBeta(fun3), minN = 3, modedist=0.2, nCores = 1)
    nmode.hi <- nmode[nmode>2]
    fun3 <- fun3[ ! featureNames(fun3) %in% names(nmode.hi), ]

    pData(fun3) <- pData2
}else{
     print("Do not remove multimodal beta-values probes")
}

# O) Create analysis tables 
betaNorm <- getBeta(fun3) 
mNorm <- getM(fun3) 

size.betaNorm <- as.data.frame(t(dim(betaNorm)))
size.mNorm <- as.data.frame(t(dim(mNorm)))

# Check for NA and -Inf values
betaRaw.na <- betaRaw[!complete.cases(betaRaw),]
dim(betaRaw.na) 
NAbetas <- rownames(betaRaw.na)
betaNAsNorm <- betaNorm[rownames(betaNorm) %in% NAbetas,] # Should be [1] 0 x ncol(betaNorm)

# Check for infinite values in m table, replace in the no infinite values m table
TestInf <- which(apply(mNorm,1,function(i) sum(is.infinite(i)))>0)
save(TestInf, file="InfiniteValueProbes.RData")
mNoInf <- mNorm
mNoInf[!is.finite(mNoInf)] <- min(mNoInf[is.finite(mNoInf)])
TestInf2 <- which(apply(mNoInf,1,function(i) sum(is.infinite(i)))>0)
TestInf2 # Should be named integer(0)

size.mNoInf <- as.data.frame(t(dim(mNoInf)))

# Write tables
write.table(targets, file="TargetsFile.csv", sep=",", col.names=NA)
write.table(betaNorm, file="NormalisedFilteredBetaTable.csv", sep=",", col.names=NA)
write.table(mNorm, file="NormalisedFilteredMTable.csv", sep=",", col.names=NA)
write.table(mNoInf, file="NormalisedFilteredMTable_noInf.csv", sep=",", col.names=NA)

# Create dimensions table
DimensionsTable <- rbind(size.RGSet, size.MSet, size.gset, size.fun, size.fun1, size.fun1a, size.fun2, 
                         size.fun3, size.betaNorm, size.mNorm, size.mNoInf)
rownames(DimensionsTable) <- c("RGset_size", "MSet_size", "gset_size", "fun_size", "detP_lost", 
                               "CrossRx_lost", "XYChr_lost", "SNP_loss", "beta_size", "m_size", "mNoInf_size")
colnames(DimensionsTable) <- c("Probe_number", "Sample_number")
write.table(DimensionsTable, file="DimensionsTable.txt", sep="\t", col.names=NA)

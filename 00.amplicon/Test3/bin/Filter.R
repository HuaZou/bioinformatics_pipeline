#!/usr/bin/R 

# library 
library(dada2)
library(dplyr)
library(argparser)

# parameter input
parser <- arg_parser("filter trimmed fq by dada2") %>%
    add_argument("-f", "--fwd", help = "trimmed fq1 directroy") %>%
    add_argument("-r", "--rev", help = "trimmed fq2 directroy") %>%
    add_argument("-t1", "--trun1", help = "truncLen 1st parameter") %>%
    add_argument("-t2", "--trun2", help = "truncLen 2nd parameter") %>%
    add_argument("-me", "--maxEE", help = "maxEE") %>%
    add_argument("-tq", "--truncQ", help = "truncQ") %>%
    add_argument("-mn", "--maxN", help = "maxN") %>%
    add_argument("-o", "--out", help = "output directory of filtered fq")
    
args <- parse_args(parser)

# prepare for function 
dir_fq1 <- args$f
dir_fq2 <- args$r
trun1   <- args$t1
trun2   <- args$t2
maxEE   <- args$me
trunQ   <- args$tq
maxN    <- args$mn
out     <- args$o

# output directory of filtered fq.gz
filtpathF <- file.path(out, "FWD")
filtpathR <- file.path(out, "REV") 

# trimmed fq.gz
fastqFs <- sort(list.files(dir_fq1))
fastqRs <- sort(list.files(dir_fq2))
    
if(length(fastqFs) != length(fastqRs)) {
    stop("Forward and reverse files do not match.")
}

# filter trimmed.fq.gz
filter_res <- filterAndTrim(fwd=file.path(dir_fq1, fastqFs),
                    filt=file.path(filtpathF, fastqFs),
                    rev=file.path(dir_fq2, fastqRs),
                    filt.rev=file.path(filtpathR, fastqRs),
                    truncLen=c(trun1, trun2),
                    maxEE=maxEE,
                    truncQ=trunQ,
                    maxN=maxN,
                    rm.phix=TRUE,
                    compress=TRUE,
                    verbose=TRUE,
                    multithread=TRUE)
save(filter_res, paste0(out, "/filter_result.RData"))

# generate picture file
plotQP <- function(fq_dir, fq, name){
  prefix <- paste(out, paste0(name, ".pdf"), sep="/")
  plotQualityProfile(file.path(fq_dir, fq))
  ggsave(prefix)
}

plotQP(fq_dir1, fastqFs, "trim_FWD")
plotQP(fq_dir2, fastqRs, "trim_REV")
plotQP(filtpathF, fastqFs, "filt_FWD")
plotQP(filtpathR, fastqRs, "filt_REV")

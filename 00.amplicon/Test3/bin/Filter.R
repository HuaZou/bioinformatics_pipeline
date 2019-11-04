#!/usr/bin/R 

# library 
library(dada2)
library(dplyr)
library(argparser)
library(ggplot2)

# parameter input
parser <- arg_parser("filter trimmed fq by dada2") %>%
    add_argument("-f", "--fwd", help = "trimmed fq1 directroy") %>%
    add_argument("-r", "--rev", help = "trimmed fq2 directroy") %>%
    add_argument("-t1", "--trun1", type = "integer", help = "truncLen 1st parameter") %>%
    add_argument("-t2", "--trun2", type = "integer", help = "truncLen 2nd parameter") %>%
    add_argument("-me", "--maxEE", type = "integer", help = "maxEE") %>%
    add_argument("-tq", "--truncQ", type = "integer", help = "truncQ") %>%
    add_argument("-mn", "--maxN", type = "integer", help = "maxN") %>%
    add_argument("-o", "--out", help = "output directory of filtered fq")
    
args <- parse_args(parser)

# prepare for function 
dir_fq1 <- args$f
dir_fq2 <- args$r
trun1   <- as.numeric(args$t1)
trun2   <- as.numeric(args$t2)
maxEE   <- as.numeric(args$me)
trunQ   <- as.numeric(args$tq)
maxN    <- as.numeric(args$mn)
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

sample.names <- sapply(strsplit(basename(fastqFs), "_"), `[`, 1)
filtqFs <- paste0(sample.names, "_F_filt.fastq.gz")
filtFs <- file.path(filtpathF, filtqFs)
filtqRs <- paste0(sample.names, "_R_filt.fastq.gz")
filtRs <- file.path(filtpathF, filtqRs)

# filter trimmed.fq.gz
filter_res <- filterAndTrim(
					fwd         = file.path(dir_fq1, fastqFs),
                    filt        = filtFs,
                    rev         = file.path(dir_fq2, fastqRs),
                    filt.rev    = filtRs,
                    truncLen    = c(trun1, trun2),
                    maxEE       = rep(maxEE, 2),
                    truncQ      = trunQ,
                    maxN        = maxN,
                    rm.phix     = TRUE,
                    compress    = TRUE,
                    verbose     = TRUE,
                    multithread = TRUE)

save(filter_res, file=paste0(out, "/filter_result.RData"))

# generate picture file
plotQP <- function(fq_dir, fq, name){
  prefix <- paste(out, paste0(name, ".pdf"), sep = "/")
  plotQualityProfile(file.path(fq_dir, fq))
  ggsave(prefix)
}

plotQP(dir_fq1, fastqFs, "trim_FWD")
plotQP(dir_fq2, fastqRs, "trim_REV")
plotQP(filtpathF, filtqFs, "filt_FWD")
plotQP(filtpathR, filtqRs, "filt_REV")

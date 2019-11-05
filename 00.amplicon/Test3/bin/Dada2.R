#!/usr/bin/R 

library(dada2)
library(dplyr)
library(argparser)
library(ggplot2)

# parameter input
parser <- arg_parser("filter trimmed fq by dada2") %>%
    add_argument("-f", "--fwd", help = "filtered fq1 directroy") %>%
    add_argument("-r", "--rev", help = "filtered fq2 directroy") %>%
    add_argument("-o1", "--out1", help = "output directory of dada2") %>%
    add_argument("-d", "--database", help = "database for 16s rDNA") %>%
    add_argument("-o2", "--out2", help = "output directory of taxonomy") 

args <- parse_args(parser)

# prepare for function 
dir_fq1  <- args$f
dir_fq2  <- args$r
out1     <- args$o1
database <- args$d
out2     <- args$o2

# filtered.fq.gz
filtFs <- file.path(dir_fq1, list.files(dir_fq1))
filtRs <- file.path(dir_fq2, list.files(dir_fq2))

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plot_errF <- paste0(out1, "/plotErrorF.pdf")
plotErrors(errF, nominalQ=TRUE)
ggsave(plot_errF)

plot_errR <- paste0(out1, "/plotErrorR.pdf")
plotErrors(errR, nominalQ=TRUE)
ggsave(plot_errR)

# Sample Inference
derepF <- derepFastq(filtFs)
dadaFs <- dada(derepF, err=errF, multithread=TRUE)
derepR <- derepFastq(filtRs)
dadaRs <- dada(derepR, err=errR, multithread=TRUE)

# merge paired reads
mergers <- mergePairs(dadaFs, derepF, dadaRs, derepR, verbose=TRUE)
save(derepF, dadaFs, derepR, dadaRs, mergers, file=paste0(out1, "/dada2.RData"))

# Construct sequence table
seqtab <- makeSequenceTable(mergers) 
saveRDS(seqtab, paste0(out1, "/seqtab.rds"))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab.nochim, paste0(out1, "/seqtab.nochimera.rds"))

# assign taxonmy
tax <- assignTaxonomy(seqtab.nochim, database, multithread=TRUE)
saveRDS(tax, paste0(out2, "/taxonomy.rds"))


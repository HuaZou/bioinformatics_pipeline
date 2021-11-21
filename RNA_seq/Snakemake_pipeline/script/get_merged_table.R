#!/bin/usr/R 
library(dplyr)
library(argparser)
library(tximport)

# parameter input
parser <- arg_parser("merge transcript abundance files") %>%
    add_argument("-s", "--sample", 
        help = "the sampleID files") %>%
    add_argument("-d", "--dir", 
        help = "the dir of output") %>%
    add_argument("-g", "--gene",
        help = "transcriptID into geneID") %>%
    add_argument("-t", "--type",
        help = "the type of software") %>%
    add_argument("-c", "--count",
        help = "the method of scale for featurecounts") %>%
    add_argument("-n", "--name",
        help = "the name of transcript file") %>%
    add_argument("-o", "--out", 
        help = "the merged files's dir and name")

args <- parse_args(parser)

# prepare for function 
sample <- args$s	 
dir    <- args$d
tx2gen <- args$g
type   <- args$t
scount <- args$c
name   <- args$n
prefix <- args$o

tx2gene_table <- read.csv(tx2gen, header=T)
tx2gene <- tx2gene_table %>% 
            #dplyr::filter(tx_biotype == "protein_coding") %>%
            dplyr::select(tx_id, gene_id, gene_name)
phen <- read.table((sample, header=T, sep="\t")
sample_name <- unique(as.character(phen$SampleID))
files <- file.path(dir, sample_name, name)
names(files) <- sample_name
result <- tximport(files, 
                   type = type, 
                   countsFromAbundance = scount,
                   tx2gene = tx2gene, 
                   ignoreTxVersion = T)
save(result, file = prefix)

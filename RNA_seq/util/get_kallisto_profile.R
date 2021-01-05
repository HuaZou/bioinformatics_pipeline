#!/bin/usr/R
library(dplyr)
library(argparser)
library(tximport)
library(tibble)

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
    add_argument("-c", "--count", default = "no",
        help = "the method of scale for featurecounts") %>%
    add_argument("-p", "--occurrence", default = 0.2,
        help = "the threshold of occurrence") %>%
    add_argument("-nu", "--number", default = 10,
        help = "the threshold of ncount")  %>%
    add_argument("-n", "--name", default = "abundance.h5", 
        help = "the name of transcript file") %>%
    add_argument("-o", "--out",
        help = "the merged files's dir and name")
args <- parse_args(parser)


transform_profile <- function(matrix, 
                        exon_length,
                        y=phen,
                        occurrence=0.2, 
                        ncount=10){
  # matrix=count_format
  # exon_length=count_length 
  # occurrence=0.2
  # ncount=10
  
  # filter with occurrence and ncount
  prf <- matrix %>% rownames_to_column("Type") %>% 
    filter(apply(dplyr::select(., -one_of("Type")), 1, 
                 function(x){sum(x > 0)/length(x)}) > occurrence) %>%
            data.frame(.) %>% 
    column_to_rownames("Type")
  prf <- prf[rowSums(prf) > ncount, ]
  
  # change sampleID 
  sid <- intersect(colnames(prf), y$SampleID)
  phen.cln <- y %>% filter(SampleID%in%sid) %>%
    arrange(SampleID)
  prf.cln <- prf %>% dplyr::select(as.character(phen.cln$SampleID))
  # determine the right order between profile and phenotype 
  for(i in 1:ncol(prf.cln)){ 
    if (!(colnames(prf.cln)[i] == phen.cln$SampleID[i])) {
      stop(paste0(i, " Wrong"))
    }
  }  
  
  # normalizate the profile using FPKM and TPM
  exon_length.cln <- exon_length[as.character(rownames(prf.cln)), , F]
  for(i in 1:nrow(exon_length.cln)){ 
    if (!(rownames(exon_length.cln)[i] == rownames(prf.cln)[i])) {
      stop(paste0(i, " Wrong"))
    }
  } 

  dat <- as.matrix(prf.cln) / exon_length.cln$Length
  dat_FPKM <- t(t(dat)/colSums(prf.cln)) * 10^9
  dat_TPM <- t(t(dat)/colSums(dat)) * 10^6
  
  dat_count <- prf.cln[order(rownames(prf.cln)), ]
  dat_fpkm <- dat_FPKM[order(rownames(dat_FPKM)), ]
  dat_tpm <- dat_TPM[order(rownames(dat_TPM)), ]
  
  res <- list(count=dat_count, fpkm=dat_fpkm, tpm=dat_tpm)
  return(res)
}

output_profile <- function(result, tag="STAR"){
  
  if(!dir.exists(out)){
    dir.create(out)
  }
  count_name <- paste0(out, tag, "_filtered_counts.tsv")
  FPKM_name <- paste0(out, tag, "_filtered_FPKM.tsv")
  TPM_name <- paste0(out, tag, "_filtered_TPM.tsv")  
  
  write.table(x = result$count, 
              file = count_name, 
              sep = '\t', 
              quote = F,
              col.names = NA)
  write.table(x = result$fpkm, 
              file = FPKM_name, 
              sep = '\t', 
              quote = F,
              col.names = NA)
  write.table(x = result$tpm, 
              file = TPM_name, 
              sep = '\t', 
              quote = F,
              col.names = NA)
}


# prepare for function
sample <- args$s
dir    <- args$d
tx2gen <- args$g
type   <- args$t
scount <- args$c
occur  <- args$p
num    <- args$nu
name   <- args$n
out    <- args$o

tx2gene_table <- read.table(tx2gen, header = T, sep = "\t")
tx2gene <- tx2gene_table %>%
            dplyr::select(ensembl_transcript_id, ensembl_gene_id, external_gene_name) 
phen <- read.table(sample, header=T, sep="\t")
sample_name <- unique(as.character(phen$SampleID))
files <- file.path(dir, sample_name, name)
names(files) <- sample_name
result <- tximport(files,
                   type = type,
                   countsFromAbundance = scount,
                   tx2gene = tx2gene,
                   ignoreTxVersion = T)

count_format <- round(result$counts) %>% data.frame()
count_length <- apply(result$length, 1, mean) %>% data.frame() %>%
  setNames("Length")

# calculate values  
final_profile <- transform_profile(count_format, count_length, occurrence=occur, ncount=num)

# output
output_profile(final_profile, tag=type)
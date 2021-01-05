#!/usr/bin/R

# clear all vectors
rm(list = ls())

# library function
library(argparser)
library(dplyr)
library(stringr)
library(data.table)
library(tibble)

# parameter input
parser <- arg_parser("get profile") %>%
    add_argument("-f", "--prof", 
        help = "profile table rownames->taxonomy; colnames->sampleID") %>%
    add_argument("-c", "--occurrence", 
        help = "the threshold of occurrence") %>%
    add_argument("-n", "--number", 
        help = "the threshold of ncount")  %>%
    add_argument("-o", "--out", 
        help = "result with director", default = "./")

args <- parse_args(parser)


transform_profile <- function(matrix, 
                        exon_length,
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
  sid <- intersect(rownames(prf), rownames(exon_length))
  prf.cln <- prf %>% rownames_to_column("geneid") %>%
    filter(geneid%in%sid) %>% arrange(geneid) %>%
    column_to_rownames("geneid")
  
  exon_length.cln <- exon_length %>% rownames_to_column("geneid") %>%
    filter(geneid%in%sid) %>% arrange(geneid) %>%
    column_to_rownames("geneid")

  # determine the right order
  for(i in 1:nrow(prf.cln)){ 
    if (!(rownames(prf.cln)[i] == rownames(exon_length.cln)[i])) {
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

output_profile <- function(result, type="STAR"){
  
  # result <- star_prf
  # type <- "STAR"
  
  if(!dir.exists(out)){
    dir.create(out)
  }
  count_name <- paste0(out, type, "_filtered_counts.tsv")
  FPKM_name <- paste0(out, type, "_filtered_FPKM.tsv")
  TPM_name <- paste0(out, type, "_filtered_TPM.tsv")  
  
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
prof <- fread(args$f)    
currence <- args$c	 
number <- args$n	 
out <- args$o	

# curation data 
count_format <- prof %>% dplyr::select(c("Geneid", ends_with("bam"))) %>%
      dplyr::rename_at(vars(ends_with("bam")), 
                       funs(str_replace(., "(\\S+align_star/)", ""))) %>%
      dplyr::rename_at(vars(ends_with("bam")), 
                       funs(str_replace(., "Aligned.sortedByCoord.out.bam", ""))) %>%
    column_to_rownames("Geneid")
count_length <- prof %>% dplyr::select(Geneid, Length) %>%
  column_to_rownames("Geneid")

# calculate values  
final_profile <- transform_profile(count_format, count_length, occurrence=currence, ncount=number)

# output
output_profile(final_profile)

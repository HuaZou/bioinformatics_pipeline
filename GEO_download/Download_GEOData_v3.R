suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(convert))

option_list <- list(
  make_option(c("-g", "--GEO"), type="character", default="GSE65858",
              help="GEO number", metavar="character"),
  make_option(c("-p", "--GPL"), type="character",
              help="GPL platform number", metavar="character"),
  make_option(c("-t", "--type"), type="character",
              help="Type of GPL platform", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output", metavar="character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
current_dir <- getwd()


# input parameters
GEO_name <- opt$GEO
GPL_number <- opt$GPL
Array_type <- opt$type
dir <- opt$out


# clinical and expression profile
gset <- getGEO(GEO = GEO_name,
               destdir = dir,
               AnnotGPL = F,
               getGPL = F)
phen <- pData(gset[[1]])
prof <- exprs(gset[[1]])
#print(head(phen))
#print(head(prof))

# probe 2 geneid
gpl <- getGEO(GEO = GPL_number, destdir = dir)
#print(head(Table(gpl)[, 1:10]))
if (Array_type == "array") {

  temp_tab <- Table(gpl)
  # ILMN
  colnames(temp_tab)[which(colnames(temp_tab) == "ILMN_Gene")] <- "Gene"

  # microarray
  colnames(temp_tab)[which(colnames(temp_tab) == "Gene Symbol")] <- "Gene"
  colnames(temp_tab)[which(colnames(temp_tab) == "GENE_SYMBOL")] <- "Gene"

  # circRNA
  colnames(temp_tab)[which(colnames(temp_tab) == "circRNA")] <- "Gene"

  # SOMAscan human proteomic assay
  colnames(temp_tab)[which(colnames(temp_tab) == "EntrezGeneSymbol")] <- "Gene"

  # GB_LIST
  #colnames(temp_tab)[which(colnames(temp_tab) == "GB_LIST")] <- "Gene"

  probe2gene <- temp_tab %>%
    dplyr::select(all_of(c("ID", "Gene"))) %>%
    stats::setNames(c("ProbeID", "Gene")) %>%
    dplyr::filter(Gene != "") %>%
    dplyr::mutate(Gene = gsub(" /// \\S+", "", Gene))
} else if (Array_type == "other") {
  library(idmap1)
  probe2gene <- getIDs(tolower(GPL_number)) %>%
    dplyr::select(c("probe_id", "symbol")) %>%
    stats::setNames(c("ProbeID", "Gene"))
} else if (Array_type == "unname") {
  probe2gene <- Table(gpl) %>%
    dplyr::select(c("ID")) %>%
    stats::setNames(c("ProbeID")) %>%
    dplyr::mutate(Gene = ProbeID)
}

##################################################################
# Creating Clinical Annotation Table
preprocess_Clinical <- function(dataset = phen) {

  # dataset = phen

  if (length(grep(":ch1$", colnames(dataset))) > 0) {
    dat <- dataset[, c(1, grep(":ch1$", colnames(dataset)))]
    colnames(dat) <- gsub(":ch1$", "", colnames(dat))
  } else {
    dat <- dataset
  }

  res <- dat %>%
    tibble::rownames_to_column("bcr_patient_barcode") %>%
    dplyr::distinct(bcr_patient_barcode, .keep_all = TRUE) %>%
    reshape::rename(c(bcr_patient_barcode = "Barcode"))

  return(res)
}

# gene expression profile
preprocess_profile <- function(dataset,
                               metadata,
                               annotation = probe2gene) {

  # dataset=prof
  # metadata=phen_post
  # annotation=probe2gene

  sid <- dplyr::intersect(colnames(dataset), metadata$Barcode)
  prof_cln <- data.frame(dataset) %>%
    dplyr::select(dplyr::all_of(sid)) %>%
    tibble::rownames_to_column("ProbeID") %>%
    dplyr::inner_join(annotation, by = "ProbeID")

  # filter features according gene symbol
  idx <- grep("Gene", colnames(prof_cln))
  prof_cln$median <- apply(prof_cln[, -idx], 1, median)
  prof_cln <- with(prof_cln, prof_cln[order(Gene, median, decreasing = T), ])
  prf_deduplicated <- prof_cln[!duplicated(prof_cln$Gene), ] %>%
    dplyr::select(-median)

  # return result: gene symbol
  rownames(prf_deduplicated) <- NULL
  prof_res <- prf_deduplicated %>%
    dplyr::select(ProbeID, Gene, dplyr::everything()) %>%
    dplyr::select(-ProbeID) %>%
    tibble::column_to_rownames("Gene")

  return(prof_res)
}

# ExpressionSet Object
get_ExprSet <- function(metadata = phen_post,
                        profile = prof_post,
                        occurrence = 0.2) {

  # metadata=phen_post
  # profile=prof_post
  # occurrence=0.2

  sid <- intersect(metadata$Barcode, colnames(profile))
  phen_cln <- metadata %>% filter(Barcode %in% sid) %>%
    tibble::column_to_rownames("Barcode")

  prof_cln <- profile %>% 
    dplyr::select(rownames(phen_cln)) %>%
    tibble::rownames_to_column("tmp") %>%
    dplyr::filter(apply(dplyr::select(., -one_of("tmp")), 1, function(x) {
      sum(x != 0) / length(x)}) > occurrence) %>%
    tibble::column_to_rownames("tmp")

  # determine the right order between profile and phenotype
  for (i in 1:ncol(prof_cln)) {
    if (!(colnames(prof_cln)[i] == rownames(phen_cln)[i])) {
      stop(paste0(i, " Wrong"))
    }
  }

  exprs <- as.matrix(prof_cln)
  adf <-  new("AnnotatedDataFrame", data=phen_cln)
  experimentData <- new("MIAME",
                        name="Hua Zou", lab="Hua",
                        contact="zouhua1@outlook.com",
                        title="GEO Expression Data",
                        abstract="The gene ExpressionSet",
                        url="www.zouhua.top",
                        other=list(notes="Created from text files"))
  expressionSet <- new("ExpressionSet", exprs=exprs,
                       phenoData=adf,
                       experimentData=experimentData)

  return(expressionSet)
}
##################################################################

phen_post <- preprocess_Clinical(dataset = phen)
prof_post <- preprocess_profile(dataset = prof,
                                metadata = phen_post,
                                annotation = probe2gene)
ExprSet_object <- get_ExprSet(metadata = phen_post,
                              profile = prof_post,
                              occurrence = 0.1)


# output
outdir <- paste0(dir, "/", GEO_name, "_process/")
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

phen_origin <- paste0(outdir, GEO_name, "_clinical_origin.csv")
phen_process <- paste0(outdir, GEO_name, "_clinical_post.csv")
write.csv(phen, file = phen_origin, row.names = F)
write.csv(phen_post, file = phen_process, row.names = F)

prof_origin <- paste0(outdir, GEO_name, "_profile_origin.tsv")
prof_process <- paste0(outdir, GEO_name, "_profile_post.tsv")
write.table(data.frame(prof) %>% rownames_to_column("GeneID"),
            file = prof_origin, row.names = F, quote = F, sep = "\t")
write.table(prof_post %>% rownames_to_column("GeneID"),
            file = prof_process, row.names = F, quote = F, sep = "\t")

probe2gene_name <- paste0(outdir, GPL_number, "_probe2gene_table.tsv")
write.table(probe2gene, file = probe2gene_name,
            row.names = F, quote = F, sep = "\t")

ExprSet_name <- paste0(outdir, GEO_name, "_GeneExprSet.RDS")
saveRDS(ExprSet_object, file = ExprSet_name)

message("Congrats, Program Ended without problems")

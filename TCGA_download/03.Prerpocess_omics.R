options(warn = -1)
library(dplyr)
library(tibble)
library(data.table)
library(SummarizedExperiment)
library(optparse)
options(warn = 0)

option_list = list(
  make_option(c("-e", "--expr"), type="character",  
        help="SummarizedExperiment Object", metavar="character"),
  make_option(c("-t", "--type"), type="character",  
        help="omic-data type", metavar="character"),
  make_option(c("-g", "--geneanno"), type="character", 
              default = "/disk/user/zouhua/pipeline/TCGA_download_v2/util/human_gene_all.tsv",  
              help="omic-data type", metavar="character"),  
  make_option(c("-p", "--prefix"), type="character", 
        help="prefix", metavar="character")
); 
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

current_dir <- getwd()

# input parameters
ExprObject <- opt$expr
Type <- opt$type
GeneAnno <- opt$geneanno
Prefix <- opt$prefix

# importing data 
SumExperData <- readRDS(ExprObject)

####################################################  
#            phenotypic information                #
#################################################### 
get_phenotype <- function(dataset=SumExperData){
  
  # dataset=SumExperData
  
  clinical_trait <- colData(dataset) %>%
    data.frame() %>%
    dplyr::select(bcr_patient_barcode,
                  gender,
                  race,
                  age_at_index,                
                  vital_status,
                  days_to_death,
                  days_to_last_follow_up,
                  tissue_or_organ_of_origin) %>%
    rownames_to_column("SampleID") %>%
    mutate(Group=ifelse(as.numeric(substr(SampleID, 14, 15)) == 11, "Normal", "Tumor")) %>%
    mutate(bcr_patient_barcode=gsub("\\-", "_", bcr_patient_barcode),
           SampleID=gsub("\\-", "_", SampleID))
  
  # dead people
  dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == "Dead") %>%
    dplyr::select(-days_to_last_follow_up) %>%
    reshape::rename(c(bcr_patient_barcode = "Barcode",
                      gender              = "Gender",
                      race                = "Race",                    
                      age_at_index        = "Age",                    
                      vital_status        = "OS",
                      days_to_death       = "OS.Time")) %>%
    mutate(OS=ifelse(OS == "Dead", 1, 0))%>%
    mutate(OS.Time = OS.Time / 365)
  #alive people
  alive_patient <- clinical_trait %>%
    dplyr::filter(vital_status == "Alive") %>%
    dplyr::select(-days_to_death) %>%
    reshape::rename(c(bcr_patient_barcode = "Barcode",
                      gender              = "Gender",
                      race                = "Race",                    
                      age_at_index        = "Age",                    
                      vital_status        = "OS",
                      days_to_last_follow_up    = "OS.Time")) %>%
    mutate(OS=ifelse(OS == "Dead", 1, 0))%>%
    mutate(OS.Time = OS.Time / 365)  

  post_clinical <- rbind(dead_patient, alive_patient)      
  
  return(post_clinical)
  
}

####################################################  
#            phenotypic information                #
#################################################### 
get_Omics <- function(dataset=SumExperData,
                      datatype=Type){
  
  dataset=SumExperData
  datatype=Type
  
  profile <- assay(dataset) %>%
    data.frame()
  colnames(profile) <- gsub("\\.", "_", colnames(profile))
  
  if(datatype == "mRNA"){
    gene_anno <- fread(GeneAnno)
    # Choosing mRNA data 
    prof <- profile %>%
      rownames_to_column("GeneID") 
    
    prof_anno <- gene_anno %>%
      dplyr::filter(transcript_biotype == "protein_coding") %>%
      dplyr::select(c(external_gene_name, ensembl_gene_id, transcript_biotype)) %>%
      dplyr::inner_join(prof, by=c("ensembl_gene_id"="GeneID")) %>%
      distinct() %>%
      dplyr::select(-c(ensembl_gene_id, transcript_biotype)) 
    
    # filter features according gene symbol
    prof_anno$median <- apply(prof_anno[, -1], 1, median)
    prof_anno <- with(prof_anno, prof_anno[order(external_gene_name, median, decreasing = T), ]) 
    prof_uniq <- prof_anno[!duplicated(prof_anno$external_gene_name), ] %>% 
      dplyr::select(-median) %>%
      column_to_rownames("external_gene_name")     
  }else if(datatype == "DNA_Methylation"){
    prof_uniq <- profile %>% na.omit()
  }  
  
  res <- prof_uniq %>% 
    rownames_to_column("Feature")
  return(res)
}

####################################################  
#            phenotypic information                #
#################################################### 
get_ExprSet <- function(metadata=phen,
                        profile=prof){
  
  # metadata=phen
  # profile=prof
  
  sid <- intersect(metadata$SampleID, colnames(prof))
  phen_cln <- metadata[metadata$SampleID%in%sid, ] %>%
    column_to_rownames("Barcode")
  prof_cln <- prof %>% dplyr::select(phen_cln$SampleID) 
  
  if(!any(colnames(prof_cln) == phen_cln$SampleID)){
    stop("The order of samplenames between phen_cln and prof_cln was wrong")
  }
  
  colnames(prof_cln) <- rownames(phen_cln)
  
  require(convert)
  exprs <- as.matrix(prof_cln)
  adf <-  new("AnnotatedDataFrame", data=phen_cln)
  experimentData <- new("MIAME",
                        name="Hua", lab="Dong gdl Lab",
                        contact="Hua@grmh.cn",
                        title="Experiment",
                        abstract="The ExpressionSet Value",
                        url="www.grmh-gdl.cn",
                        other=list(notes="TCGA Data"))
  expressionSet <- new("ExpressionSet", exprs=exprs,
                       phenoData=adf, 
                       experimentData=experimentData)
  
  return(expressionSet) 
  
}

####################################################  
#             Copy number variation                #
#################################################### 
get_CNV <- function(dataset=SumExperData){
  
  # dataset=SumExperData
  
  hg_marker_file <- read.delim(GeneAnno)
  
  # If you are using Masked Copy Number Segment for GISTIC analysis, 
  # please only keep probesets with freqcnv = FALSE
  hg_marker_file_false <- hg_marker_file %>% 
    filter(freqcnv == FALSE) %>%
    dplyr::select(probeid, chr, pos)
  tumor_seg <- dataset[substr(dataset$Sample, 14, 15) == "01", ]
  untumor_seg <- dataset[substr(dataset$Sample, 14, 15) != "01", ]
  
  # put seg file and marker file into GISTIC software for further analysis
  res <- list(hg=hg_marker_file_false,
              tumor=tumor_seg,
              untumor=untumor_seg)
  return(res)
}


####################################################  
#      extracting  miRNA  profile                  #
####################################################
get_miRNA <- function(dataset=SumExperData){
  
  # dataset=SumExperData
  
  # count 
  count_name <- grep("read_count", colnames(dataset), value = T)
  count_name_new <- gsub("-", "_", gsub("read_count_", "", count_name))
  dat_count <- dataset %>%
    dplyr::select(all_of(count_name))
  rownames(dat_count) <- dataset$miRNA_ID
  colnames(dat_count) <- count_name_new
  dat_count_cln <- dat_count %>% rownames_to_column("Feature")
  
  # F/RPKM
  Per_name <- grep("reads_per_million_miRNA_mapped_", colnames(dataset), value = T)
  Per_name_new <- gsub("-", "_", gsub("reads_per_million_miRNA_mapped_", "", Per_name))
  dat_per <- dataset %>%
    dplyr::select(all_of(Per_name))
  rownames(dat_per) <- dataset$miRNA_ID
  colnames(dat_per) <- Per_name_new
  dat_per_cln <- dat_per %>% rownames_to_column("Feature")
  
  # put seg file and marker file into GISTIC software for further analysis
  res <- list(count=dat_count_cln,
              FPKM=dat_per_cln)
  return(res)
}

####################################################  
#                     process                      #
#################################################### 
 if(!dir.exists(paste0(current_dir, "/Clean/"))){
    dir.create(paste0(current_dir, "/Clean/"))
  }
# mRNA & Methylation
if(Type %in%c("mRNA", "DNA_Methylation")){
  phen <- get_phenotype(dataset=SumExperData)
  phen_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_clinical.csv")
  write.csv(phen, phen_filename, row.names = F)  

  prof <- get_Omics(dataset=SumExperData, datatype=Type)
  prof_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_profile.tsv")
  write.table(prof, prof_filename, row.names = F, quote = F, sep = "\t")
  
  ExprSet <- get_ExprSet(metadata=phen, profile=prof)
  ExprSet_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_ExprSet.RDS")
  saveRDS(ExprSet, ExprSet_filename, compress = TRUE)

  message("ExprSet object of Omics has been successfully transformed")
}

# CNV 
if(Type == "CNV"){
  CNV_list <- get_CNV(dataset=SumExperData)
  marker_filename <- paste0(current_dir, "/Clean/", "snp6.na35.remap.hg38.tsv")
  write.table(CNV_list$hg, marker_filename, row.names = F, quote = F, sep = "\t")
  
  tumor_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_tumor.tsv")
  write.table(CNV_list$tumor, tumor_filename, row.names = F, quote = F, sep = "\t") 
  
  untumor_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_untumor.tsv")
  write.table(CNV_list$untumor, untumor_filename, row.names = F, quote = F, sep = "\t")
  
  message("CNV has been successfully transformed")
  
}

# miRNA
if(Type == "miRNA"){
  miRNA_list <- get_miRNA(dataset=SumExperData)

  count_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_count.tsv")
  write.table(miRNA_list$count, count_filename, row.names = F, quote = F, sep = "\t") 
  
  per_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_PKM.tsv")
  write.table(miRNA_list$FPKM, per_filename, row.names = F, quote = F, sep = "\t")

  message("miRNA has been successfully transformed")
  
}

#  
if(Type == "DNA_Methylation"){
  phen <- get_phenotype(dataset=SumExperData)
  if(!dir.exists(paste0(current_dir, "/Clean/"))){
    dir.create(paste0(current_dir, "/Clean/"))
  }
  phen_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_clinical.csv")
  write.csv(phen, phen_filename, row.names = F)  
  
  prof <- get_Omics(dataset=SumExperData, datatype=Type)
  prof_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_profile.tsv")
  write.table(prof, prof_filename, row.names = F, quote = F, sep = "\t")
  
  ExprSet <- get_ExprSet(metadata=phen, profile=prof)
  ExprSet_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_ExprSet.RDS")
  saveRDS(ExprSet, ExprSet_filename, compress = TRUE)

  message("ExprSet object of Omics has been successfully transformed")
}

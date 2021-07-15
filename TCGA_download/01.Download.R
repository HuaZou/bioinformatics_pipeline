options(warn = -1)
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(optparse)
options(warn = 0)


option_list <- list(
  make_option(c("-t", "--tumor"), type="character", default="TCGA-KIRC", 
        help="tumor type from TCGA", metavar="character"),
  make_option(c("-o", "--out"), type="character",  
        help="output directory", metavar="character")
); 
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

current_dir <- getwd()

# input parameters
tumor_type <- opt$tumor
data_type <- opt$out


# clinical information output
if(!dir.exists(paste0(current_dir, "/Clinical/"))){
  dir.create(paste0(current_dir, "/Clinical/"))
}
origin_clinical_name <- paste0(current_dir, "/Clinical/", tumor_type, "_clinical_origin.csv")
if(!file.exists(origin_clinical_name)){
  clinical <- GDCquery_clinic(project = tumor_type, type = "clinical")
  write.csv(clinical, origin_clinical_name, row.names = F)
}
message("Clinical information has been successfully downloaded")


####################################################  
#   omics-data layer: SummarizedExperiment Object  #
#################################################### 
get_OmicsData <- function(project  = tumor_type,
                          Outdir   = "mRNA"){
  if(Outdir == "mRNA"){
    query_Data <- GDCquery(project = project,
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type = "HTSeq - FPKM") # HTSeq - FPKM-UQ; HTSeq - Counts; STAR - Counts
  }else if(Outdir == "miRNA"){
    query_Data <- GDCquery(project = project,
                           data.category = "Transcriptome Profiling",
                           data.type = "miRNA Expression Quantification",
                           workflow.type = "BCGSC miRNA Profiling")     
  }else if(Outdir == "CNV"){
    query_Data <- GDCquery(project = project,
                           data.category = "Copy Number Variation",
                           data.type = "Copy Number Segment")     
  }else if(Outdir == "DNA_Methylation"){
    query_Data <- GDCquery(project = project,
                           data.category = "DNA methylation",
                           legacy = TRUE,
                           platform = "Illumina Human Methylation 450")     
  }
  
  GDCdownload(query = query_Data,
              method = "api",
              files.per.chunk = 60,
              directory = Outdir)
  
  expdat <- GDCprepare(query = query_Data,
                       directory = Outdir)
  return(expdat)
}

OmicsData <- get_OmicsData(project = tumor_type, Outdir = data_type)
if(!dir.exists(paste0(current_dir, "/Omics/"))){
  dir.create(paste0(current_dir, "/Omics/"))
}
OmicsData_filename <- paste0(current_dir, "/Omics/", tumor_type, "_", data_type, ".RDS")                         
saveRDS(OmicsData, file = OmicsData_filename)
message("Omics data also has been successfully downloaded")

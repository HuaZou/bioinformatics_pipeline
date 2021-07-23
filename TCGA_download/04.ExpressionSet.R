options(warn = -1)
library(dplyr)
library(tibble)
library(data.table)
library(convert)
library(optparse)
options(warn = 0)

option_list = list(
  make_option(c("-p", "--phen"), type="character",  
              help="clinical information", metavar="character"),
  make_option(c("-e", "--expr"), type="character",  
              help="expression data", metavar="character"),  
  make_option(c("-t", "--type"), type="character",  
              help="omic-data type", metavar="character"),
  make_option(c("-o", "--prefix"), type="character", 
              help="prefix", metavar="character")
); 
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

current_dir <- getwd()

# input parameters
Phenotype_para <- opt$phen
Profile_para <- opt$expr
Type <- opt$type
Prefix <- opt$prefix

# importing data
Phenotype <- read.csv(Phenotype_para)
Profile <- fread(Profile_para)


####################################################  
#            phenotypic information                #
#################################################### 
get_ExprSet <- function(metadata=Phenotype,
                        profile=Profile){
  
  metadata=Phenotype
  profile=Profile
  
  colnames(profile) <- substr(colnames(Profile), 1, 12)
  
  sid <- intersect(metadata$Barcode, colnames(profile))
  phen_cln <- metadata %>% filter(Barcode%in%sid) %>%
    column_to_rownames("Barcode")
  prof_cln <- profile %>% column_to_rownames("Feature") %>%
    dplyr::select(rownames(phen_cln)) 
  
  if(!any(colnames(prof_cln) == rownames(phen_cln))){
    stop("The order of samplenames between phen_cln and prof_cln was wrong")
  }
  
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

ExprSet <- get_ExprSet(metadata=Phenotype, profile=Profile)
ExprSet_filename <- paste0(current_dir, "/Clean/", Prefix, "_", Type, "_ExprSet_clinical.RDS")
saveRDS(ExprSet, ExprSet_filename, compress = TRUE)
message("ExpressionSet has been successfully transformed")

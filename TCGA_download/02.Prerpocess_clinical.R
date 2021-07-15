options(warn = -1)
library(dplyr)
library(tibble)
library(optparse)
options(warn = 0)

option_list = list(
  make_option(c("-p", "--phen"), type="character",  
        help="clinical information", metavar="character"),
  make_option(c("-t", "--type"), type="character",  
        help="tumor_type", metavar="character")
); 
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
current_dir <- getwd()

# input parameters
pheno <- opt$phen
tumor_type <- opt$type

# import data 
clinical <- read.csv(pheno)

#################################################
#                    TCGA-PAAD                  #
#################################################
preprocess_TCGA_PAAD <- function(dataset = clinical){

  #dataset = clinical

  clinical_trait <- dataset  %>% 
    dplyr::select(bcr_patient_barcode,
                  gender,
                  race,
                  age_at_index,                
                  vital_status,
                  pack_years_smoked,
                  cigarettes_per_day,                
                  days_to_death,
                  days_to_last_follow_up,
                  primary_diagnosis,
                  tissue_or_organ_of_origin,
                  ajcc_pathologic_stage,             
                  ajcc_pathologic_t,
                  ajcc_pathologic_n,
                  ajcc_pathologic_m) %>%
    distinct(bcr_patient_barcode, .keep_all = TRUE) %>%
    mutate(bcr_patient_barcode=gsub("-", "_", bcr_patient_barcode))
  # dead people
  dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == "Dead") %>%
    dplyr::select(-days_to_last_follow_up) %>%
    reshape::rename(c(bcr_patient_barcode = "Barcode",
                      gender              = "Gender",
                      race                = "Race",                    
                      age_at_index        = "Age",                    
                      vital_status        = "OS",
                      days_to_death       = "OS.Time",
                      tissue_or_organ_of_origin = "Tissue_Origin",
                      ajcc_pathologic_stage = "Stage",
                      ajcc_pathologic_t     = "T",
                      ajcc_pathologic_n     = "N",
                      ajcc_pathologic_m     = "M")) %>%
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
                      days_to_last_follow_up    = "OS.Time",
                      tissue_or_organ_of_origin = "Tissue_Origin",
                      ajcc_pathologic_stage = "Stage",
                      ajcc_pathologic_t     = "T",
                      ajcc_pathologic_n     = "N",
                      ajcc_pathologic_m     = "M")) %>%
    mutate(OS=ifelse(OS == "Dead", 1, 0))%>%
    mutate(OS.Time = OS.Time / 365)

  # merger data 
  post_clinical <- rbind(dead_patient, alive_patient)
  return(post_clinical)  
}


#################################################
#                    TCGA-UCS                   #
#################################################
preprocess_TCGA_UCS <- function(dataset = clinical){
  
  # dataset = clinical
  
  clinical_trait <- dataset  %>% 
    dplyr::select(bcr_patient_barcode,
                  gender,
                  height,
                  weight,
                  bmi,
                  race,
                  ethnicity,
                  age_at_index,                
                  vital_status,
                  cigarettes_per_day,                
                  days_to_death,
                  days_to_last_follow_up,
                  primary_diagnosis,
                  tissue_or_organ_of_origin,
                  figo_stage) %>%
    distinct(bcr_patient_barcode, .keep_all = TRUE) %>%
    mutate(bcr_patient_barcode=gsub("-", "_", bcr_patient_barcode))
  # dead people
  dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == "Dead") %>%
    dplyr::select(-days_to_last_follow_up) %>%
    reshape::rename(c(bcr_patient_barcode = "Barcode",
                      gender              = "Gender",
                      race                = "Race",                    
                      age_at_index        = "Age",                    
                      vital_status        = "OS",
                      days_to_death       = "OS.Time",
                      tissue_or_organ_of_origin = "Tissue_Origin",
                      ajcc_pathologic_stage = "Stage")) %>%
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
                      days_to_last_follow_up    = "OS.Time",
                      tissue_or_organ_of_origin = "Tissue_Origin",
                      ajcc_pathologic_stage = "Stage")) %>%
    mutate(OS=ifelse(OS == "Dead", 1, 0))%>%
    mutate(OS.Time = OS.Time / 365)
  
  # merger data 
  post_clinical <- rbind(dead_patient, alive_patient)
  return(post_clinical)  
}

if(!dir.exists(paste0(current_dir, "/Clinical/"))){
  dir.create(paste0(current_dir, "/Clinical/"))
}
post_clinical_name <- paste0(current_dir, "/Clinical/", tumor_type, "-post_clinical.csv")
if(tumor_type == "TCGA-PAAD"){
  postclinical <- preprocess_TCGA_PAAD(dataset = clinical)
  write.csv(postclinical, post_clinical_name, row.names = F) 
}else if(tumor_type == "TCGA-UCS"){
  postclinical <- preprocess_TCGA_UCS(dataset = clinical)
  write.csv(postclinical, post_clinical_name, row.names = F)  
}

message("Clinical information has been successfully downloaded")

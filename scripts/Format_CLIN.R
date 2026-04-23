## load libraries
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
#input_dir <- args[1]
#output_dir <- args[2]
#annot_dir <- args[3]

input_dir <- "data/input"
output_dir <- "data/output"
annot_dir <- "data/annot"

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" , header=TRUE )
rownames(clin) <- clin$geo_accession
cols <- c('sample_name', 'gender', 'patient_age_at_diagnosis', 'Clinical_Stage',
          'Response', 'OS', 'OS_time', 'PFI_1', 'PFI_time_1')

clin <- clin[, c(cols, colnames(clin)[!colnames(clin) %in% cols])]
colnames(clin)[colnames(clin) %in% cols] <- c('patient', 'sex', 'age', 'stage', "recist",
                                              'os', 't.os', 'pfs', 't.pfs')

new_cols <- c( "primary" , "response.other.info" , "response" , 'histo', "drug_type" , "dna" , "rna", "rna.info")
clin[new_cols] <- NA

clin$drug_type <- 'PD-1/PD-L1'
clin$rna <- 'rnaseq'
clin$rna.info <- 'tpm'
clin$primary <- 'Kidney'
clin$sex <- ifelse(clin$sex == 'male', 'M', 'F')
clin$stage <- toupper(gsub("_clinical_stage", "", clin$stage))

recist_clean <-  sub("_response$", "", clin$recist)
clin$recist <- dplyr::case_when(
  grepl("partial", recist_clean) ~ "PR",
  grepl("stable", recist_clean) ~ "SD",
  grepl("progressive", recist_clean) ~ "PD",
  TRUE ~ NA_character_
)

clin$t.os <- clin$t.os / 30.44
clin$t.pfs <- clin$t.pfs / 30.44

clin$response = Get_Response( data=clin )

# Tissue and drug annotation
annotation_tissue <- read.csv(file=file.path(annot_dir, 'curation_tissue.csv'))
clin <- annotate_tissue(clin=clin, study='Choueiri', annotation_tissue=annotation_tissue, check_histo=FALSE)

annotation_drug <- read.csv(file=file.path(annot_dir, 'curation_drug.csv'))
clin <- add_column(clin, treatmentid=clin$drug_type, .after='tissueid')

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )


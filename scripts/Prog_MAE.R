## load libraries

library(Biobase)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(stringr)

source_location <- "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main"
# source_location <- "~/Documents/GitHub/Pachyderm/PredictIO_ICB/ICB_Common"
load("C:/ICB_Choueiri/data/output/Gencode.v40.annotation.RData")

get_MultiAssayExp <- function(study, input_dir, expr_with_counts_isoforms=FALSE){ 
  source(paste(source_location, "code", "CNA_Functions.R", sep = "/"))
  source(paste(source_location, "code", "SNV_Functions.R", sep = "/"))
  source(paste(source_location, "code", "Create_SummarizedExp.R", sep = "/"))
  data = data.frame(study = study,
                    expr_bool = TRUE,
                    snv_bool = FALSE,
                    cna_bool = FALSE,
                    cin_bool = FALSE,
                    coverage = FALSE,
                    indel_bool = FALSE,
                    mutsig_bool = FALSE,
                    WES.comment = NA)
  data <- data[data$study == study, ]
  
  se_list <- Create_SummarizedExperiments( 
    input_dir=input_dir,
    study= data$study, 
    expr_bool= data$expr_bool, 
    snv_bool= data$snv_bool, 
    cna_bool= data$cna_bool, 
    cin_bool= data$cin_bool, 
    coverage= data$coverage, 
    indel_bool= data$indel_bool,
    expr_with_counts_isoforms=expr_with_counts_isoforms
  )
  
  cols <-list()
  for(assay_name in names(se_list)){
    cols[[assay_name]]<- data.frame(colData(se_list[[assay_name]]))
  }
  
  # Format and merge coldata
  allcols <- lapply(cols, function(col){
    return(rownames(col))
  })
  allcols <- unique(unlist(allcols))
  
  coldata <- NA
  for(col in cols){
    if(!is.data.frame(coldata)){
      coldata <- col
    }
    missing <- allcols[!allcols %in% rownames(coldata)]
    filtered <- col[rownames(col) %in% missing, ]
    if(nrow(filtered) > 0){
      coldata <- rbind(coldata, filtered)
    }
  }
  coldata <- coldata[order(rownames(coldata)), ]
  
  return(MultiAssayExperiment(experiments=se_list, colData=coldata))
}

## get MAE object

dat <- get_MultiAssayExp(study = "Choueiri", output_dir, expr_with_counts_isoforms=FALSE)
saveRDS(dat, file=file.path(output_dir, "ICB_Choueiri.rds"))


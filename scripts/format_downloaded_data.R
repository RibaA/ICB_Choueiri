## load libraries

library(data.table)
library(readxl)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
#work_dir <- args[1]
work_dir <- "data"

# load data from iATLAS 
dat <- readRDS(file.path(work_dir, 'Choueiri_CCR_2016.rds'))

# CLIN.txt
clin <- dat$clin
clin <- as.data.frame(clin)
clin <- clin[order(clin$sample_name), ]

# EXP_TPM.tsv 
expr <- dat$expr
expr <- as.data.frame(expr)
expr <- expr[order(expr$Run_ID), ]
rownames(expr) <- expr$Run_ID
expr <- expr[, -1]
expr <- expr[, order(colnames(expr))]

# samples with both clinical and expression data
int <- intersect(clin$sample_name, rownames(expr))
clin <- clin[clin$sample_name %in% int, ]
expr <- expr[rownames(expr) %in% int, ]

write.table(clin, file=file.path(work_dir, 'CLIN.txt'), sep = "\t" , quote = FALSE , row.names = FALSE)

expr <- t(expr)
expr_updated <- apply(expr, 2, function(x) as.numeric(trimws(x)))
rownames(expr_updated) <- rownames(expr)
expr <- expr_updated

write.table(expr, file=file.path(work_dir, 'EXP_TPM.tsv'), sep = "\t" , quote = FALSE , row.names = TRUE, col.names=TRUE)

## load libraries

library(data.table)
library(R.utils)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
#input_dir <- args[1]
#output_dir <- args[2]

input_dir <- "data/input"
output_dir <- "data/output"

tpm = as.matrix( read.csv( file.path(input_dir, "EXP_TPM.tsv"), stringsAsFactors=FALSE , sep="\t" , header=TRUE ) )
colnames(tpm) <- gsub("\\.", "-", colnames(tpm))

#############################################################################
#############################################################################

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
tpm = log2( tpm[ , colnames(tpm) %in% case[ case$expr %in% 1 , ]$patient ] + 0.001 )

write.table( tpm , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )

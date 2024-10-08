### Code for analyzing RNA seq data using DESeq2 

#load libraries 
library(DESeq2)
library(ggplot2)

setwd("~/projects/eco_genomics/transcriptomics/")

#Import counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)

countsTableRound <- round(countsTable) #bc DESeq2 doesnt like decimals 
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

dss <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds,
                              design = ~ DevTemp + FinalTemp)








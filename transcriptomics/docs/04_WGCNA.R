###
###Script for analyzing and visualizing gene correlation networks 

#check out of library
library(DESeq2)
library(ggplot2)
library(WGCNA); options(stringAsFactors = FALSE);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/transcriptomics/")


##########################################################################

## STEP 1: import our data

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)


countsTableRound <- round(countsTable) #bc DESeq2 doesnt like decimals 
tail(countsTableRound)

#read in conditions table 
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, row.names = 1)

# Filter the matrix to just BASE data bc those are the data for which we have traits measured
filtered_count_matrix_BASEonly <- countsTable[, conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE", ]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)


## STEP 2: Detecting outliers 

#detect outlier genes
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)

table(gsg$goodGenes)
# FALSE  TRUE 
# 37235 82203 

table(gsg$goodSamples)
#TRUE 
#7 

# FIlter out bad genes 
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes == TRUE,]
dim(data_WGCNA) #good job :))

# use clustering with tree dendrogram to identify outlier samples 
htree <- hclust(dist(t(data_WGCNA)), method = 'average')
plot(htree)

# PCA - outlier detection method 
pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x
# make a data frame
pca_data <- as.data.frame(pca_data)

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca_data, aes(PC1, PC2))+
  geom_point()+
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2]))


## STEP 3 - Normalization 

colData <- row.names(filtered_sample_metadata_BASEonly)

#run deseq2 with no matrix defined 
dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sample_metadata_BASEonly,
                                    design = ~1 ) #~ there are no specified groups

dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >=6,]
nrow(dds_WGCNA_75)  #filtered down to 29,559 transcripts 


dds_norm <- vst(dds_WGCNA_75) #perform Variance stableization

#get and save normalized counts to use below 
norm.counts <- assay(dds_norm) %>% 
  t()

### STEP 4: Network construction 

# Choose a set of soft-threshholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# call the network topology analysis function (takes a couple min to run)
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed", #to focus on transcripts that are positively corrilated with 
                         verbose = 5)

sft.data <- sft$fitIndices

#plot to pick power 

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'tomato') +
  labs(x= "Power", y= "Scale free topology model fit, signed R^2" ) +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'tomato') +
  labs(x= "Power", y= "Mean Conectivity" ) +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)


















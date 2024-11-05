
###########################
# Sage Sularz 
# Ecological Genomics
# Transcriptomics Homework
###########################

library(DESeq2)
library(ggplot2)
library(WGCNA); options(stringAsFactors = FALSE);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/transcriptomics/")


##########################################################################

## STEP 1: Import data ##
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)

countsTableRound <- round(countsTable)  
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, row.names = 1)

filtered_count_matrix_BASEonly <- countsTable[, conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE", ]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)


## STEP 2: Detecting outliers ##

gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)

table(gsg$goodGenes)
# FALSE  TRUE 
# 37235 82203 

table(gsg$goodSamples)
#TRUE 
#7 

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


## STEP 3 - Normalization ##

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

## STEP 4: Network construction ##

# Choose a set of soft-threshholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# call the network topology analysis function 
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

###################################################

## Using power thresholds of 26 and 14 compare...
     # the number of modules
     # number of genes per module
     # the strength of correlations with traits.
 
soft_power <- 14
temp_cor <- cor
cor <- WGCNA::cor  

norm.counts[] <- sapply(norm.counts, as.numeric)

# creates the network and identifies modules based on the perameters (14 vs 26) that I chose
bwnet14 <- blockwiseModules(norm.counts,
                            maxBlockSize = 30000,
                            TOMType = "signed",
                            power = soft_power, 
                            mergeCutHeight = 0.25,
                            numericLabels = FALSE,
                            randomSeed = 1234,
                            verbose = 3)


cor <- temp_cor # resets the temp_cor func to base R's cor func. instead of using WCGNA's 

## Step 5: Explore Module ##

module_eigengenes <- bwnet14$MEs

head(module_eigengenes)
dim(module_eigengenes) # Number of Modules

#get the number of genes for each module
table(bwnet14$colors)

# Plot the dendrogram and module colors 
plotDendroAndColors(bwnet14$dendrograms[[1]], cbind(bwnet14$unmergedColors, bwnet14$colors),
                    c("unmerged", "merged"),
                    addGuide = TRUE,
                    hang = 0.03, 
                    guideHang = 0.05)

saveRDS(bwnet14, file = "hmk/bwnet14.rds")

# To load the bwent file in later use:
#bwent26 <- readRDS("outputs/bwent26.rds") 

## Step 6: Correlations of modules with traits ##

#define the numbers of genes and samples 
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Test for correlation between module eigengenes and trait data
module.trait.corr <- cor(module_eigengenes, traitData, use = 'p')

# calculate pval for each correlation 
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Visualize module-trait association as a heat map
heatmap.data <- merge(module_eigengenes, traitData, by ='row.names')
head(heatmap.data)

# adress error of row.names not being numeric 
heatmap.data <- heatmap.data %>%
  column_to_rownames(var = 'Row.names')

names(heatmap.data)

# Make pretty heatmap of correlations 
CorLevelPlot(heatmap.data, 
             x = names(heatmap.data)[42:44],  # these vals will change based on 
             y = names(heatmap.data)[1:41],   # number of eigengenes
             col = c("blue2", "skyblue", "beige", "pink", "tomato"))































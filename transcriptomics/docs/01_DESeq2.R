### Code for analyzing RNA seq data using DESeq2 

#load libraries 
library(DESeq2)
library(ggplot2)

options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/transcriptomics/")

#Import counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)


countsTableRound <- round(countsTable) #bc DESeq2 doesnt like decimals 
tail(countsTableRound)

#read in conditions table 
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)


########################## Explore counts matrix ########################## 


#how many reads do we have from each sample?
colSums(countsTableRound)
mean(colSums(countsTableRound)) #average number of reads per sample #18454529
 #aim for 20 million reads and after filtering maybe a bit lower

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound), 
        cex.names = 0.5, las = 2, ylim = c(0,30000000))
abline(h=mean(colSums(countsTableRound)), col = "blue4", lwd=2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound))  #3244.739
median(rowSums(countsTableRound))  #64 #evidence of over dispersion 

apply(countsTableRound,2,mean) #can see on avg across samples 
                               #gives sense of variation in sequencing effort across samples


######################## Start analysis in DESeq2 ###############################


dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, #
                              design = ~ DevTemp + FinalTemp)
dim(dds)

#start filtering based on avg number of transcripts
dds <- dds[rowSums(counts(dds) >= 10) >= 15,]
nrow(dds) #35,527 #number of transcripts with more than 10 reads in
          #more than or equal to 15

#run the DESeq model to test for global differential gene expression 
dds <- DESeq(dds)

#list the results you have generated
resultsNames(dds) 
#   [1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

# visualize our global gene expression patterns using PCA
#first wee need to transform the data for plotting using variance stabilization 

vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup=c("DevTemp", "FinalTemp"), returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))

final_temp_colors <- c("BASE" = "grey", "A28" = "hotpink", "A33" = "tomato")
shapes_choose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp))+
  geom_point(size = 5) +
  scale_shape_manual(values = shapes_choose) +
  scale_color_manual(values = final_temp_colors) +
  labs(x = paste0('PC1: ', percentVar[1], '%'),
       y = paste0('PC2: ', percentVar[2], '%')) +
  theme_bw(base_size = 16)

p






























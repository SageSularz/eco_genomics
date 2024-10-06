#Ecological Genomics 
#Homework 1 
#Effect of Indv Level Missingness 
#Sage Sularz 


#load libraries 
library(vcfR)
library(SNPfiltR)
library(tidyverse)
library(qqman)

#set WD and call files 
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
head(vcf)
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta") #reference genome 
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="") #annotated genome

#associate files 
chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)
plot(chr1)

#extract and analyis DP 
DP <- extract.gt(vcf, 
                 element="DP", 
                 as.numeric=T, 
                 convertNA=T)
quantile(DP)

DP[DP==0] <- NA 

quantile(DP, na.rm=T) 

dim(DP) 

heatmap.bp(DP[1:1000,], rlabels = F, clabels = F)

#read and prepare metadata 
meta <- read.csv("metadata/meta4vcf.csv", header=T)
head(meta)
meta2 <- meta[,c(1,4)]  #pull ID and Region
head(meta2)

names(meta2) <- c("id", "pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

#Filter min and max dp 
vcf.filt <- hard_filter(vcf, 
                        depth=3) 

vcf.filt <- filter_allele_balance(vcf.filt,
                                  min.ratio = 0.15,
                                  max.ratio=0.85)
vcf.filt <- max_depth(vcf.filt, 
                      maxdepth=60) 

#Filter indv missingness
vcf.filt.indMiss <- missing_by_sample(vcf.filt, 
                                      popmap = meta2,
                                      cutoff=0.75) 
meta75 <- meta[meta$id %in% colnames(vcf.filt.indMiss@gt),]

#vcf.filt.indMiss <- missing_by_sample(vcf.filt, 
#                                      popmap = meta2,
#                                      cutoff=0.5) 
#meta50 <- meta[meta$id %in% colnames(vcf.filt.indMiss@gt),]

#vcf.filt.indMiss <- missing_by_sample(vcf.filt, 
#                                      popmap = meta2,
#                                      cutoff=0.9) 
#meta90 <- meta[meta$id %in% colnames(vcf.filt.indMiss@gt),]

# gets rid of monomrophic or multi-allelic sites
vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss) 

#filtering SNP missingness
vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, 
                                      cutoff=0.5)

vcf.filt.indSNPMiss <- min_mac(vcf.filt.indSNPMiss,
                               min.mac=2)

DP2 <- extract.gt(vcf.filt.indSNPMiss, 
                  element="DP", 
                  as.numeric=T, 
                  convertNA=T)


heatmap.bp(DP2[1:5000,], rlabels=F, clabels=F)

write.vcf(vcf.filt.indSNPMiss, 
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

######## PART TWO === Diversity #######






library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA) 

options(bitmapType = "cairo") #if you need it

setwd("~/projects/eco_genomics/population_genomics/")
#path picks up where WD left off 
vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

#we need to thin the SNPs for LD (linkage disequalibrium) before we run
#PCA and admixture to satify the assumption of independence among loci

vcf.thin <- distance_thin(vcf, min.distance = 500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

#need to subset meta file
meta2 <- meta[meta$id %in% colnames (vcf.thin@gt[, -1]),]

dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

#hide the uncompressed VCF file too big for github outside our repo 

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf") # opens -leaves alone

geno <- vcf2geno(input.file="/gpfs1/home/s/s/ssularz/vcf_final.filtered.thinned.vcf",
                 output.file="outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)

#lots of things just happened in not sure about some of them but we will go back into it next class 

plot(CentPCA$projections,
     col=as.factor(meta2$region))
legend("bottomright", legend = as.factor(unique(meta2$region),
                                         fill=as.factor(unique(meta2$region))))
#somethings not right probs just need to proofread 

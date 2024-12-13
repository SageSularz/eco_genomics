
# Ecological Genomics - Homework 1
# Effect of Individual-Level Missingness - Sage Sularz
library(vcfR)
library(SNPfiltR)

options(bitmapType = "cairo")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")

chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna,ann=gff)

DP <- extract.gt(vcf, element = "DP", as.numeric = T)
DP[DP==0] <- NA

vcf.filt <- hard_filter(vcf, depth = 3) 
vcf.filt <- max_depth(vcf.filt, maxdepth = 60) 

meta <- read.csv("metadata/meta4vcf.csv", header=T)
meta2 <- meta[,c(1,4)]  #pull ID and Region
names(meta2) <- c("id", "pop")
meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

filter_vals <- c(0.75, 0.5, 0.3)

# Loop through the cutoff values
for (i in 1:3) {
  filter_val <- filter_vals[i]
  
  vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                        popmap=meta2,
                                        cutoff=filter_val) # INDV MISSINGNESS
  
  vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
  vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1) #how many times do you need to see a allele to filter it?
  
  vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5) #50% is standard for completness its a good balance 
  
  DP2 <- extract.gt(vcf.filt.indSNPMiss,
                    element = "DP",
                    as.numeric = T)
  
  heatmap.bp(DP2[1:5000, ],
             rlabels = F, clabels = F)
  print("loop")
  
  write.vcf(vcf.filt.indSNPMiss, 
            file=paste0("~/projects/eco_genomics/population_genomics/hmk/vcf_final.filtered_",filter_val,".vcf.gz"))
}




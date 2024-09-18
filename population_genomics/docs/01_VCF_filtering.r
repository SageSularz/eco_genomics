library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

list.files()
list.files("variants/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
vcf #gives quick sum sample #, chromo #, size, etc 
head(vcf) # peak into the organization of the file...meta, fixed,

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")

chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna,ann=gff)

plot(chr1)

pdf(file="~/projects/eco_genomics/population_genomics/figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()

### 09-17-24
DP <- extract.gt(vcf, element = "DP", as.numeric = T)
dim(DP)
DP[1:5,1:10] #how many rows and individuals 

quantile(DP, na.rm = T) #what does the dist. look like?

DP[DP==0] <- NA

#Visualize the matrix of DP and missingness in 

heatmap.bp(DP[1:1000,], rlabels = F, clabels = F)

library(SNPfiltR) #visualize filtering step and apply filters then export vcf file so you dont have to repeat things

vcf.filt <- hard_filter(vcf, depth = 3) #what val do we want to assgin to hard filter?

                                        #could explore more values DP 5, 10....
vcf.filt <- max_depth(vcf.filt, maxdepth = 60) #filter out genos with >60 reads/SNPs

meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2 <- meta[,c(1,4)]  #pull ID and Region
head(meta2)

names(meta2) <- c("id", "pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=0.75) #0 wouldnt filter anyone 0.5 would filter half
vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1) #how many times do you need to see a allele to filter it?

vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5) #50% is standard for completness its a good balance 

DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element = "DP",
                  as.numeric = T)

heatmap.bp(DP2[1:5000, ],
           rlabels = F, clabels = F)

write.vcf(vcf.filt.indSNPMiss, 
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")






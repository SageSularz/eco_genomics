#Ecological Genomics 
#Homework 1 
#Effect of Indv Level Missingness 
#Sage Sularz 


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

vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=0.75) # INDV MISSINGNESS
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

######## PART TWO === Diversity #######
X11.options(type="cairo")

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] #%in% "also found in"
dim(meta2)

# calculate diversity stats using the genetic_diff fxn in vcfR
#function calculates both Gst (a measure of genetic differentiation) 
# and Hs (expected heterozygosity or genetic diversity within populations) for each region
vcf.div <- genetic_diff(vcf,
                        pops=as.factor(meta2$region),
                        method = "nei")
#sd(vcf.div$Hs)
str(vcf.div)

chr.main <- unique(vcf.div$CHROM)[1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))

vcf.div.MHplot <- vcf.div.MHplot %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main,"_",POS))

str(vcf.div.MHplot)

vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)

vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)

manhattan(vcf.div.MHplot,
          chr="V2" ,
          bp="POS",
          p="Gst",
          col=c("blue4","orange3") ,
          logp=F ,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))

write.csv(vcf.div.MHplot, "~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.csv",
          quote=F,
          row.names=F)

#Divsersity within groups

names(vcf.div.MHplot)

vcf.div.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + #why switch from %>% to +? i get we are working on the graph now but still why?
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title="Genome-Wide Expected Heterozygosity (Hs)",fill="Regions",
       x="Gene Diversity within Regions", y="Counts of SNPs")

ggsave("Histograme_GenomDIversity_byRegion.pdf",
       path="~/projects/eco_genomics/population_genomics/figures/")

vcf.div.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + #why switch from %>% to +? i get we are working on the graph now but still why?
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title="Genome-Wide Expected Heterozygosity (Hs)",fill="Regions",
       x="Gene Diversity within Regions", y="Counts of SNPs")

ggsave("Histograme_GenomDIversity_byRegion.pdf",
       path="~/projects/eco_genomics/population_genomics/figures/")



######Part Three- PCA
setwd("~/projects/eco_genomics/population_genomics/")
#path picks up where WD left off 
vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

vcf.thin <- distance_thin(vcf, min.distance = 500)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

#hide the uncompressed VCF file too big for github outside our repo 

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf") # opens -leaves alone

geno <- vcf2geno(input.file="/gpfs1/home/s/s/ssularz/vcf_final.filtered.thinned.vcf",
                 output.file="outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)










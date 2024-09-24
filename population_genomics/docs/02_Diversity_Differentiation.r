# Estimating diversity and genetic differentiation in the filtered centaurea data

library(vcfR)
library(tidyverse)
library(qqman)

#helps solve plotting issues
X11.options(type="cairo")

#read in our VCF file 

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#read in our meta data - info on population of orgin what region the pops came from, what contenetint etc
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) # vcf file has 595 samples 
dim(meta) #meta has 629 inds

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] #%in% "also found in"
dim(meta2)

#calculate diversity stats using genertic_diff fxn in vcfR
vcf.div <- genetic_diff(vcf,
                         pops=as.factor(meta2$region),
                         method = "nei")

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

###09-24-24  Finishing up Diversity 
write.csv(vcf.div.MHplot, "~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.csv",
          quote=F,
          row.names=F)

names(vcf.div.MHplot)

#allows plotting to work on R on VACC
options(bitmapType = "cairo")
#want to visulize a bit more like a density plot 
#use tidy verse to paste each colunm on top of eachother into one long column 
vcf.div.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + #why switch from %>% to +? i get we are working on the graph now but still why?
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title="Genome-Wide Expected Heterozygosity (Hs)",fill="Regions",
       x="Gene Diversity within Regions", y="Counts of SNPs")

#save the last plot you made
ggsave("Histograme_GenomDIversity_byRegion.pdf",
       path="~/projects/eco_genomics/population_genomics/figures/")


#what does playing with the filter mean in words (what are results)
vcf.div.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0 & value<0.25) %>% #take out the zero values changes sample size and tidy it up a bit 
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())









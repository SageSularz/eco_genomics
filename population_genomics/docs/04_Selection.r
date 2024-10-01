library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)

vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf",
                    type = "vcf")
vcfR <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz")
#vcfR <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz")
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcfR@gt[,1]),]

pcadapt.pca <- pcadapt(vcf,
                       K=2,
                       method = "componentwise",
                       min.maf = 0.01, 
                       LD.clumping = list(size=500,thr=0.2))

summary(pcadapt.pca)
plot(pcadapt.pca, options="stat.distribution",
     pop=meta2$region,
     i=1, j=2,
     K=2)

View(head(vcfR@fix))

vcfR.fix <- as.data.frame(vcfR@fix[,1:2])  #chromo position

chr.main <- unique(vcfR.fix$CHROM)[1:8]
chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

Pval <- pcadapt.pca$pvalues

pcadapt.MHplot <- cbind(vcfR.fix, Pval) #####NOT RUNNING
pcadapt.MHplot <- left_join(chrnum, pcadapt.MHplot, join_by(chr.main==CHROM))

pcadapt.MHplot <- pcadapt.MHplot %>%
  mutate(SNP=paste0(chr.main,"_",POS))

str(pcadapt.MHplot)

pcadapt.MHplot$V2= as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS= as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$pPC1= as.numeric(pcadapt.MHplot[,4])
pcadapt.MHplot$pPC2= as.numeric(pcadapt.MHplot[,5])

#M plot will crash if there are any 'NA' vals this can remove them 
pcadapt.MHplot <- pcadapt.MHplot %>% drop_na(pPC1)

#for PC1
manhattan(pcadapt.MHplot,
          chr = "V2",
          bp="POS",
          p="pPC1",
          col = c("blue4", "orange3"),
          logp = T,
          ylab="-log10p-value",
          genomewideline=F,
          main="PCAdapt Genome Scan for Selection (PC1)")


#for PC2
manhattan(pcadapt.MHplot,
          chr = "V2",
          bp="POS",
          p="pPC2",
          col = c("blue4", "orange3"),
          logp = T,
          ylab="-log10p-value",
          genomewideline=F,
          main="PCAdapt Genome Scan for Selection (PC2)")

View(pcadapt.MHplot %>% 
       filter(pPC1<quantile(pcadapt.MHplot$pPC1, 0.001)) %>%
        select(chr.main,POS,pPC1))
#look at chromo 8... peak of outliers can indecate selection 

summary(pcadapt.pca)
#loadings= 
#zscore and chi2 test selection for tests
#allele freq and minor allele freq =af and maf








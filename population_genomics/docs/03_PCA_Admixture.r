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
meta2 <- meta[meta$id %in% colnames (vcf@gt[, -1]),]

dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

#hide the uncompressed VCF file too big for github outside our repo 

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf") # opens -leaves alone

geno <- vcf2geno(input.file="/gpfs1/home/s/s/ssularz/vcf_final.filtered.thinned.vcf",
                 output.file="outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)

#lots of things just happened in not sure about some of them but we will go back into it next class 
#code filled in just theory


CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")
#use this to load pca from working script 

show(CentPCA) #

plot(CentPCA) #axis that pca put through cloud of points is "igan value" the first one 
#the most info content will be pc 1 then pc 2 etc levels off fast so last 600+ not very usefull 

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
       geom_point(alpha=1) +
       labs(title = "Centaurea genetic PCA",x="PC1",y="PC2",color="Region",shape="continent")
#close () after aes then you can use + to layer on after 
##make some edits like this...
#ggplot(as.data.frame(CentPCA$projections),
      # aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
 # geom_point(alpha=0.5) + #spread i think
 # labs(title = "Centaurea genetic PCA",x="PC1",y="PC2",color="Region",shape="continent")+
#  xlim(-10,10) +
 # ylim(-10,10)

ggsave("figures/CentPCA_PC1vPC2.pdf", width = 6, height = 6, units = "in")

#this was from the end of the class before 
plot(CentPCA$projections,
     col=as.factor(meta2$region))
legend("bottomright", legend = as.factor(unique(meta2$region),
                                         fill=as.factor(unique(meta2$region))))
#somethings not right probs just need to proofread 





#####10-1-24

#now we will run admixture analysis and create plots 
#for Admixture we will use LEA R package
#the funct inside LEA is called 'snmf' 

CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno", 
                  K=1 : 10, #tell it what k val you want but hard to find exact val so give it a range
                  entropy = T,
                  repetitions = 3,  #multiple algorithims for val of K and pick the rep with the most cross values
                  project = "new") # if your adding to analysis later you could choose project = "continue"
par(mfrow=c(2,1))                  
plot(CentAdmix, col="blue4", main="SNMF") # result of cross entropy (prediction error) and K val
                #plots cross entropy score we can use for selecting models
plot(CentPCA$eigenvalues[1:10], ylab="Eigenvalues", xlab="Number of PCs", col="blue4")
dev.off() #turn off par otherwise it will keep stacking your plots 

myK=5 

CE = cross.entropy(CentAdmix, K=myK)
best=which.min(CE)

myKQ = Q(CentAdmix, K=myK, run = best) # indv rows the sames but runs in random order

#combine with cbind **make sure rows are in the same order!!
myKQmeta = cbind(myKQ, meta2)

my.colors = c("blue4", "gold", "tomato","lightblue", "olivedrab")

myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region, pop, .by_group = TRUE) # I am not sure what happened here


pdf("figures/Admixture_K5.pdf", width=10, height=5)
barplot(as.matrix(t(myKQmeta[ , 1:myK])),
        border=NA,
        space=0,
        col = my.colors[1:myK],
        xlab="Geographic regions",
        ylab="Ancestory Proportions",
        main=paste0("Ancestry matrix K=",myK))
axis (1,
      at=1:length(myKQmeta$region),
      labels = myKQmeta$region,
      tick=F,
      cex.axis=0.5,
      las=3)
dev.off()


                  
                  
                  
                  
                  

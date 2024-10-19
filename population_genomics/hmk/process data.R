# Ecological Genomics - Homework 1
# Effect of Individual-Level Missingness - Sage Sularz

################# Process Data ##########################

library(vcfR)
library(tidyverse)
library(qqman)
library(SNPfiltR)
library(LEA) 

#helps solve plotting issues
X11.options(type="cairo")

filter_vals <- c(0.75, 0.5, 0.3)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
source("~/projects/eco_genomics/population_genomics/hmk/filter_driver.R")  

for (i in 1:3) {
  filter_val <- filter_vals[i]
  print("start script")
  
  setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
  
  vcf_file <- paste0("~/projects/eco_genomics/population_genomics/hmk/vcf_final.filtered_",filter_val,".vcf.gz")
  
  # Read in the VCF file
  vcf <- read.vcfR(vcf_file)
  
  meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

  meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] #%in% "also found in"

  #calculate diversity stats using genertic_diff fxn in vcfR
  vcf.div <- genetic_diff(vcf,
                          pops=as.factor(meta2$region),
                          method = "nei")


  chr.main <- unique(vcf.div$CHROM)[1:8]

  chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

  vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))

  vcf.div.MHplot <- vcf.div.MHplot %>%
    filter(Gst>0) %>%
    mutate(SNP=paste0(chr.main,"_",POS))
 

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

  write.csv(vcf.div.MHplot, 
            file = paste0("~/projects/eco_genomics/population_genomics/hmk/Genetic_Diff_byRegion_",filter_val,".csv"),
            quote=F,
            row.names=F)

  names(vcf.div.MHplot)


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
  ggsave(file=paste0("Histogram_GenomeDiversity_byRegion_", filter_val, ".pdf"),
         path="~/projects/eco_genomics/population_genomics/hmk/")
  
  #what does playing with the filter mean in words (what are results) #######
  diversity_summary <- vcf.div.MHplot %>%
    as_tibble() %>%
    pivot_longer(c(4:9)) %>%
    group_by(name) %>%
    filter(value!=0 & value<0.25) %>% #take out the zero values changes sample size and tidy it up a bit 
    summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())
  
  write.csv(diversity_summary, 
            file = paste0("~/projects/eco_genomics/population_genomics/hmk/Diversity_Summary_", filter_val, ".csv"),
            row.names = FALSE)
  print("Diversity Loop")

  ######## PCA #########
  print("Starting PCA")

  vcf_thin <- distance_thin(vcf, min.distance = 500)

  meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

  #need to subset meta file
  meta2 <- meta[meta$id %in% colnames (vcf@gt[, -1]),]

  full_path <- "/gpfs1/home/s/s/ssularz"
  
  write.vcf(vcf_thin, 
            file = paste0(full_path, "/projects/eco_genomics/population_genomics/hmk/vcf_final_filtered_thinned_", filter_val, ".vcf.gz"))
  
  # Uncompress the VCF file
  system(paste0("gunzip -c ", full_path, "/projects/eco_genomics/population_genomics/hmk/vcf_final_filtered_thinned_", 
                filter_val, ".vcf.gz > ", full_path, "/vcf_final_filtered_thinned_", filter_val, ".vcf"))
  
  # Check if the uncompressed file exists
  if (!file.exists(paste0(full_path, "/vcf_final_filtered_thinned_", filter_val, ".vcf"))) {
    stop(paste("Error: VCF file not found after uncompression for filter_val =", filter_val))
  } else {
    print(paste("File successfully uncompressed for filter_val =", filter_val))
  }
  
  # Convert the uncompressed VCF to geno format
  geno <- vcf2geno(input.file = paste0(full_path, "/vcf_final_filtered_thinned_", filter_val, ".vcf"),
                   output.file = paste0(full_path, "/projects/eco_genomics/population_genomics/hmk/vcf_final_filtered_thinned_", filter_val, ".geno"))
  
  print(paste("Completed vcf2geno for filter_val =", filter_val))


  setwd("/gpfs1/home/s/s/ssularz")
  
  CentPCA <- LEA::pca(paste0(full_path, "/projects/eco_genomics/population_genomics/hmk/vcf_final_filtered_thinned_", filter_val, ".geno"), scale = TRUE)
  

  plot(CentPCA) #axis that pca put through cloud of points is "igan value" the first one 

  ggplot(as.data.frame(CentPCA$projections),
        aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
    geom_point(alpha=1) +
    labs(title = "Centaurea genetic PCA",x="PC1",y="PC2",color="Region",shape="continent")

  ggsave(file=paste0("CentPCA_PC1vPC2_",filter_val,".pdf"), 
         path="~/projects/eco_genomics/population_genomics/hmk/",
         width = 6, height = 6, units = "in")
  print("PCA Loop")
}

# Ecological Genomics - Homework 1
# Effect of Individual-Level Missingness - Sage Sularz

################# Process Data ##########################

library(vcfR)
library(tidyverse)
library(qqman)
library(SNPfiltR)
library(LEA) 

#helps solve plotting issues
X11.options(type="cairo")

filter_vals <- c(0.75, 0.5, 0.3)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
source("~/projects/eco_genomics/population_genomics/hmk/filter_driver.R")  

for (i in 1:3) {
  filter_val <- filter_vals[i]
  print("start script")
  
  setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
  
  vcf_file <- paste0("~/projects/eco_genomics/population_genomics/hmk/vcf_final.filtered_",filter_val,".vcf.gz")
  
  # Read in the VCF file
  vcf <- read.vcfR(vcf_file)
  
  meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
  
  meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] #%in% "also found in"
  
  #calculate diversity stats using genertic_diff fxn in vcfR
  vcf.div <- genetic_diff(vcf,
                          pops=as.factor(meta2$region),
                          method = "nei")
  
  
  chr.main <- unique(vcf.div$CHROM)[1:8]
  
  chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))
  
  vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))
  
  vcf.div.MHplot <- vcf.div.MHplot %>%
    filter(Gst>0) %>%
    mutate(SNP=paste0(chr.main,"_",POS))
  
  
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
  
  write.csv(vcf.div.MHplot, 
            file = paste0("~/projects/eco_genomics/population_genomics/hmk/Genetic_Diff_byRegion_",filter_val,".csv"),
            quote=F,
            row.names=F)
  
  names(vcf.div.MHplot)
  
  
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
  ggsave(file=paste0("Histogram_GenomeDiversity_byRegion_", filter_val, ".pdf"),
         path="~/projects/eco_genomics/population_genomics/hmk/")
  
  #what does playing with the filter mean in words (what are results) #######
  diversity_summary <- vcf.div.MHplot %>%
    as_tibble() %>%
    pivot_longer(c(4:9)) %>%
    group_by(name) %>%
    filter(value!=0 & value<0.25) %>% #take out the zero values changes sample size and tidy it up a bit 
    summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())
  
  write.csv(diversity_summary, 
            file = paste0("~/projects/eco_genomics/population_genomics/hmk/Diversity_Summary_", filter_val, ".csv"),
            row.names = FALSE)
  print("Diversity Loop")
  
  ######## PCA #########
  print("Starting PCA")
  
  vcf_thin <- distance_thin(vcf, min.distance = 500)
  
  meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
  
  #need to subset meta file
  meta2 <- meta[meta$id %in% colnames (vcf@gt[, -1]),]
  
  full_path <- "/gpfs1/home/s/s/ssularz"
  
  write.vcf(vcf_thin, 
            file = paste0(full_path, "/projects/eco_genomics/population_genomics/hmk/vcf_final_filtered_thinned_", filter_val, ".vcf.gz"))
  
  # Uncompress the VCF file
  system(paste0("gunzip -c ", full_path, "/projects/eco_genomics/population_genomics/hmk/vcf_final_filtered_thinned_", 
                filter_val, ".vcf.gz > ", full_path, "/vcf_final_filtered_thinned_", filter_val, ".vcf"))
  
  # Check if the uncompressed file exists
  if (!file.exists(paste0(full_path, "/vcf_final_filtered_thinned_", filter_val, ".vcf"))) {
    stop(paste("Error: VCF file not found after uncompression for filter_val =", filter_val))
  } else {
    print(paste("File successfully uncompressed for filter_val =", filter_val))
  }
  
  # Convert the uncompressed VCF to geno format
  geno <- vcf2geno(input.file = paste0(full_path, "/vcf_final_filtered_thinned_", filter_val, ".vcf"),
                   output.file = paste0(full_path, "/projects/eco_genomics/population_genomics/hmk/vcf_final_filtered_thinned_", filter_val, ".geno"))
  
  print(paste("Completed vcf2geno for filter_val =", filter_val))
  
  
  setwd("/gpfs1/home/s/s/ssularz")
  
  CentPCA <- LEA::pca(paste0(full_path, "/projects/eco_genomics/population_genomics/hmk/vcf_final_filtered_thinned_", filter_val, ".geno"), scale = TRUE)
  
  
  plot(CentPCA) #axis that pca put through cloud of points is "igan value" the first one 
  
  ggplot(as.data.frame(CentPCA$projections),
         aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
    geom_point(alpha=1) +
    labs(title = "Centaurea genetic PCA",x="PC1",y="PC2",color="Region",shape="continent")
  
  ggsave(file=paste0("CentPCA_PC1vPC2_",filter_val,".pdf"), 
         path="~/projects/eco_genomics/population_genomics/hmk/",
         width = 6, height = 6, units = "in")
  print("PCA Loop")
}







  



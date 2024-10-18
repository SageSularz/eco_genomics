
# load library
library(eulerr)

# set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)
    # [1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28" "group_D22A28_vs_D18A28" 
    # [5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

# 1. compare baseline gene expression between developmental treatment groups
res_D18_BASE_D22_BASE <- results(dds, contrast = c("group", "D18BASE", "D22BASE"), alpha = 0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE)

#make a list of which genes in our comparisons of interest are differentially exposed (list of DEGs)
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < 0.05,])

plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))



############ some wack shit just happened but it was right around 48 min


# 2. compare baseline gene expression between developmental treatment groups at A28
res_D18_BASE_D22_A28 <- results(dds, contrast = c("group", "D18BASE", "D22A28"), alpha = 0.05)
res_D18_BASE_D22_A28 <- res_D18_BASE_D22_A28[!is.na(res_D18_BASE_D22_A28$padj),]
res_D18_BASE_D22_A28 <- res_D18_BASE_D22_A28[order(res_D18_BASE_D22_A28$padj),]
head(res_D18_BASE_D22_A28)
summary(res_D18_BASE_D22_A28)

#make a list of which genes in our comparisons of interest are differentially exposed (list of DEGs)
degs_D18_BASE_D22_A28 <- row.names(res_D18_BASE_D22_A28[res_D18_BASE_D22_A28$padj < 0.05,])

plotMA(res_D18_BASE_D22_A28, ylim=c(-4,4))















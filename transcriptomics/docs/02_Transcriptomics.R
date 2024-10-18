#load in Day 1 transcriptomic script and run DESseq object 
library(pheatmap)
options(bitmapType = "cairo")

resultsNames(dds) 
# [1] "Intercept"      "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

# pull out the results for Developmental temp 22 vs 18 
res_D22vsD18 <- results(dds, name="DevTemp_D22_vs_D18", alpha = 0.05)

#order by significance (specifically difference in xpression between samples)
#what is the gene expression in 22 compared to 18
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18)

# look at counds of a specific top gene that we are interested in to validate that the model is working 
d <- plotCounts(dds,gene="TRINITY_DN140854_c0_g5_i2", int=(c("DevTemp", "FinalTemp")), returnData = TRUE)
d

p <- ggplot(d,aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) +
  theme_minimal() + theme(text=element_text(size=20), panel.grid.major=element_line(color = "grey")) 

p <- p+ geom_point(position=position_jitter(w=0,h=0),size=3)
p

# MAplot amount of logfold change 
plotMA(res_D22vsD18, ylim=c(-4,4))
#far right is highly xpressed but not as many 
#good example of overdispersion 


# Volcano plot 

# Convert our deseq results into a dataframe to plot 
res_df <- as.data.frame(res_D22vsD18)

# add a column to dataframe to denote if a gene is significantly differently xpressed 
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >1, "Significant", "Not Significant")

#plot 
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color= Significant))+
  geom_point(alpha=0.8) +
  scale_color_manual(values=c("slateblue", "tomato"))+
  labs(x="Log2 Fold Change", y="-log10 Adj P-Value", title = "Volcano Plot")+
  theme_minimal()+
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color="orange") +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="orange")


#heat map 
vsd <- vst(dds, blind = FALSE)

#all the genes would be too much to analyze
topgenes <- head(rownames(res_D22vsD18), 20)
mat <- assay(vsd)[topgenes,]
df <- as.data.frame(colData(dds)[c("DevTemp", "FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cols=T, cluster_rows=T)
#columns are samples rows are genes
#tree: phenogram trying to group samples by gene expression patterns




















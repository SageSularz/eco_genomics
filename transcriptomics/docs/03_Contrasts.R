###
### 

# check out of the library 
library(eulerr)
library(dplyr)
library(tidyr)
library(gridExtra)

# set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)
    # [1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28" "group_D22A28_vs_D18A28" 
    # [5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

# 1. compare baseline gene expression between developmental treatment groups
res_D18_BASE_D22_BASE <- results(dds, contrast=c("group", "D18BASE", "D22BASE"), alpha = 0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),] 
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < 0.05,])

plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))

# 2. compare A28 gene expression between developmental treatment groups
res_D18_A28_D22_A28 <- results(dds, contrast=c("group", "D18A28", "D22A28"), alpha = 0.05)
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),] 
head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < 0.05,])

plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))

# 3. compare A33 gene expression between developmental treatment groups
res_D18_A33_D22_A33 <- results(dds, contrast=c("group", "D18A33", "D22A33"), alpha = 0.05)
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),]
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),] 
head(res_D18_A33_D22_A33)
summary(res_D18_A33_D22_A33)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj < 0.05,])

plotMA(res_D18_A33_D22_A33, ylim=c(-4,4))

#tota;
length(degs_D18_BASE_D22_BASE) #1935
length(degs_D18_A28_D22_A28) #296
length(degs_D18_A33_D22_A33) #78

#look at the overlaps in which genes are differentially expressed in multiple contrasts
length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A28_D22_A28)) #107
length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A33_D22_A33)) #44
length(intersect(degs_D18_A33_D22_A33, degs_D18_A28_D22_A28)) #29

nested_intersection <- intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)
length(intersect(degs_D18_A33_D22_A33, nested_intersection)) #23


################# 10-22-24 ######################
##
# Making Euler Plots 
##

# calculate the number of unique genes in each portion of the euler plot 
1935-107-44+23 #1807 genes differentially expressed uniquely at baseline btwn 18  vs 22
296-107-29+23 #183 uniquly expressed when exposed to 28 
78-44-29+23 #28 uniquely expressed when exposed to 33

107-23 #84 unique to baseline and and 28 
44-23 #21 unique to baseline and A33
29-23 #6 unique to A28 and A33

# we now have all the pieces we need to plot!
myEuler <- euler(c("BASE"=1807, "A28"=183, "A33"=28, "BASE&A28"=84, 
                   "BASE&A33"=21, "A28&A33"=6, "BASE&A28&A33"=23))
plot(myEuler, lty=1:3, quantities=TRUE)


################ 
#
#Make a scatterplot of the responses to A28/A33 when coepads develop at 18 vs 22
#
##############


# contrast D18_A28vsBASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group", "D18BASE", "D18A28"), alpha=0.05))

# contrast D22_A28vsBASE

res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group", "D22BASE", "D22A28"), alpha=0.05))

# merge dataframes
res_df28 <- merge(res_D18_BASEvsA28, res_D22_BASEvsA28, by="row.names", suffixes=c(".18", ".22"))
View(res_df28)
rownames(res_df28) <- res_df28$Row.names
res_df28 <- res_df28[,-1]


# make col collumn that will col based on values in data fram 

# Define color mapping logic with the mutate fuction 
res_df28 <- res_df28 %>%
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2",
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta3",
    padj.22 < 0.05 & stat.22 < 0 ~ "blue",
    padj.22 < 0.05 & stat.22 > 0 ~ "tomato",
  ))

# count the number of color points per fill color 
color_counts <- res_df28 %>%
  group_by(fill) %>%
  summarise(count=n())
label_positions <- data.frame(
  fill=c("blue", "magenta3", "tomato", "turquoise2"),
  x_pos = c(1, 5, 0, -7.5),
  y_pos = c(-5, 0, 9, 3)
)

label_data <- merge(color_counts, label_positions, by = "fill")

# PLOT
##Add ins
plot28 <- ggplot(res_df28, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y=y_pos, label=count, color=fill),   ##
            size=5) +    
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  ##
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x="Log2FoldChange 28 vs. BASE at 18",
       y="log2FoldChange 28 vs. BASE at 22",
       title = "How does response to 28 C vary by DevTemp?")+
  theme_minimal()

plot28
##################################
#
# 10-24-24
#
##################################

### Repeat for A33 ###

# contrast D18_A33vsBASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group", "D18BASE", "D18A33"), alpha=0.05))

# contrast D22_A33vsBASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group", "D22BASE", "D22A33"), alpha=0.05))

# merge dataframes
res_df33 <- merge(res_D18_BASEvsA33, res_D22_BASEvsA33, by="row.names", suffixes=c(".18", ".22"))
View(res_df33)
rownames(res_df33) <- res_df33$Row.names
res_df33 <- res_df33[,-1]


# Define color mapping logic with the mutate fuction 
res_df33 <- res_df33 %>%
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2",
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta3",
    padj.22 < 0.05 & stat.22 < 0 ~ "blue",
    padj.22 < 0.05 & stat.22 > 0 ~ "tomato",
  ))

# count the number of color points per fill color 
color_counts <- res_df33 %>%
  group_by(fill) %>%
  summarise(count=n())

label_positions <- data.frame(
  fill=c("blue", "magenta3", "tomato", "turquoise2"),
  x_pos = c(1, 5, 0, -7.5),
  y_pos = c(-5, 0, 9, 3)
)

label_data <- merge(color_counts, label_positions, by = "fill")

# PLOT
plot33 <- ggplot(res_df33, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y=y_pos, label=count, color=fill),   ##
            size=5) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  ##
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x="Log2FoldChange 33 vs. BASE at 18",
       y="log2FoldChange 33 vs. BASE at 22",
       title = "How does response to 33 C vary by DevTemp?")+
  theme_minimal()

plot33

#put the two plots together in a two panel plot 

combined_plot <- grid.arrange(plot28,plot33, ncol = 2)
ggsave("~/projects/eco_genomics/transcriptomics/figures/combined_scatter_plot.png", combined_plot, width = 12, height = 6)











































## Transcriptomics

### Easy Refs:

1.  Issues with plots on Rstudio on VACC -\> options(bitmapType = "cairo")

### 

### 10-09-24 Intro

Physiological meachanisms of development plasticity "How do animals acclimate?"

potential experimental questions

-   does the temp that they expereince growing up affect ULT?

-   How does gene expression response differ between 28 and 33 degrees and does this differ with baseline?

-   What genes are differentially expressed at DT 22 compared to DT 18?

factors

-   development temp

    -   levels = 18 & 22

-   Final temp

    -   levels = baseline (BASE) & 28 (A28) & 33 (A33)

### 10-10-24 '

explore counts matrix

-   average number of reads per sample 18,454,529

-   the average number of counts per gene: mean= 3,244,739 med= 64

line 60 DESeq

-   estimating size factors

-   estimating dispersions

-   gene-wise dispersion estimates

-   mean-dispersion relationship

-   final dispersion estimates fitting model and testing

### 10-15-24

log2foldchange- measure of difference in gene expression

-   does everything in doublings

-   make sure you check the direction!!

    -   confirm by plotting

<!-- -->

-   everything is a comparason between two catagories

looking at the result names and compared developmental temps

### 10-17-24 Contrasting

notes:

Qs:

### 10-21-24

Notes:

Plots review:

1.  Bar plot of number of reads per sample

    -   to visualize the sucess and variation of sequencing effort

2.  PCA

    -   visualize the variation among groups

3.  MA plot

    -   x = average counts y = LogFoldChange

    -   difference in expression between two groups

    -   and how that difference related to how highly expressed that gene is

    -   positive vals upreg in 22 aka down in 18

    -   negative val upreg in 18 aka down in 22

4.  Point plot of a differentially expressed gene

5.  Volcano plot

6.  Euler plot

    -   shared vs unique genes among different contrasts

7.  Heat map

    -   see by color and saturation diff among genes or matrix of data (ex. correl. WGCNA)

    -   we can visualize almost any matrix of information

8.  Scatterplots

    -   can be colored by contrast and significance

### 10-24-24

Plan for the day:

1.  Review plots
2.  Finish scatter plot
3.  WGCNA (fill in steps form notes) (GOT HERE TODAY)
4.  Two ways to approach gene ontology (GO) enrichment analysis
5.  Make files to prepare for GO analysis

Notes:

-   LFC \>\> pos val = pos test stat

-   wrapped up "03_Contrasts" started "04_WGCNA"

-   cluster dendrogram plot

-   selecting a power...

    -   plot will help you select the best biological threshhold

-   next week we can idetify the moduals then test for corilation

Q's:

-   grid.arrange funct.\>\> play with this more, how to rotate? arrange in different ways?

### 10-29-24

*Finishing up WGCNA*

Notes:

-   pick power that has more bio topograph without getting too low on conectivity

Qs:

# Coding and data notes for the population genomics module

#Author: Sage Sularz

### Easy Refs:

1.  Issues with plots on Rstudio on VACC -\> options(bitmapType = "cairo")
2.  Goal is always to find the least complicated model you can

### 09-10-2024 - Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS data from 3 regions (eu, NE, PNW) starting today with variant call format files (VCF's)

#this was a virtual day - light notes #dear god dont forget the head command again #sffr if you take one thing away from today let it be that \^\^

### 09-12-24 Viewing VCF Files and talking about filtering

### 09-17-24 Finishing up with filtering and VCFs then getting into diversity stats

0/0=ref hom, 0/1=het, 1/1= alt hom, ./.= NA

#grouping-- df[rows,columns] df[ ,c(1,4)]

Filtering Stratigies== /\_\\ find middle zone with pros and trade offs of each

#depth: low dp (limited accuracy of GF call)/// use hard_filter funct High dp (assembly error, paralogy-reads from duplicated genes mapping to same position) use max_depth func #missingness - indevidual level and SNP level

#low frequency alleles

#Questions? #\>\> how does file types effect head or similar commands?

-   start to visualize values associatied with each "triangle point"

### 09-19-24 Diversity

Notes:

-   matching data sets using %in%

-   Manhat plot

    -   what regions are doing something different? Higher FST in the outer edges of chromos \> Telomeres? I was distracted but i think this seems right

    -   this is cool follow up on this later

-   colnames([vcf\@gt](mailto:vcf@gt){.email}[,-1]) \> minus first column alter amount of info

1.  Pulled vcf from 01 script

Questions:

-   what is the difference between View and view (caps why?)

Objectives (from BS):

Step 1: Read in our filtered VCF file (in your repo)

Step 2: Read in our metadata file (in the class datashare on /gpfs1)

Step 3: Make sure they're the same \# of samples (remember...we filtered some out! This is only something you have to dbl-check)

Step 4: Run genetic_diff() on the vcf file with the sample metadata provided for grouping -- which groups should we run? We have options...

Step 5: Process the resulting output to get it ready for plotting using tidyverse commands...

Step 6.a: Visualize gene diversity (Hs) across the genome and among our different groups.\
*Note that Hs is equivalent to Hexp (expected heterozygosity) when loci are biallelic (max 2 alleles)*

Step 6.b: Visualize genetic differentiation (Gst) across the genome and among our different groups.\
*Note that Gst is equivalent to Fst when loci are biallelic*

Step 7: Tweak plots -- save interesting ones to our outputs/ directories in our repo!

### 09-24-24 Diversity and starting PCA

Notes:

-   %\>% use as a pipline in tidy
-   missed some context on ggplot follow up and think about this a bit
-   Manhat plot
    -   outliers may be genes of interest maybe bc of selection
-   

Questions:

-   line 63 \> still kinda confused on this flow and the %\>% funct.
-   pivot_longer make big column, how to make big row? does it even matter?
-   how to determine min distance

### 09-26-24 PCA/Admixture

Notes

-   PCA plot

    -   doesnt make assumption about where samples taken

    -   nice to look at genetic pop struc in a different way from fst without having to define groups

    -   

-   how to discribe pc vals

    -   each snp have a corrilation with each axis

    -   if you sum all the eigan values you get vectors\*

-   Admixture analysis

    -   also good way to look at pop struct with groups

    -   has a genetic model behind it (unlike pca which is just maths) -\> HWE

    -   HWE is good to prediction het content from 1=2pq

Groups are

1.  choose value of k (k=# groups) (1-10)
2.  assign indvs to one of k groups (up to you )
3.  calc allele freq in each group (p, q)
4.  calc 2pq from that & compare to obv freq of het
5.  swing back to step two and try to regroup to try and increase similarity to obv vals

-   Q= fractional ansestory (whaa)

    -   make a matrix outlineing how much each indv fits into each group, ex 0.5 equal between two groups

-   cross-validation in model training

    -   connected to Q and k

-   Screeplot: shows the mag of the eigenvalues in decending order, starting with the first pc value. we only want the first few pc values they contain most the info

Questions

-   how can i view the PCA plot in multiple dementions i need to see it all together to make sense of it i think

-   still unclean on how to pick the right K

how not to over interpret structure plots \<- try this paper to fill in some blanks

### 10-01-24  Admixture/selection

Agenda

1.  admixture analysis
    -   calculation

    -   plotting
2.  micro lecture on selection outliers
3.  calculate selection outliers
4.  homework

NOTES:

-   What to plug in for K range

    -   every extra K you run will take longer and you will see the fit of the data will improve as K increases until a certain point then level out (See admix plot)

    -   get sense of pca and eigan vals to use as rough guide

    -   "look for elbow on plot" and select that as range

-   Admix plot

    -   scale is less important compared to proportion between points

    -   compared to last weeks plot that looked like 'U' this plot would eventually move back up with error and create 'U'

-   look at the admix and PCA plots...

    -   is there overlap in conclusions I can pull?

    -   when when I use one vs the other

    -   do you understand the primary goals in each

    -   what features am I specifically look at in each

-   

QUESTIONS:

-   Monte carlo analysis was mention I would like to understand the math behind it

-   used cbind to merge on the y axis how can i merge on the x axis

SELECTION NOTES

-   NCBI NLM genome data viewer from eod find link on BS

    -   view annotated genome data

    -   what gene is that outlier on? Whats happening there? Is this signiicant?

-   Kegg or GO are good enrichment pathways

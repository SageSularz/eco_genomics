# Coding and data notes for the population genomics module

#Author: Sage Sularz

### Easy Refs:

1.  Issues with plots on Rstudio on VACC -\> options(bitmapType = "cairo")

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

### 09-24-24 Diversity 

Notes:

-   %\>% use as a pipline in tidy
-   missed some context on ggplot follow up and think about this a bit
-   Manhat plot
    -   outliers may be genes of interest maybe bc of selection
-   

Questions:

-   line 63 \> still kinda confused on this flow and the %\>% funct.
-   pivot_longer make big column, how to make big row? does it even matter?
-   

# Coding and data notes for the population genomics module

#Author: Sage Sularz

##09-10-2024 - Intro to Centaurea GBS data and working with VCF files #We'll be analyzing the GBS data from 3 regions (eu, NE, PNW) starting today with variant call format files (VCF's)

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

-   colnames([vcf\@gt](mailto:vcf@gt){.email}[,-1]) \> minus first column alter amount of info

1.  Pulled vcf from 01 script

Questions:

-   what is the difference between View and view (caps why?)

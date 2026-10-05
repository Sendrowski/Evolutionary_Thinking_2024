---
title: "The Origins and Molecular Evolution of SARS-CoV-2"
output:
  html_notebook:
    fig_caption: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
  pdf_document:
    toc: yes
    toc_depth: '3'
  html_document:
    keep_md: TRUE
---

## Author: Michi Tobler, Division of Biology Professor, Kansas State University

### Curated by: Calin Pantea

## **Quick note: This exercise focuses on multiple themes; I left most of the ones you have already encountered in, for extra context. It looks pretty well-made and comprehensible; but you are also able to jump straight to the selection part (4.2) if you wish, without having to run anything above.**

------------------------------------------------------------------------

# 0. Installing Required Packages

**Note: Once you have successfully installed all components, I recommend that you delete this section, including the code snippet below. Alternatively, you can silence the code using hashtags.**


``` r
install.packages("ape")
install.packages("dplyr")
install.packages("phytools")
install.packages("tidyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

```
## Warning: package(s) not installed when version(s) same as or greater than current; use `force = TRUE` to
##   re-install: 'ggtree'
```

```
## Update all/some/none? [a/s/n]: Update all/some/none? [a/s/n]: Update all/some/none? [a/s/n]: Update all/some/none? [a/s/n]: Update all/some/none? [a/s/n]: Update all/some/none? [a/s/n]:
```

```
## Warning in readLines(n = 1): type 29 is unimplemented in 'type2char'
```

```
## Error in readLines(n = 1): INTEGER() can only be applied to a 'integer', not a 'unknown type #29'
```

## **Note: If any of the packages fail to install or to load due to outdated dependencies, try uninstalling the dependency and then reinstalling the required version. Direct update attempts seem to sometimes fail or make RStudio hang/loop.**

------------------------------------------------------------------------

# 1. Initiate the Project


``` r
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)
library(phytools)
library(tidyr)

# set your working directory to the path leading to the data files:
setwd("")
```

```
## Error in setwd(""): cannot change working directory
```

``` r
getwd()
```

```
## [1] "/Users/au732936/PycharmProjects/Evolutionary_Thinking_2024/week41/Friday"
```

## **Note: you do not have to go through any of this with a technical focus; the main focus of the exercise is to get a feel for how to interpret - or draw conclusions based on - trees, the molecular clock, and signatures of selection. I have included some basic guideline answers at the end for the interpretative questions, but left everything else in. This is not an R course, so you don't have to know how to write any of the code. But do try to follow the way everything is set up and understand how we obtain our results :)

------------------------------------------------------------------------

# 2. Where Did COVID-19 Come From?

Coronaviruses are common pathogens in humans and mostly have benign effects (e.g., the common cold). But some past coronavirus strains have caused substantial illness in humans, including MERS (Middle East Respiratory Syndrome), which caused a limited outbreak in 2012, and SARS (Severe Acute Respiratory Syndrome), which caused a limited outbreak in 2002/2003. The ubiquity of coronaviruses raises questions about the origins of SARS-CoV-2, which is causing the current COVID-19 pandemic.

------------------------------------------------------------------------

## 2.1. Hypothesis

Based on the information above and what you know about the [evolution of flu viruses](https://michitobler.github.io/primer-of-evolution/evolution-of-dna-sequences.html#the-evolution-of-the-seasonal-flu), phrase a hypothesis about the potential evolutionary origins of SARS-CoV-2.

*Your answer goes here*

------------------------------------------------------------------------

## 2.2. Identifying the Evolutionary Origins of SARS-CoV-2

Phylogenetic analyses allow us to examine how organisms are related, and in the case of SARS-CoV-2, where its evolutionary origins might lie. In the modern age of DNA sequencing, virus genomes from the current pandemic have become rapidly available. In the section below, you will visualize a phylogenetic tree that is based on whole genome sequencing of different SARS-CoV-2 isolates, isolates of SARS-CoV (which caused the outbreak in 2002/2003), as well as other SARS-like coronaviruses that were isolated from a variety of host species. We will walk you through the visualization and annotation of the tree step by step in the following section.

------------------------------------------------------------------------

### 2.2.1. Importing and Plotting a Phylogenetic Tree

Phylogenetic trees come in a wide variety of formats that are not easily readable with standard programs. In this exercise, we will use trees that are encoded in the Newick format (nwk file ending), and R can read and interpret such trees with the `read.tree()` function from the `ape` package. To load the phylogeny of the different coronaviruses, you can simply execute the code below, and a object `cov.tree` should show up in your Global Environment.


``` r
#This line simply reads in the phylogenetic tree (akin to the read.csv function for reading in tables)
cov.tree <- read.tree("coronaviruses_tree.nwk")
```

```
## Warning in file(file, "r"): cannot open file 'coronaviruses_tree.nwk': No such file or directory
```

```
## Error in file(file, "r"): cannot open the connection
```

Once your tree is imported into R, you can visualize it using the `ggtree()` function from the `ggtree` package. `ggtree` is an extension of `ggplot2`, which you have already worked with. You only have to specify the data frame containing the phylogenetic tree (which we named `cov.tree` above). Like with `ggplot()`, you can add graphical elements using the plus sign. For example, the names of the different viral isolates are added to the tips of the phylogeny using `geom_tiplab()`.


``` r
#You can plot phylogenetic trees with the ggtree function (it works just like ggplot). Note that "size=2" associated with geom_tiplab just modifies the font size
p1 <- ggtree(cov.tree) + geom_tiplab(size=2)
```

```
## Error: object 'cov.tree' not found
```

``` r
p1
```

```
## Error: object 'p1' not found
```

At this point, we can see the tree structure. However, the tip labels in this case have limited information, because they are mostly just weird acronyms that researchers used to designate their samples. To address our question above, we need to add some additional information to our tree.

------------------------------------------------------------------------

### 2.2.2. Importing Metadata

Files containing phylogenetic trees only include the names of the taxa included and the length and topology of the branches that connect them. Additional information about each sample (the so-called metadata) is stored in a separate file. You will see that the first column in the metadata file corresponds to the tip labels in the phylogeny. The metadata for this tree includes the country of origin, the host species, virus type, the authors that collected the original data, and the accession (a unique code that allows you to find the specific sample in designated databases). The metadata is stored in a simple table that you can import with the `read.csv()` function.


``` r
#This line reads the metadata of samples included in the phylogenetic tree
cov.tree.meta<- read.csv("coronaviruses_tree_meta.csv")
```

```
## Warning in file(file, "rt"): cannot open file 'coronaviruses_tree_meta.csv': No such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

``` r
cov.tree.meta
```

```
## Error: object 'cov.tree.meta' not found
```

------------------------------------------------------------------------

### 2.2.3. Annotating the Phylogenetic Tree

We can use the metadata to annotate the phylogenetic tree in a way that allows us to address our central question. In the example below, I am using the variable `virus.type` to color the tip of the phylogeny (adding `geom_tippoint()`), and rather than labeling the tips with the name of the strain, I am labeling them with the host species the virus was isolated from (using `geom_tiplab()`).


``` r
#This line of code simply plots the same tree as above, but without the sample label
p2 <- ggtree(cov.tree)
```

```
## Error: object 'cov.tree' not found
```

``` r
#This codes adds the metadata to the tree
p2 %<+% cov.tree.meta + 
  geom_tiplab(aes(label=host), size=2) +
  geom_tippoint(aes(color=virus.type), na.rm = TRUE)
```

```
## Error: object 'cov.tree.meta' not found
```

------------------------------------------------------------------------

## 2.3. Interpretation

What does the phylogeny tell you about the evolutionary origins of SARS-CoV-2, which causes the current pandemic? Discuss the findings in the context of your hypothesis.

*Your answer goes here.*

------------------------------------------------------------------------

# 3. Molecular Clock: How Fast Do Mutations Accumulate?

Unlike most previous disease outbreaks, our ability to sequence a large number of virus genomes has allowed us to literally observe the evolution of the pathogen live. Visit the [the Nextstrain graphical interface](https://nextstrain.org/ncov/global) to see how the pathogen spread across the globe and how mutations accumulate as different viral clades diversify.

The [neutral theory](https://michitobler.github.io/primer-of-evolution/evolution-of-dna-sequences.html#the-neutral-theory) predicts that mutations steadily accumulate during the course of evolution. In this section, you will analyze data from 2,625 whole-genome sequences of SARS-CoV-2 that researchers have collected all over the world. The data set includes the original isolate from December 26 2019 and representative isolates up until July 6 2022.

------------------------------------------------------------------------

## 3.1. COVID-19 Genomes

In the first step, let's visualize the phylogenetic relationships among all of these strains in the context of where they were collected. Hold your breath... this will be a mess:) No surprise when you plot over 2,500 data points at once!

Note this uses the same approach as above. You import the tree and the metadata. Then you plot the tree and add the metadata information to the tree. The code below is complete and should run automatically. Note that the scale bar represents the number of mutations between isolates.


``` r
#Load the phylogenetic tree
cov19.tree <- read.tree("nextstrain_ncov_open_global_6m_tree.nwk")
```

```
## Error in file(file, "r"): cannot open the connection
```

``` r
#Load the corresponding metadata
cov19.tree.meta <- read.csv("nextstrain_ncov_open_global_6m_metadata.csv")
```

```
## Error in file(file, "rt"): cannot open the connection
```

``` r
#Visualize the tree and color the tips by region
#Note that I am using a fan layout here rather than the traditional layout, because it spreads out the tips more
#Depending on your monitor size, you may want to adjust the size of the tip points using the size parameter (currently set to 1.5)
p3 <- ggtree(cov19.tree, layout="fan", open.angle=10)  + geom_treescale(offset=20)
```

```
## Error: object 'cov19.tree' not found
```

``` r
p3 %<+% cov19.tree.meta +
  geom_tippoint(aes(color=region), na.rm = TRUE, size=1.5)
```

```
## Error: object 'cov19.tree.meta' not found
```

------------------------------------------------------------------------

## 3.2. Quantifying Mutations

To examine how mutations accumulate during virus evolution, we need to tally up the number of nucleotide substitutions for each isolate relative to the original virus. Using the phylogenetic data we have available, this is really simple, because we can sum up the length of branches that separates each isolate from the root of the tree (i.e., the original virus sequence).

The code you need to do this is really simple; it requires a single function (`distRoot()`). However, I do not recommend that you run this on your computer! For 2,625 genomes, this calculation takes quite a while, and depending on your hardware, it might just crash your computer. So, I already ran this analysis for you, and you can simply import the results in the next section.


``` r
##YOU DO NOT NEED TO RUN THIS! THIS TAKES A LONG TIME TO RUN. THIS CODE IS HERE FOR ILLUSTRATION ONLY, IN CASE YOU ARE INTERESTED HOW IT WORKS

#library(adephylo)
#dist <- distRoot(cov19.tree, method = "patristic")
#ncov_mutations <- merge (cov19.tree.meta, dist, by.x = "strain", by.y = 0)
#ncov_mutations$date2 <- as.Date(ncov_mutations$date, format = "%Y-%m-%d")
#ref <- as.Date("2019-12-26", format = "%Y-%m-%d")
#ncov_mutations$time <- ncov_mutations$date2-ref
```

------------------------------------------------------------------------

## 3.3 Plot the Data

If we want to know how fast SARS-CoV-2 accumulates mutations over time, we can estimate the genome-wide molecular clock (i.e., the rate at which mutations occur over time). This can be done by plotting the number of mutations separating any isolate from the root as a function of time and then calculating the slope of a regression line between time (in days) and the number of mutations. We will do this using virus isolates for the first two years of the pandemic.

First, let's import and plot the data using a scatter plot. The data file "ncov_mutations.csv" includes the number of mutations separating any isolate from the root (calculated from the phylogeny as described above) as well as the number of days that passed between the collection of the original isolate. Remember, to draw a regression line, you can add `geom_smooth(method=lm, formula= y~0+x, se=FALSE)` to the code that makes a scatter plot.


``` r
#Import data
cov19.mutations <- read.csv("ncov_mutations.csv")
```

```
## Warning in file(file, "rt"): cannot open file 'ncov_mutations.csv': No such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

``` r
#Select first two years of data
cov19.mutations2 <- cov19.mutations[which(cov19.mutations$DaysSinceRoot<730),]
```

```
## Error: object 'cov19.mutations' not found
```

``` r
#Plot data
ggplot(cov19.mutations2, aes(x=DaysSinceRoot, y=Mutations)) + 
  geom_point() +
  geom_smooth(method="lm", se=FALSE, formula=y~x-1) + 
  theme_classic() + 
  xlab("Days") + 
  ylab("Mutations")
```

```
## Error: object 'cov19.mutations2' not found
```

------------------------------------------------------------------------

## 3.4. Estimate the Slope

The regression line essentially represents the average at which rate mutations accumulate. Numerically, the number of mutations per day (i.e., the slope of the regression line in the graph above) can be calculated with the code below.


``` r
#You can run a regression with the lm function, designating the relationship between the variables as "y ~ 0 + x". The 0 indicates that the regression line must go through the origin.
regression <- lm(cov19.mutations2$Mutations ~ 0 + cov19.mutations2$DaysSinceRoot)
```

```
## Error in eval(predvars, data, env): object 'cov19.mutations2' not found
```

``` r
#The slope of the regression (in number of mutations per day) can be found by calling regression$coefficients[1]. Multiplying that number by 365 gives us a substitution rate per year.
mutations.per.year <- 365*regression$coefficients[1]
```

```
## Error: object 'regression' not found
```

``` r
mutations.per.year
```

```
## Error: object 'mutations.per.year' not found
```

------------------------------------------------------------------------

## 3.5. Interpretation

------------------------------------------------------------------------

### 3.5.1. General Patterns

As you evaluate the analyses in this section, what do you notice? What happened to all the viruses that were highly similar to the original samples and common in the first 100 days of the pandemic?

*Your answer goes here*

------------------------------------------------------------------------

### 3.5.2. Application of the Molecular Clock

The calculation of the rate of mutation accumulation was really straightforward for the first two years of the pandemic. However, let's plot the data for the full range of the dataset. I am color coding the data for the first two years separately:


``` r
#Here, I am just making a new variable for first two years or not
cov19.mutations$group <- ifelse(cov19.mutations$DaysSinceRoot<730, "Early", "Late")
```

```
## Error: object 'cov19.mutations' not found
```

``` r
#Plot data
ggplot(cov19.mutations, aes(x=DaysSinceRoot, y=Mutations, color=group)) + 
  geom_point(aes()) +
  geom_smooth(data=subset(cov19.mutations,group=="Early"), method="lm", se=FALSE, formula=y~x-1, color="black") + 
  theme_classic() + 
  xlab("?") + 
  ylab("?")
```

```
## Error: object 'cov19.mutations' not found
```

With your knowledge of the mutation rate, how do you interpret the results? Why are there apparent jumps in the number of mutations of isolates past day 730?

*Your answer goes here*

------------------------------------------------------------------------

# 4. Distribution of Mutations and Positive Selection

Now you saw how SARS-CoV-2 has continuously accumulated mutations as it was spreading across the globe. The next question of course is where in the viral genomes do these mutations occur?

------------------------------------------------------------------------

## 4.1. Distribution of Mutations

Using the 2,625 genomes that were sequenced, we can literally tally up where in the genome we see the variable sites. Here, we are providing you with a data set that measured the genetic diversity at every site in the SARS-CoV-2 genome. Diversity is measured as the number of observed mutations per codon (events).

To properly visualize this type of data, it requires some R wizardry. You don't need to know the details of how this code works, just press play and focus on interpreting the results of the graph. Individual genes are alternatively coloured red and black for easier visualisation.


``` r
#Import data
diversity <- read.csv("nextstrain_ncov_open_global_6m_diversity_aa.csv")
```

```
## Error in file(file, "rt"): cannot open the connection
```

``` r
diversity2 <- diversity %>% 
  group_by(gene) %>% 
  complete(position = 1:max(position), fill = list(events = 0))
```

```
## Error: object 'diversity' not found
```

``` r
#This is so
diversity3 <- diversity2 %>% 
  
  # Compute chromosome size
  group_by(gene) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(diversity, ., by=c("gene"="gene")) %>%
  
  # Add a cumulative position of each SNP
  arrange(gene, position) %>%
  mutate( BPcum=position+tot)
```

```
## Error: object 'diversity' not found
```

``` r
axisdf = diversity3 %>% group_by(gene) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
```

```
## Error: object 'diversity3' not found
```

``` r
ggplot(diversity3, aes(x=BPcum, y=events)) + 
  geom_point( aes(color=as.factor(gene))) + scale_color_manual(values = rep(c("red", "black"), 7 )) + 
  scale_x_continuous(label = axisdf$gene, breaks= axisdf$center) +
  coord_flip()+
  theme_classic() +
  xlab("Gene") +
  ylab("Diversity") +
  theme(legend.position = "none")
```

```
## Error: object 'diversity3' not found
```

As you evaluate the graph, what do you notice? How do you interpret these results in the context of neutral evolution?

*Your answer goes here.*

------------------------------------------------------------------------

## 4.2. Detecting the Footprint of Selection

Applying the principles of the neutral theory, we can also identify what parts of the viral genome are likely under selection. In the data set below, researchers quantified ⍵ for some of the larger genes in the SARS-CoV-2 genome. [⍵ is defined as dN/dS](https://michitobler.github.io/primer-of-evolution/evolution-of-dna-sequences.html#detecting-signatures-of-selection), so it is the ratio of rate of non-synonymous mutations to the rate of synonymous mutations. Plot the values of ⍵ for the different genes using a bar graph.


``` r
#Import data:
dn.ds <- read.csv("dnds_ex/ncov_dnds.csv")

# Filter data for a specific gene
gene_data <- dn.ds[dn.ds$gene == "orf7a",]

# And plot the ratio for the single gene you selected, with a line corresponding to neutrality (dN/dS = 1):
ggplot(gene_data, aes(x=gene, y=w)) +
    geom_bar(stat="identity") +
    theme_classic() +
    geom_hline(yintercept=1) +
    xlab("Gene") +
    ylab("dN/dS Ratio")+
    ylim(0,2)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

``` r
# You are also able to visualise all ratios at once:

ggplot(dn.ds, aes(x=gene, y=w)) +
    geom_bar(stat="identity") +
    theme_classic() +
    geom_hline(yintercept=1) +
    xlab("Gene") +
    ylab("dN/dS Ratio")+
    ylim(0,2)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-2.png)

## 4.3. Interpretation

Based on the information, what genes do you think are under positive selection? What do the results mean for our understanding of COVID-19?

*Your answer goes here.*

------------------------------------------------------------------------

# 5. Answers

------------------------------------------------------------------------


``` r
# For 2.1:

# - flu virus strains evolve over time, with selection happening on their mechanisms of targeting a specific protein

# - coronaviruses likely evolve in the same general way, backed by the diversity of effects that they have (some cause mild illensses, others severe)

# - considering these, SARS-CoV-2 might have originated as a diverging coronavirus strain that specialised in a specific target protein, with rapid evolution as a consequence of host-pathogen evolutionary arms race driving its infection capability
```


``` r
# For 2.3:

# previous hypothesis seems likely; SARS-CoV-2 strains are very closely related to SARS-CoV-like strains present in bats and pangolins, and appear to have originated in a single human infection event and then radiated from there, likely due to host infectivity optimisation as the virus infected more and more individuals and became more efficient at binding its target protein(s).
```


``` r
# For 3.5.1

# It is evident that most early strains were extinguished relatively early, likely as a result to the combined effect of host immune adaptation and competition from more fit strains. The latter radiated/diversified into many different substrains, and the process kept repeating. Essentially, we see how the 'best' strain/few strains at infecting humans took over most infection cases and reduced all other strains to minimal infectivity success. This suggests powerful selection pressures acting in favour of any 'better' strain.
```


``` r
# For 3.5.2

# There have been two noticeable radiation events that we can follow on the tree, and which seem to correspond to the two jumps in mutation rate. We see that diversification was driven not just by selection on pathogen infectivity alone, but also by an increase in mutation rate - which we would expect to help the virus evolve faster and be more successful in the aforementioned arms race with its host.
```


``` r
# For 4.1:

# We can observe that the highest density of mutations is confined to some specific genes or regions, with others showing a significantly decreased mutation density. We could likely assume that the regions with reduced mutation densities are subject to negative selection, as mutations affecting the structure of the proteins coded by these genes would likely have stronger deleterious effects and decrease the fitness of the virus.
```


``` r
# For 4.3:

# We can see that all the genes are under purifying selection, with the strongest effects on the genes/exons "M", "orf1b_reg1", and "orf1b_reg2". This makes these genes likely candidates for the proteins used by the virus to infect its host.
```


------------------------------------------------------------------------

# 6. Resources

------------------------------------------------------------------------

## 6.1. Data References

Phylogenetic trees and genetic diversity data of SARS-CoV-2 and other betacoronaviruses were retrieved from [Nextstrain](https://nextstrain.org/):

-   [Phylogenetic tree of SARS-CoV-2 and SARS-related betacoronaviruses](https://nextstrain.org/groups/blab/sars-like-cov) (accessed: 9/23/2020)
-   [SARS-CoV-2 phylogenetic tree and diversity data](https://nextstrain.org/ncov/global) (accessed: 9/23/2020)

For additional information on Nextstrain, see:

-   Hadfield, J, Megill C, Bell SM, Huddleston J, Potter B, Callender C, Sagulenko P, Bedford T, Neher RA. 2018. [Nextstrain: real-time tracking of pathogen evolution](https://academic.oup.com/bioinformatics/article/34/23/4121/5001388). *Bioinformatics* 34: 4121-4123.

Data for D~N~/D~S~ ratios were kindly provided by Manuela Sironi and Diego Forni based on the following work:

-   Cagliani R, Forni D,Clerici M, Sironi. 2020. [Coding potential and sequence conservation of SARS-CoV-2 and related animal viruses](https://www.sciencedirect.com/science/article/pii/S1567134820301842). *Infection, Genetics and Evolution* 83: 104353.

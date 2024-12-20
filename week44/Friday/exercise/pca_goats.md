---
title: "goats_evolutionary_thinking_pca"
author: "Bjarke Meyer Pedersen"
format: html
editor: visual
---

This in not a normal R markdown document, but instead a Quarto document. A Quarto file, is much like a normal R markdown, but it has some extra features, which you can read about here.

This exercise is build around a dataset from Colli et al. Genet Sel Evol (2018) 50:58 https://doi.org/10.1186/s12711-018-0422-x and lecture notes from Genomics Boot Camp, an online course by Gábor Mészáros.

## Quarto

Quarto enables you to weave together content and executable code into a finished document. To learn more about Quarto see <https://quarto.org>.

## Population Structure within Goats

Today we will do a principal component analysis (PCA) based on quality-controlled genotype data. And as you probably guessed from the title we will look a goat data today. But before jumping into the exciting world of goats. Let's look at the PCA in general, but try to keep it on a high level of abstraction (read; let's try to avoid the math for now).

## The Principal Component Analysis

Principal Component Analysis (PCA) is a powerful statistical technique employed in population genetics to unravel genetic variation within and among populations. By reducing the complexity of genetic data into a few orthogonal axes, named Principal Components. These axes are then ordered by the amounts of variation explained by that linear combination of variables. The math behind it are, I geuss a bit beyond the scope of this course, but it will be covered in Statistical and Machine Learning in Bioinformatics in the spring.

## The Goats

The PCA itself is a way to visualize complex systems in a simple way. In our case, we want to show relationships between the worldwide goat populations genotyped in the ADAPTmap project. After the quality control, we have 4532 animals left. If we compute genetic distances (with PLINK), we get a matrix of 4532 by 4532 animals, with more than 10 million pairwise combinations. So a "rather" long list to scroll through, let alone make sense of it. Fortunately, we can apply some clever statistical methods to simplify it for us. We will lose some precision and nuances about the relationships between breeds, but in exchange, we can see everything in one figure.

## PLINK installation and merging of files

To process the data, we will use a program called PLINK, which is a free, open-source whole genome association analysis toolset, designed to perform a range of basic, large-scale analyses in a computationally efficient manner. You can read more about it here:\[https://zzz.bwh.harvard.edu/plink/\].

Let us start by assigning the working directory to the folder that contains the data

```{r}
setwd('../pca_goat_data_et/') #put your own path here
```

We will also have to install PLINK in our system, if you have conda or mamba create a new IF YOU have mamba/conda use ... conda create -n plink -c bioconda plink in the terminal and then activate the environment in the R terminal.

If this fails, you can download the program ''in the normal way'' from this link https://www.cog-genomics.org/plink/

## The .bed -bam .fam

We have different file types, and in bioinformatics there is a lot of different files types... and most of the work we do is actually to jump from and to specific files that each program we want to use requires. Luckily Wikipedia is a great resource in understanding these different file types and formats.

Here is a small introduction:

1.  **BED File (.bed)** A BED file is a format for representing genomic features, widely used in genomics and bioinformatics. It includes information like chromosome, start, and end positions, making it easy to visualize and analyze genomic regions.

    ![](images/bed_format.png)

2.  **BAM File (.bam)** A BAM file is a compressed binary format for storing high-throughput sequencing data. It's essential for read mapping, variant calling, and genetic variation analysis in genomics.

    ![](BAM_format.png)

3.  **FAM File (.fam)** FAM files store family and genetic data, crucial for genetic association studies. They contain pedigree and phenotypic information used in genome-wide association studies (GWAS) to identify genetic variants associated with specific traits or diseases.

    ![](images/fam_format.png)

## Quality control

Now we have our data, we have merged our dataset and are almost ready to jump into the fun bit.

But before we jump into the fun stuff, we will have to do some QC on the data, otherwise how can we trust it?

With SNP data, as ours, we can take some easy steps in ensuring correctness of our results. Let's make some very simple filtering requirements:

-   if more than 10% of an individuals genotype at the SNPS are unknown/missing, we will remove the individual

    **Does this seem reasonable? and reassure yourselves that you know what is meant by genotype at an SNP**

-   If an SNP has more than 10% missing genotypes?

    **What does this mean?**

-   If an SNP has a minor allele frequency (MAF) below 5% we remove it

    **Explain what the MAF means and why would it be reasonable to remove snps with MAF below 5%?**

-   Last but not least we will remove unplaced mitochondrial and sex chromosome snps.

    ```{r}
    system("../../Downloads/plink_mac_20231018/plink --bfile data/snps --mind 0.1 --geno 0.1 --maf 0.05 --chr-set 29 --make-bed --out afterQC")
    ```

\--mind 0.1 removes individuals with more than 10% missing genotypes

\--geno 0.1 removes SNPs with more than 10% missing genotypes

\--maf 0.05 removes SNPs with minor allele frequency of less than 5%

\--autosome removes unplaced, sex chromosome and mitochondrial SNPs (keeps only autosomal SNPs)

**How many SNPS were removed?**

**How many individuals were removed?**

## Perform PCA

There are multiple options how to create a PCA plot (distance matrix type, software) with possible small differences in the final picture

The \--pca option of PLINK is demonstrated here, but other equivalent solutions are possible

```{r}
system("../../Downloads/plink_mac_20231018/plink --bfile ADAPTmap_genotypeTOP_20160222_full --chr-set 29 --autosome --pca --out pcaResult")
```

## Visualize the results

Load data and the very very very important tidyverse

```{r}
library(tidyverse)

# read in result files
eigenValues <- read_delim("pcaResult.eigenval", delim = " ", col_names = F)
eigenVectors <- read_delim("pcaResult.eigenvec", delim = " ", col_names = F)
```

\## Proportion of variation captured by each vector

```{r}
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2) 
```

**Plot how much the different PC explain, the vector above is sorted so pc1 is entry 1 and pc2 entry 2 and so forth,**

## PCA plot

Plot the eigenvectors, X1 in eigenvectors corresponds to

```{r}
ggplot(data = eigenVectors, aes(x = X3, y = X4, lable = X1, col = X1)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = "PCA of the goat populations",
       x = paste0("Principal component 1 (", eigen_percent[1, 1], " %)"),
       y = paste0("Principal component 2 (", eigen_percent[2, 1], " %)"),
       colour = "Goat breeds") +
  geom_text(aes(label=X1),hjust=0, vjust=0)
```

**OPTIONAL** make the plot more identical with the one of the paper by making a region column

The figure above is from [Colli et al. (2018), published in Genetics, Selection, Evolution](https://gsejournal.biomedcentral.com/track/pdf/10.1186/s12711-018-0422-x) (open access), figure 4a.

You see that the general shape of the two graphs is virtually identical, which is kind of expected, as we used the same initial data and same methods. There are a few minor differences though, worth pointing out.

Maybe the most striking one is the difference in explained variance shown on the axes. This might be due to the differences in quality control. If you check out the Methods section of Colli et al. (2018), you will see that they looked at the data from a variety of different angles, not just missingness. Also, this is a good example of what other factors could be taken into account in QC, so it is worthwhile to check it out anyway.

The other issue is that in the paper there are much fewer data points. This is because only the middle point, i.e. one point for each breed is visualized. Additionally, they have avoided the necessity of the legend by putting the breed abbreviations straight in the figure. Still, many of them are overlapping and unreadable, but arguably better than a legend with hundreds of entries.

**What can we learn about goats and their population structure from this PCA?**

**How does the plot change if you use more strict filtering parameters?**

**What evolutionary explanations causes populations structure? and Can you see an overall pattern to the goat populations?**

**Additional information: Because PCAs is such a fast and easy way to gain valuable insights, one should be careful with the interpretations based in PCAs alone:** https://www.nature.com/articles/s41598-022-14395-4

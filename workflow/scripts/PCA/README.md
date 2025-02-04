# Script for making a PCA from mixed-ploidy VCF files

## Introduction

This R script is based on scripts from [Magdalena Bohutínskás Github](https://github.com/mbohutinska/PolyplChapter/tree/main/PartB). I added support for higher and uneven ploidies (up to 6x) and a custom plotting function with canonical colour coding for different ploidies.

## Requirements

R with [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html) and [vcfR](https://cran.r-project.org/web/packages/vcfR/index.html)


## Usage

```
Rscript make_PCA.r input.vcf pop_map.tsv pca_PCs.tsv pca_Eigen.tsv pca_plot.pdf
```

Where:

### Input files:

- input.vcf is your vcf input file

- pop_map.tsv is a tab separated file that assigns individuals to populations and sets their ploidy. The file should look as follows, assuming population1 has three diploid individuals and population2 has three tetraploid individuals:

individual1&lt;TAB&gt;population1&lt;TAB&gt;2

individual2&lt;TAB&gt;population1&lt;TAB&gt;2

individual3&lt;TAB&gt;population1&lt;TAB&gt;2

individual4&lt;TAB&gt;population2&lt;TAB&gt;4

individual5&lt;TAB&gt;population2&lt;TAB&gt;4

individual6&lt;TAB&gt;population2&lt;TAB&gt;4

### Output files:

- pca_PCs.tsv table with the results of the PCA

- pca_Eigen.tsv table with Eigenvalues from the PCA

- pca_plot.pdf plot of the PCA results

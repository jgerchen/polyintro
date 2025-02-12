# Script for estimating Nei's genetic distances from mixed-ploidy VCF files

## Introduction

This R script uses StAMPP to estimate Nei's genetic distances between populations and individuals with variable ploidies. It then uses the resulting distance matrices to build Neighborjoining trees and Neighbornets.
## Requirements

R with [StAMPP](https://cran.r-project.org/web/packages/StAMPP/index.html), [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html), [vcfR](https://cran.r-project.org/web/packages/vcfR/index.html), [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html) and [viridis](https://cran.r-project.org/web/packages/viridis/index.html).


## Usage

```
Rscript neis_distance.r input.vcf pop_map.tsv heatmap_pop.pdf heatmap_ind.pdf nj_pop.pdf nj_ind.pdf neighbornet_pop.pdf neighbornet_ind.pdf nei_ind.phy

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

- heatmap_pop.pdf is a population level heatmap of Neis D

- heatmap_ind.pdf is an individual level heatmap of Neis D

- nj_pop.pdf is a population level neighbor-joining tree based on Neis D

- nj_ind.pdf is an individual level neighbor-joining tree based on Neis D

- neighbornet_pop.pdf is a population level neighbornet based on Neis D

- nj_ind.pdf is an individual level neighbornet based on Neis D

-nei_ind.phy is an individual level pairwise distance matrix based on Neis D in phylip format, which can be loaded into secondary programs, for example [SplitsTree](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)

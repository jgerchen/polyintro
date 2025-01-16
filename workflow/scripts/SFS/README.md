# Scripts for generating and plotting side frequency spectra from mixed-ploidy VCF files
## Introduction
This script generates simple and 2d site frequency spectra (SFS) from VCF files, allowing for variation in ploidy. In addition, it deals with the presence of missing data by subsampling sides per population to a number of haplotypers defined by the user and only sides, counting only sides that are above the haplotype threshold. It can generate regular and folded SFSs.
In addition, R scripts for plotting the results are provided.
## Requirements
SFS script (SFS.py):
python3 (uses only standard library)
plotting script for simple SFS (plot_SFS.R):
R (uses only base R)
plotting script for 2d SFS SFS (plot_SFS2d.R):
R with [viridis](https://cran.r-project.org/web/packages/viridis/index.html) or [viridisLite](https://cran.r-project.org/web/packages/viridisLite/index.html)
##Usage

###SFS script
```
python3 SFS.py --vcf your_vcf.vcf.gz --pop_map pop_map.tsv --out2d out2d.txt --folded_out2d out2d_folded.txt --output out.txt --folded_output folded_out.txt
```
Where:

--vcf is your vcf file

--pop_map is a tab separated file that assigns individuals to populations and that defines to how many haploid genome copies each population should be subsampled to to account for missing data (it only uses the first entry for each population). The file should look as follows:

individual1<TAB>population1<TAB>8

individual2<TAB>population1<TAB>8

individual3<TAB>population1<TAB>8

individual4<TAB>population2<TAB>4

individual5<TAB>population2<TAB>4

individual6<TAB>population2<TAB>4

Let's assume individuals in population1 are tetraploid and individuals in population2 are diploid. In both cases we would subsample our data to two individuals, note that population1 has double the number of haploid genome copies because its tetraploid. This will allow us to have sites with missing data for one individual in each population, sites with more missing data will be ignored.

--out and --out_folded are normal and folded SFS for all populations defined in your pop_map

--out2d and --out2d_folded are normal and folded 2D SFS for all pairs of populations defined in your pop map
##plotting script for simple SFS
```
Rscript plot_SFS.R SFS.txt SFS_plot.pdf
```
Where SFS.txt is your simple SFS (output from SFS.py) and SFS_plot.pdf is your output pdf containing SFS plots for all populations in SFS.txt

##plotting script for 2d SFS
```
Rscript plot_SFS2d.R SFS2d.txt SFS_plot
```
Where SFS2d.txt is your 2d SFS (output from SFS.py) and SFS_plot is the prefix of output pdf files. The script will create individual pdfs of 2d SFS plots for all combinations of populations in SFS2d.txt


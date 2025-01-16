


python3 SFS.py --vcf your_vcf.vcf.gz --pop_map pop_map.tsv --out2d out2d.txt --folded_out2d out2d_folded.txt --output out.txt --folded_output folded_out.txt

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


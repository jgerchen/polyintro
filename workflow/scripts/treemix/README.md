# Script for making TreeMix input files

## Make treemix input

make_treemix.py

"-v", "--vcf", input VCF file (gzipped)

"-o", "--output", output treemix file

"-p", "--populations", population map (format: individual   population, further columns ignored)

Note: Currently the population map must include all individuals in the VCF file. If you don't want to incldue individuals in the TreeMix input file, set their population to xxx .

## Running TreeMix

To run a basic TreeMix analysis you can do 
```
treemix -i treemix_input -root root_population -o treemix_output
```

Here treemix_input is the previously generated treemix input file, root is one of the populations defined in the population map to be used as outgroup and treemix output is the output file. If you want to add migration edges you can add the -m option with the number of migration edges you want to add.

## Plotting Treemix results

You can get the R script for plotting treemix output from the [official treemix repository](https://bitbucket.org/nygcresearch/treemix/src/master/src/plotting_funcs.R) and run it using

```
Rscript -e "source('plotting_funcs.R'); plot_tree(treemix_output)"
```

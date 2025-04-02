# Scripts for running Structure and Clumpping and plotting results

## Make structure input

vcf_to_structure.py

"-v", "--vcf", input VCF file (gzipped)

"-o", "--output", output structure file

"-p", "--populations", population map (format:individual population, further columns ignored)

"-m", "--max_ploidy", highest ploidy in the VCF file

"-s", "--structure_mainparams", output mainparams file for running structure

"-e", "--extraparams", output extraparams file for running structure

## Clumpp structure results

make_clumpp_input.py

"-p", "--parameter_file", output parameter file

"-i", "--ind_file", output ind file

"-s", "--structure_files", common string that can be matched in all structure files

## Plot structure results

plot_structure_ind.R input_file pop_map_file output_image

input_file=CLUMPP output

pop_map_file:

individual1 population ploidy

individual2 population ploidy

individual3 population ploidy

output_image:

-> pdf file

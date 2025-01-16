import argparse
import gzip
import numpy


parser = argparse.ArgumentParser(description="Prune VCF file")
parser.add_argument("-v", "--vcf_input")
parser.add_argument("-o", "--out_vcf")
parser.add_argument("-p", "--population_input")
parser.add_argument("-d", "--min_distance")
parser.add_argument("-f", "--min_freq")

#also look at missing data?
parser.add_argument("-r", "--max_ld")


args=parser.parse_args()
vcf_input=args.vcf_input
population_input=args.population_input
out_file=args.out_file

vcf_file=gzip.open(vcf_input ,'rt')

pop_dictionary={}
ind_dictionary={}
#only use populations you want to include
#estimate rsquared for each pop and take highest rsq? Different rsq for diploids and tetraploids?
with open(population_input) as population_file:
	for pop_line in population_file:
		ind_name=pop_line.strip().split("\t")[0]
		ind_pop=pop_line.strip().split("\t")[1]
		ind_ploidy=pop_line.strip().split("\t")[2]
		if ind_pop not in pop_dictionary:
			pop_dictionary.update()
#Check by scaffold...	
#First check min distance to prev variant,  then check min freq, then calculate rsquares for all pops....

for line in vcf_file:
	if line[0]=="#":
		if line[1]!="#":
			individuals=line.strip().split("\t")[9:]		
	else:
		#zeros=0
		#ones=0
		#missing_data=0
		columns=line.strip().split("\t")
		genotypes=columns[9:]
vcf_file.close()

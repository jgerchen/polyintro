import argparse
import gzip
parser = argparse.ArgumentParser(description="Make treemix input")
parser.add_argument("-v", "--vcf_input")
parser.add_argument("-o", "--out_file")
parser.add_argument("-p", "--population_input")
parser.add_argument("-m", "--max_missing")

args=parser.parse_args()
vcf_input=args.vcf_input
population_input=args.population_input
out_file=args.out_file
max_missing=float(args.max_missing)

vcf_file=gzip.open(vcf_input ,'rt')
population_file=open(population_input)
output_file=open(out_file, 'w')

pop_dictionary={}
count_dictionary={}

for pop_line in population_file:
	ind_name=pop_line.strip().split("\t")[0]
	ind_pop=pop_line.strip().split("\t")[1]
	ind_ploidy=int(pop_line.strip().split("\t")[2])
	pop_dictionary.update({ind_name:ind_pop})
	if ind_pop not in count_dictionary:
		count_dictionary.update({ind_pop:[0, ind_ploidy, 0]})
	else:
		count_dictionary[ind_pop][1]+=ind_ploidy
pop_list=[i for i in count_dictionary]
output_file.write("\t".join(["Chr", "Position"]+pop_list)+"\n")
for line in vcf_file:
	if line[0]=="#":
		if line[1]!="#":
			individuals=line.strip().split("\t")[9:]		
	else:
		columns=line.strip().split("\t")
		contig=columns[0]
		position=columns[1]
		genotypes=columns[9:]
		#reset count dict
		for count_pop in pop_list:
			count_dictionary[count_pop][0]=0
			count_dictionary[count_pop][2]=0
		for genotype_i in range(len(genotypes)):
			individual=individuals[genotype_i]
			if individual in pop_dictionary:
				individual_population=pop_dictionary[individual]
				alleles=genotypes[genotype_i].split(":")[0].split("/")
				if len(alleles)<=1:
					alleles=genotypes[genotype_i].split(":")[0].split("|")
				count_dictionary[individual_population][0]=alleles.count("1")
				count_dictionary[individual_population][2]=alleles.count(".")
		write_site=True
		alt_count=0
		for count_pop_write in pop_list:
			alt_count+=count_dictionary[count_pop_write][0]
			if count_dictionary[count_pop_write][2]>0:
				if float(count_dictionary[count_pop_write][2])/count_dictionary[count_pop_write][1]>max_missing:
					write_site=False
				else:
		#round up missing data
					count_dictionary[count_pop_write][0]=int(round((float(count_dictionary[count_pop_write][0])/(count_dictionary[count_pop_write][1]-count_dictionary[count_pop_write][2]))*count_dictionary[count_pop_write][1], 0))
		if alt_count==0:
			write_site=False
		if write_site==True:
			output_file.write("\t".join([contig, position]+[str(count_dictionary[u][0])+"/"+str(count_dictionary[u][1]) for u in pop_list])+"\n")

output_file.close()
vcf_file.close()
population_file.close()

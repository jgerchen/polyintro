import argparse
import gzip
parser = argparse.ArgumentParser(description="Make treemix input")
parser.add_argument("-v", "--vcf_input")
parser.add_argument("-o", "--out_file")
parser.add_argument("-p", "--population_input")

args=parser.parse_args()
vcf_input=args.vcf_input
population_input=args.population_input
out_file=args.out_file

vcf_file=gzip.open(vcf_input ,'rt')
population_file=open(population_input)
output_file=gzip.open(out_file, 'wt')

pop_dictionary={}

#pop_list=[]
count_dictionary={}
population_output=[]

for pop_line in population_file:
	ind_name=pop_line.strip().split("\t")[0]
	ind_pop=pop_line.strip().split("\t")[1]
	pop_dictionary.update({ind_name:ind_pop})
	if ind_pop not in count_dictionary:
		count_dictionary.update({ind_pop:[0,0]})
		if ind_pop!="xxx":
			population_output.append(ind_pop)
output_file.write(" ".join(population_output)+"\n")
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
		for genotype_i in range(len(genotypes)):
			alleles=genotypes[genotype_i].split(":")[0].split("/")
			if len(alleles)<=1:
				alleles=genotypes[genotype_i].split(":")[0].split("|")
			individual=individuals[genotype_i]
			individual_population=pop_dictionary[individual]
			count_dictionary[individual_population][0]+=alleles.count("0")
			count_dictionary[individual_population][1]+=alleles.count("1")
		write_site=True
		for l_count in count_dictionary:
			if l_count!="xxx" and count_dictionary[l_count][0]==0 and count_dictionary[l_count][1]==0:
				write_site=False
		if write_site==True:
			output_file.write(" ".join([str(count_dictionary[u][0])+","+str(count_dictionary[u][1]) for u in count_dictionary if u!="xxx"])+"\n")
		for l_count in count_dictionary:
			count_dictionary[l_count][0]=0
			count_dictionary[l_count][1]=0
output_file.close()
vcf_file.close()
population_file.close()

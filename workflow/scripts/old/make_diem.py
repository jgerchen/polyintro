import argparse
import gzip
parser = argparse.ArgumentParser(description="Make treemix input")
parser.add_argument("-v", "--vcf_input")
parser.add_argument("-o", "--out_file")
parser.add_argument("-i", "--individual_file")
parser.add_argument("-e", "--excluded_loci")
parser.add_argument("-k", "--kept_loci")
parser.add_argument("-d", "--min_distance")

args=parser.parse_args()
vcf_input=args.vcf_input
individual_file=args.individual_file
out_file=args.out_file
excluded_loci=args.excluded_loci
kept_loci=args.kept_loci
min_dist=int(args.min_distance)

vcf_file=gzip.open(vcf_input ,'rt')
individual_list=open(individual_file)
out_included_loci=open(kept_loci, 'w')
out_excluded_loci=open(excluded_loci, 'w')

include_individuals=set()
incl_count=0
line_count=0
curr_contig=""
prev_position=0

for ind_line in individual_list:
	ind_inp=ind_line.strip().split("\t")[0]
	if ind_inp not in include_individuals:
		include_individuals.add(ind_inp)
individual_list.close()
for line in vcf_file:
	if line[0]=="#":
		if line[1]!="#":
			individuals=line.strip().split("\t")[9:]		
	else:
		#zeros=0
		#ones=0
		#missing_data=0
		line_count+=1
		output_array=["S"]
		columns=line.strip().split("\t")
		contig=columns[0]
		loc=columns[1]
		genotypes=columns[9:]
		for genotype_i in range(len(genotypes)):
			alleles=genotypes[genotype_i].split(":")[0].split("/")
			if len(alleles)<=1:
				alleles=genotypes[genotype_i].split(":")[0].split("|")
			individual=individuals[genotype_i]
			if individual in include_individuals:
				unique_alleles=set(alleles)
				if len(unique_alleles)==1:
					if "0" in unique_alleles:
						output_array.append("0")
					elif "1" in unique_alleles:	
						output_array.append("2")
					else:
						output_array.append("U")
				else:
						output_array.append("1")
		#Test that not all markers are only one letter, and that both homozygotes are contained
		o_array_set=set(output_array[1:])
		if contig!=curr_contig:
			if curr_contig!="":
				output_file.close()
			output_file=open(out_file+"_"+contig+".out", 'w')
			curr_contig=contig
			prev_position=0
		if len(o_array_set)>1 and "0" in o_array_set and "2" in o_array_set and (int(loc)-prev_position)>=min_dist:
			incl_count+=1
			output_file.write("".join(output_array)+"\n")
			out_included_loci.write("%s\t%s\t%s\n" % (incl_count,contig,loc))
			prev_position=int(loc)
		else:
			out_excluded_loci.write("%s\t%s\t%s\n" % (line_count,contig,loc))
			
output_file.close()
vcf_file.close()
out_excluded_loci.close()
out_included_loci.close()

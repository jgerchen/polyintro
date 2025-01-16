import gzip
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="input VCF file (gzipped)")
parser.add_argument("-o", "--output", help="output structure file")
parser.add_argument("-p", "--populations", help="population map (format:individual<tab>population)")

args=parser.parse_args()
input_gzip=args.vcf
output_structure=args.output
pop_map_file=args.populations

with open(pop_map_file) as pop_map:
	pop_counter=1
	pop_count_dict={}
	pop_dict={}
	for p_line in pop_map:
		p_ind=p_line.strip().split("\t")[0]	
		p_pop=p_line.strip().split("\t")[1]	
		if p_pop not in pop_count_dict:
			pop_count_dict.update({p_pop:str(pop_counter)})
			pop_counter+=1
		pop_dict.update({p_ind:pop_count_dict[p_pop]})

#{i.strip().split("\t")[1]:i.strip().split("\t")[0] for i in pop_map}
ind_list=[]
genotype_list=[]
genotype_dictionary_l1={"./.":"-9", "././.":"-9","./././.":"-9","0/0":"0", "0/1":"1", "1/1":"1", "0/0/0":"0", "0/0/1":"1","0/1/1":"1","1/1/1":"1","0/0/0/0":"0", "0/0/0/1":"1","0/0/1/1":"1","0/1/1/1":"1","1/1/1/1":"1"}
genotype_dictionary_l2={"./.":"-9", "././.":"-9","./././.":"-9", "0/0":"0", "0/1":"0", "1/1":"1", "0/0/0":"0", "0/0/1":"0","0/1/1":"1","1/1/1":"1","0/0/0/0":"0", "0/0/0/1":"0","0/0/1/1":"1","0/1/1/1":"1","1/1/1/1":"1"}
genotype_dictionary_l3={"./.":"-9", "././.":"-9","./././.":"-9", "0/0":"-9", "0/1":"-9", "1/1":"-9", "0/0/0":"0", "0/0/1":"0","0/1/1":"0","1/1/1":"1","0/0/0/0":"0", "0/0/0/1":"0","0/0/1/1":"0","0/1/1/1":"1","1/1/1/1":"1"}
genotype_dictionary_l4={"./.":"-9", "././.":"-9","./././.":"-9", "0/0":"-9", "0/1":"-9", "1/1":"-9", "0/0/0":"-9", "0/0/1":"-9","0/1/1":"-9","1/1/1":"-9","0/0/0/0":"0", "0/0/0/1":"0","0/0/1/1":"0","0/1/1/1":"0","1/1/1/1":"1"}

genotypes_l1=[]
genotypes_l2=[]
genotypes_l3=[]
genotypes_l4=[]

with gzip.open(input_gzip, 'rt') as gzVCF:
	for line in gzVCF:
		if line[0]=="#":
			if line[1]!="#":
				ind_list=line.strip().split("\t")[9:]
		else:
			pos_cats=line.strip().split("\t")
			genotype_list.append(pos_cats[0]+"_"+pos_cats[1])
			genotypes_l1.append([genotype_dictionary_l1["/".join(sorted(i.split(":")[0].replace("|", "/").split("/")))] for i in pos_cats[9:]])
			genotypes_l2.append([genotype_dictionary_l2["/".join(sorted(i.split(":")[0].replace("|", "/").split("/")))] for i in pos_cats[9:]])
			genotypes_l3.append([genotype_dictionary_l3["/".join(sorted(i.split(":")[0].replace("|", "/").split("/")))] for i in pos_cats[9:]])
			genotypes_l4.append([genotype_dictionary_l4["/".join(sorted(i.split(":")[0].replace("|", "/").split("/")))] for i in pos_cats[9:]])
with open(output_structure, "w") as struct_out:
	struct_out.write("\t\t"+"\t".join(genotype_list)+"\n")
	for ind_i in range(len(ind_list)):
		l1_gen_list=[]
		l2_gen_list=[]
		l3_gen_list=[]
		l4_gen_list=[]
		for genotype_i in range(len(genotypes_l1)):
		#line 1:
			l1_gen_list.append(genotypes_l1[genotype_i][ind_i])
		#line 2
			l2_gen_list.append(genotypes_l2[genotype_i][ind_i])
		#line 3
			l3_gen_list.append(genotypes_l3[genotype_i][ind_i])
		#line 4
			l4_gen_list.append(genotypes_l4[genotype_i][ind_i])
		struct_out.write(ind_list[ind_i]+"\t"+pop_dict[ind_list[ind_i]]+"\t"+"\t".join(l1_gen_list)+"\n")
		struct_out.write(ind_list[ind_i]+"\t"+pop_dict[ind_list[ind_i]]+"\t"+"\t".join(l2_gen_list)+"\n")
		struct_out.write(ind_list[ind_i]+"\t"+pop_dict[ind_list[ind_i]]+"\t"+"\t".join(l3_gen_list)+"\n")
		struct_out.write(ind_list[ind_i]+"\t"+pop_dict[ind_list[ind_i]]+"\t"+"\t".join(l4_gen_list)+"\n")

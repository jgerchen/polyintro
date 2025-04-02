import gzip
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="input VCF file (gzipped)")
parser.add_argument("-o", "--output", help="output structure file")
parser.add_argument("-p", "--populations", help="population map (format:individual<tab>population, further columns ignored)")
parser.add_argument("-m", "--max_ploidy", help="highest ploidy in the VCF file")
parser.add_argument("-s", "--structure_mainparams", help="output mainparams file for running structure")
parser.add_argument("-e", "--extraparams", help="output extraparams file for running structure")

args=parser.parse_args()
input_gzip=args.vcf
output_structure=args.output
pop_map_file=args.populations
max_ploidy=int(args.max_ploidy)
mainparams=args.structure_mainparams
extraparams=args.extraparams
ind_list=[]
with open(pop_map_file) as pop_map:
	pop_counter=1
	pop_count_dict={}
	pop_dict={}
	for p_line in pop_map:
		p_ind=p_line.strip().split("\t")[0]	
		ind_list.append(p_ind)
		p_pop=p_line.strip().split("\t")[1]	
		if p_pop not in pop_count_dict:
			pop_count_dict.update({p_pop:str(pop_counter)})
			pop_counter+=1
		pop_dict.update({p_ind:pop_count_dict[p_pop]})

#{i.strip().split("\t")[1]:i.strip().split("\t")[0] for i in pop_map}
genotype_list=[]
#genotype_dictionary_l1={"./.":"-9", "././.":"-9","./././.":"-9","0/0":"0", "0/1":"1", "1/1":"1", "0/0/0":"0", "0/0/1":"1","0/1/1":"1","1/1/1":"1","0/0/0/0":"0", "0/0/0/1":"1","0/0/1/1":"1","0/1/1/1":"1","1/1/1/1":"1"}
#genotype_dictionary_l2={"./.":"-9", "././.":"-9","./././.":"-9", "0/0":"0", "0/1":"0", "1/1":"1", "0/0/0":"0", "0/0/1":"0","0/1/1":"1","1/1/1":"1","0/0/0/0":"0", "0/0/0/1":"0","0/0/1/1":"1","0/1/1/1":"1","1/1/1/1":"1"}
#genotype_dictionary_l3={"./.":"-9", "././.":"-9","./././.":"-9", "0/0":"-9", "0/1":"-9", "1/1":"-9", "0/0/0":"0", "0/0/1":"0","0/1/1":"0","1/1/1":"1","0/0/0/0":"0", "0/0/0/1":"0","0/0/1/1":"0","0/1/1/1":"1","1/1/1/1":"1"}
#genotype_dictionary_l4={"./.":"-9", "././.":"-9","./././.":"-9", "0/0":"-9", "0/1":"-9", "1/1":"-9", "0/0/0":"-9", "0/0/1":"-9","0/1/1":"-9","1/1/1":"-9","0/0/0/0":"0", "0/0/0/1":"0","0/0/1/1":"0","0/1/1/1":"0","1/1/1/1":"1"}

#genotypes_l1=[]
#genotypes_l2=[]
#genotypes_l3=[]
#genotypes_l4=[]

with gzip.open(input_gzip, 'rt') as gzVCF:
	for line in gzVCF:
		if line[0]=="#":
			if line[1]!="#":
				ind_cats=line.strip().split("\t")[9:]
				#every ind gets list of lists[["1","1","-9","-9"] ,...]
				ind_genotypes=[[] for i in ind_cats]
				#TODO: make dict instead!! -> ind name to column i
				ind_names={ind_cats[i]:i for i in range(len(ind_cats))}
				
		else:
			pos_cats=line.strip().split("\t")
			pos_gens=pos_cats[9:]
			genotype_list.append(pos_cats[0]+"_"+pos_cats[1])
			for gen_i in range(len(pos_gens)):
				gen_list=sorted(pos_gens[gen_i].split(":")[0].replace("|", "/").replace(".","-9").split("/"))
				if len(gen_list)<max_ploidy:
					gen_list+=["-9"]*(max_ploidy-len(gen_list))
				ind_genotypes[gen_i].append(gen_list)
print(ind_names)
with open(output_structure, "w") as struct_out:
	struct_out.write("\t\t"+"\t".join(genotype_list)+"\n")
	for ind_i in range(len(ind_list)):
		for gen_i in range(max_ploidy):
			struct_out.write(ind_list[ind_i]+"\t"+pop_dict[ind_list[ind_i]]+"\t"+"\t".join([u[gen_i] for u in ind_genotypes[ind_names[ind_list[ind_i]]]])+"\n")
with open(mainparams, "w") as mainparams_out:
	mainparams_out.write("#define PLOIDY %s\n" % max_ploidy)
	mainparams_out.write("#define NUMINDS %s\n" % len(ind_list))
	mainparams_out.write("#define NUMLOCI %s\n" % len(genotype_list))
	mainparams_out.write("#define INFILE %s\n" % output_structure.split("/")[-1])
	mainparams_out.write("#define OUTFILE %s.out\n" % output_structure.split("/")[-1].split(".")[0])
	mainparams_out.write("#define MISSING -9\n")
	mainparams_out.write("#define ONEROWPERIND 0\n")
	mainparams_out.write("#define LABEL 1\n")
	mainparams_out.write("#define POPDATA 1\n")
	mainparams_out.write("#define POPFLAG 0\n")
	mainparams_out.write("#define LOCDATA 0\n")
	mainparams_out.write("#define PHENOTYPE 0\n")
	mainparams_out.write("#define EXTRACOLS 0\n")
	mainparams_out.write("#define MARKERNAMES 1\n")
	mainparams_out.write("#define RECESSIVEALLELES 0\n")
	mainparams_out.write("#define MAPDISTANCES 0\n")
	mainparams_out.write("#define PHASED 0\n")
	mainparams_out.write("#define PHASEINFO 0\n")
	mainparams_out.write("#define MARKOVPHASE 0\n")
	mainparams_out.write("#define NOTAMBIGUOUS -999\n")
	mainparams_out.write("#define MAXPOPS 2\n")
	mainparams_out.write("#define BURNIN 10000\n")
	mainparams_out.write("#define NUMREPS 100000\n")
with open(extraparams, "w") as extraparams_out:
	extraparams_out.write("#define NOADMIX 0\n")
	extraparams_out.write("#define LINKAGE 0\n")
	extraparams_out.write("#define USEPOPINFO 0\n")
	extraparams_out.write("#define LOCPRIOR 0\n")
	extraparams_out.write("#define FREQSCORR 1\n")
	extraparams_out.write("#define ONEFST 0\n")
	extraparams_out.write("#define INFERALPHA 1\n")
	extraparams_out.write("#define POPALPHAS 0\n")
	extraparams_out.write("#define ALPHA 1.0\n")
	extraparams_out.write("#define INFERLAMBDA 0\n")
	extraparams_out.write("#define POPSPECIFICLAMBDA 0\n")
	extraparams_out.write("#define LAMBDA 1.0\n")
	extraparams_out.write("#define FPRIORMEAN 0.01\n")
	extraparams_out.write("#define FPRIORSD 0.05\n")
	extraparams_out.write("#define UNIFPRIORALPHA 1\n")
	extraparams_out.write("#define ALPHAMAX 10.0\n")
	extraparams_out.write("#define ALPHAPRIORA 1.0\n")
	extraparams_out.write("#define ALPHAPRIORB 2.0\n")
	extraparams_out.write("#define LOG10RMIN -4.0 \n")
	extraparams_out.write("#define LOG10RMAX 1.0\n")
	extraparams_out.write("#define LOG10RPROPSD 0.1\n")
	extraparams_out.write("#define LOG10RSTART -2.0\n")
	extraparams_out.write("#define GENSBACK 2\n")
	extraparams_out.write("#define MIGRPRIOR 0.01\n")
	extraparams_out.write("#define PFROMPOPFLAGONLY 0\n")
	extraparams_out.write("#define LOCISPOP 1\n")
	extraparams_out.write("#define LOCPRIORINIT 1.0\n")
	extraparams_out.write("#define MAXLOCPRIOR 20.0\n")
	extraparams_out.write("#define PRINTNET 1\n")
	extraparams_out.write("#define PRINTLAMBDA 1\n")
	extraparams_out.write("#define PRINTQSUM 1\n")
	extraparams_out.write("#define SITEBYSITE 0\n")
	extraparams_out.write("#define PRINTQHAT 0\n")
	extraparams_out.write("#define UPDATEFREQ 100\n")
	extraparams_out.write("#define PRINTLIKES 0\n")
	extraparams_out.write("#define INTERMEDSAVE 0\n")
	extraparams_out.write("#define ECHODATA 1\n")
	extraparams_out.write("#define ANCESTDIST 0\n")
	extraparams_out.write("#define NUMBOXES 1000\n")
	extraparams_out.write("#define ANCESTPINT 0.90\n")
	extraparams_out.write("#define COMPUTEPROB 1\n")
	extraparams_out.write("#define ADMBURNIN 500\n")
	extraparams_out.write("#define ALPHAPROPSD 0.025\n")
	extraparams_out.write("#define STARTATPOPINFO 0\n")
	extraparams_out.write("#define RANDOMIZE 1\n")
	extraparams_out.write("#define SEED 2245\n")
	extraparams_out.write("#define METROFREQ 10\n")
	extraparams_out.write("#define REPORTHITRATE 0\n")

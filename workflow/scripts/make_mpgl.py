import argparse
import gzip
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="input VCF file")
parser.add_argument("-p", "--ploidy", help="list of sample ploidies (Must be the same number as individuals in input file)")
parser.add_argument("-o", "--outfile", help="output file in mpgl format")
args = parser.parse_args()



input_vcf=args.vcf
input_ploidies=args.ploidy
output_file=args.outfile

with open(input_ploidies) as ploidy_file:
	i_ploidies=[int(line.strip()) for line in ploidy_file]
individuals=[]
#i_ploidies=[]

genotypes=[]
gen_chroms=[]
#record_ploidies=True

with gzip.open(input_vcf, 'rt') as i_vcf:
	for vcf_line in i_vcf:
		if vcf_line[0]=="#":
			if vcf_line[1]!="#":
				individuals=vcf_line.strip().split("\t")[9:]
				assert len(i_ploidies)==len(individuals),"Error: ploidies and individuals must have the same length!"
		else:
			gen_cats=vcf_line.strip().split("\t")
			gen_chroms.append(gen_cats[0])
			gen_indexes=gen_cats[8].split(":")
			
			gen_PL_index=gen_indexes.index("PL")
			gen_gts=gen_cats[9:]
			genotypes.append([])
			for gen_gt_i in range(len(gen_gts)):
				gen_gt_cats=gen_gts[gen_gt_i].split(":")
				if gen_PL_index>len(gen_gt_cats):
					gen_gt_PL="."
				else:
					gen_gt_PL=gen_gt_cats[gen_PL_index]
				if gen_gt_PL==".":
					for no_PT in range(i_ploidies[gen_gt_i]+1):
						genotypes[-1].append("0")
				else:
					gen_gt_PLs=gen_gt_PL.split(",")
					for PL in gen_gt_PLs:
						genotypes[-1].append(PL)
#				if record_ploidies==True:
#					i_ploidies.append(len(gen_gt_PLs))
#			if record_ploidies==True:
#				record_ploidies=False
#for i in range(len(i_ploidies)):
#	print("%s %s" % (individuals[i], i_ploidies[i]))
with open(output_file, "w") as mpgl_out:
	mpgl_out.write(str(len(individuals))+" "+str(len(genotypes))+"\n"+" ".join(individuals)+"\n")
	for genotype_i in range(len(genotypes)):
		mpgl_out.write(gen_chroms[genotype_i]+" "+" ".join(genotypes[genotype_i])+"\n")

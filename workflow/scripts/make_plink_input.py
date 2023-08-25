import gzip
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="input VCF file (gzipped)")
parser.add_argument("-o", "--output", help="output VCF file (gzipped)")
parser.add_argument("-i", "--include", help="list of individuals to include")

args=parser.parse_args()
input_gzip=args.vcf
output_gzip=args.output
input_keep_inds=args.include

with open(input_keep_inds) as keep_inds:
	kinds_set=set([i.strip() for i in keep_inds])
keepind_i=[]

with gzip.open(output_gzip, 'wt') as out_gzVCF:
	with gzip.open(input_gzip, 'rt') as gzVCF:
		for line in gzVCF:
			if line[0]=="#":
				if line[1]=="#":
					out_gzVCF.write(line)
				else:
					all_cats=line.strip().split("\t")
					ind_cats=all_cats[9:]
					for i_cat in range(len(ind_cats)):
						if ind_cats[i_cat] in kinds_set:
							keepind_i.append(i_cat)
					out_gzVCF.write("\t".join(all_cats[0:9])+"\t"+"\t".join([ind_cats[i] for i in keepind_i])+"\n")
			else:
				pos_cats=line.strip().split("\t")
				out_gzVCF.write("scaf_"+pos_cats[0]+"\t"+pos_cats[1]+"\t"+pos_cats[0]+"_"+pos_cats[1]+"\t"+"\t".join(pos_cats[3:9])+"\t"+"\t".join([pos_cats[9:][i] for i in keepind_i])+"\n")
				
			

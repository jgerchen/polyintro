import gzip
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="input VCF file (gzipped)")
parser.add_argument("-o", "--output", help="output VCF file (gzipped)")
parser.add_argument("-i", "--include", help="list of loci to include (format: contig_position)")

args=parser.parse_args()
input_gzip=args.vcf
output_gzip=args.output
input_keep_loci=args.include

with open(input_keep_loci) as keep_loci:
	kloci_set=set([i.strip() for i in keep_loci])

with gzip.open(output_gzip, 'wt') as out_gzVCF:
	with gzip.open(input_gzip, 'rt') as gzVCF:
		for line in gzVCF:
			if line[0]=="#":
				out_gzVCF.write(line)
			else:
				pos_cats=line.strip().split("\t")
				if pos_cats[0]+"_"+pos_cats[1] in kloci_set:
					out_gzVCF.write(line)
				
			

import gzip
import argparse
import glob
parser=argparse.ArgumentParser()
parser.add_argument("-p", "--parameter_file", help="output parameter file")
parser.add_argument("-i", "--ind_file", help="output ind file")
parser.add_argument("-s", "--structure_files", help="common string that can be matched in all structure files")

args=parser.parse_args()
parameter_file=args.parameter_file
ind_file=args.ind_file
structure_files=args.structure_files

#glob output files
struct_glob=glob.glob(structure_files+"*")
n_reps=len(struct_glob)
k_values=[]
n_inds=[]
with open(ind_file,"w") as ind_out:
	for struct_output in struct_glob:
		with open(struct_output) as struct_out:
			write_to_indfile=False
			for struct_line in struct_out:
				if write_to_indfile==True:
					if len(struct_line.strip().split())==0:
						write_to_indfile=False
					elif struct_line.strip().split()[0]!="Label":
						ind_counter+=1
						ind_out.write(str(ind_counter)+"\t"+str(ind_counter)+"\t"+"\t".join(struct_line.strip().split()[2:])+"\n")
						if column_counter==0:
							column_counter=len(struct_line.strip().split())
				if struct_line.strip()=="Inferred ancestry of individuals:":
					write_to_indfile=True
					ind_counter=0
					column_counter=0
			ind_out.write("\n")
			k_values.append(column_counter-5)
			n_inds.append(ind_counter)
assert len(set(k_values))==1, "Error: different k values in input files."
assert len(set(n_inds))==1, "Error: different number of individuals in input files."

with open(parameter_file, "w") as parameter_out:
	parameter_out.write("DATATYPE 0\n")
	parameter_out.write("INDFILE %s\n" % ind_file)
	parameter_out.write("OUTFILE %s.out\n" % structure_files)
	parameter_out.write("MISCFILE %s.misc\n" % structure_files)
	parameter_out.write("K %s\n" % k_values[0])
	parameter_out.write("C %s\n" % n_inds[0])
	parameter_out.write("R %s\n" % n_reps)
	parameter_out.write("M 1\n")
	parameter_out.write("W 1\n")
	parameter_out.write("S 1\n")
	parameter_out.write("OVERRIDE_WARNINGS 1\n")
	parameter_out.write("ORDER_BY_RUN 0\n")
	parameter_out.write("PRINT_EVERY_PERM 0\n")
	parameter_out.write("PRINT_PERMUTED_DATA 0\n")

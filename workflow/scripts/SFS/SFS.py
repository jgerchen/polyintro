import argparse
import gzip
import random
import itertools

parser = argparse.ArgumentParser(description='Estimate SFS from VCF file.')
parser.add_argument('-v', '--vcf')
parser.add_argument('-p', '--pop_map')
parser.add_argument('-d', '--out2d')
parser.add_argument('-e', '--folded_out2d')
parser.add_argument('-o', '--output')
parser.add_argument('-f', '--folded_output')

args = parser.parse_args()

sample_list=[]
vcf_input=args.vcf
pop_map_input=args.pop_map
output=args.output
folded_output=args.folded_output
out2d=args.out2d
folded_out2d=args.folded_out2d
down_sample=True

with open(pop_map_input) as pop_map:
    pop_dict={}
    ind_pop_dict={}
    for pop_map_entry in pop_map:
        pop_map_entry_ind=pop_map_entry.strip().split("\t")[0]
        pop_map_entry_pop=pop_map_entry.strip().split("\t")[1]
        if pop_map_entry_pop not in pop_dict:
            pop_map_entry_subsample=int(pop_map_entry.strip().split("\t")[2])
            pop_dict.update({pop_map_entry_pop:pop_map_entry_subsample})
        #else:
        #    pop_dict[pop_map_entry_pop].append(pop_map_entry_ind)
        ind_pop_dict.update({pop_map_entry_ind:pop_map_entry_pop})
pop_SFS_dict={}
pop_SFS_dict_folded={}
for pop_i in pop_dict:
    pop_SFS_dict.update({pop_i:[0]*(pop_dict[pop_i]+1)})
    pop_SFS_dict_folded.update({pop_i:[0]*(pop_dict[pop_i]+1)})

pop_pairwise=[i_pw for i_pw in itertools.combinations(list(pop_dict),2)]
pop_SFS_dict_pairwise={pop_pair[0]+"_"+pop_pair[1]:[[0]*(pop_dict[pop_pair[1]]+1) for u in range(pop_dict[pop_pair[0]]+1)] for pop_pair in pop_pairwise}
pop_SFS_dict_pairwise_folded={pop_pair[0]+"_"+pop_pair[1]:[[0]*(int(pop_dict[pop_pair[1]]/2)+1) for u in range(int(pop_dict[pop_pair[0]]/2)+1)] for pop_pair in pop_pairwise}
if vcf_input[-3:]==".gz":
    inp_file=gzip.open(vcf_input, 'rt')
else:
    inp_file=open(vcf_input)

for line in inp_file:
    if line[0]=="#":
        if line[1]!="#":
            sample_list=line.strip().split("\t")[9:]
#with open("pop.map", "w") as pop_map_out:
#    for sample in sample_list:
#        pop_map_out.write("%s\t%s\n" % (sample, sample.split("_")[0]))
    else:

        loc_snp_cats=line.strip().split("\t")[9:]
        loc_snp_genotypes=[i.split(":")[0] for i in loc_snp_cats]
        #loc_pop_count_dict={i:0 for i in pop_dict}
        pop_gen_dict={u:[] for u in pop_dict}
        for sample_i in range(len(sample_list)):
            sample_alleles=loc_snp_genotypes[sample_i].split("/")
            if len(sample_alleles)==1:
                sample_alleles=loc_snp_genotypes[sample_i].split("|")
            if sample_list[sample_i] in ind_pop_dict:
                pop_gen_dict[ind_pop_dict[sample_list[sample_i]]]+=sample_alleles
        for locus_pop in pop_gen_dict:
            keep_site=True
            if "." in pop_gen_dict[locus_pop]:
                pop_gen_dict[locus_pop]=[i for i in pop_gen_dict[locus_pop] if i!="."]
            #downsample chromosomes
            if down_sample==True:
                if len(pop_gen_dict[locus_pop])<pop_dict[locus_pop]:
                    keep_site=False
                else:
                    pop_gen_dict[locus_pop]=random.sample(pop_gen_dict[locus_pop], pop_dict[locus_pop])
            #make regular SFS
            if keep_site==True:
#            for keep_locus_pop in pop_gen_dict:
                one_count=pop_gen_dict[locus_pop].count("1")
                pop_SFS_dict[locus_pop][one_count]+=1
                #make folded SFS
                zero_count=pop_gen_dict[locus_pop].count("0")
                if one_count>zero_count:
                    pop_SFS_dict_folded[locus_pop][zero_count]+=1
                    #if zero_count not in pop_SFS_dict_folded[keep_locus_pop]:
                    #    pop_SFS_dict_folded[keep_locus_pop].update({zero_count:1})
                    #else:
                    #    pop_SFS_dict_folded[keep_locus_pop][zero_count]+=1
                else:
                    pop_SFS_dict_folded[locus_pop][one_count]+=1
                    #if one_count not in pop_SFS_dict_folded[keep_locus_pop]:
                    #    pop_SFS_dict_folded[keep_locus_pop].update({one_count:1})
                    #else:
                    #    pop_SFS_dict_folded[keep_locus_pop][one_count]+=1
        for pop_pw_comparison in pop_pairwise:
            #print(pop_pw_comparison)
            #Do we check here for sufficient data or before?
            if len(pop_gen_dict[pop_pw_comparison[0]])==pop_dict[pop_pw_comparison[0]] and len(pop_gen_dict[pop_pw_comparison[1]])==pop_dict[pop_pw_comparison[1]]:
                pop_one_one_count=pop_gen_dict[pop_pw_comparison[0]].count("1")
                pop_one_zero_count=pop_gen_dict[pop_pw_comparison[0]].count("0")
                if pop_one_one_count>pop_one_zero_count:
                    pop_one_out_count=pop_one_zero_count
                else:
                    pop_one_out_count=pop_one_one_count

                #print(pop_one_one_count)
                pop_two_one_count=pop_gen_dict[pop_pw_comparison[1]].count("1")
                pop_two_zero_count=pop_gen_dict[pop_pw_comparison[1]].count("0")
                if pop_two_one_count>pop_two_zero_count:
                    pop_two_out_count=pop_two_zero_count
                else:
                    pop_two_out_count=pop_two_one_count
                #print(pop_two_one_count)
                pop_SFS_dict_pairwise[pop_pw_comparison[0]+"_"+pop_pw_comparison[1]][pop_one_one_count][pop_two_one_count]+=1
                pop_SFS_dict_pairwise_folded[pop_pw_comparison[0]+"_"+pop_pw_comparison[1]][pop_one_out_count][pop_two_out_count]+=1

with open(output, "w") as out_file:
    for out_pop in pop_SFS_dict:
        out_file.write(out_pop+"\t"+"\t".join([str(u) for u in pop_SFS_dict[out_pop]])+"\n")
with open(folded_output, "w") as out_file_folded:
    for out_pop_folded in pop_SFS_dict_folded:
        folded_n_zeros=int((len(pop_SFS_dict_folded[out_pop_folded])-1)/2)
        out_file_folded.write(out_pop_folded+"\t"+"\t".join([str(u) for u in pop_SFS_dict_folded[out_pop_folded][:-folded_n_zeros]])+"\n")

with open(out2d, "w") as out_file_paired:
    with open(folded_out2d, "w") as out_file_paired_folded:
        for out_pairwise_combination in pop_pairwise:
            out_file_paired.write("#%s_%s\n" % (out_pairwise_combination[0], out_pairwise_combination[1]))
            out_file_paired_folded.write("#%s_%s\n" % (out_pairwise_combination[0], out_pairwise_combination[1]))
            for paired_line in pop_SFS_dict_pairwise[out_pairwise_combination[0]+"_"+out_pairwise_combination[1]]:
                out_file_paired.write("\t".join([str(p) for p in paired_line])+"\n")
            for paired_line_folded in pop_SFS_dict_pairwise_folded[out_pairwise_combination[0]+"_"+out_pairwise_combination[1]]:
                out_file_paired_folded.write("\t".join([str(f) for f in paired_line_folded])+"\n")

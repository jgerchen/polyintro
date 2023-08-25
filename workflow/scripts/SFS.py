import argparse
import gzip
import random

parser = argparse.ArgumentParser(description='Estimate SFS from VCF file.')
parser.add_argument('-v', '--vcf')
parser.add_argument('-p', '--pop_map')
parser.add_argument('-d', '--down_projection')
parser.add_argument('-o', '--output')
parser.add_argument('-f', '--folded_output')

args = parser.parse_args()

sample_list=[]
vcf_input=args.vcf
pop_map_input=args.pop_map
output=args.output
folded_output=args.folded_output
down_projection=int(args.down_projection)
down_sample=True

with open(pop_map_input) as pop_map:
    pop_dict={}
    ind_pop_dict={}
    for pop_map_entry in pop_map:
        pop_map_entry_ind=pop_map_entry.strip().split("\t")[0]
        pop_map_entry_pop=pop_map_entry.strip().split("\t")[1]
        if pop_map_entry_pop not in pop_dict:
            pop_dict.update({pop_map_entry_pop:[pop_map_entry_ind]})
        else:
            pop_dict[pop_map_entry_pop].append(pop_map_entry_ind)
        ind_pop_dict.update({pop_map_entry_ind:pop_map_entry_pop})
pop_SFS_dict={}
pop_SFS_dict_folded={}
pop_SFS_freq_dict={}
pop_SFS_freq_dict_folded={}
for pop_i in pop_dict:
    pop_SFS_dict.update({pop_i:{}})
    pop_SFS_dict_folded.update({pop_i:{}})
    pop_SFS_freq_dict.update({pop_i:{}})
    pop_SFS_freq_dict_folded.update({pop_i:{}})


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
            pop_gen_dict[ind_pop_dict[sample_list[sample_i]]]+=sample_alleles
        keep_site=True
        for locus_pop in pop_gen_dict:
            if "." in pop_gen_dict[locus_pop]:
                pop_gen_dict[locus_pop].remove(".")
            #downsample chromosomes
            if down_sample==True:
                if len(pop_gen_dict[locus_pop])<down_projection:
                    keep_site=False
                    break
                else:
                    pop_gen_dict[locus_pop]=random.sample(pop_gen_dict[locus_pop], down_projection)
            #make regular SFS
        if keep_site==True:
            for keep_locus_pop in pop_gen_dict:

        if down_sample==True:
            one_count=pop_gen_dict[keep_locus_pop].count("1")
            if one_count not in pop_SFS_dict[keep_locus_pop]:
                pop_SFS_dict[keep_locus_pop].update({one_count:1})
            else:
                pop_SFS_dict[keep_locus_pop][one_count]+=1
            #make folded SFS
            zero_count=pop_gen_dict[keep_locus_pop].count("0")
            if one_count>zero_count:
                if zero_count not in pop_SFS_dict_folded[keep_locus_pop]:
                pop_SFS_dict_folded[keep_locus_pop].update({zero_count:1})
                else:
                pop_SFS_dict_folded[keep_locus_pop][zero_count]+=1
            else:
                if one_count not in pop_SFS_dict_folded[keep_locus_pop]:
                pop_SFS_dict_folded[keep_locus_pop].update({one_count:1})
                else:
                pop_SFS_dict_folded[keep_locus_pop][one_count]+=1
#make naive SFS based on frequencies
        else:
            one_count=pop_gen_dict[keep_locus_pop].count("1")
            zero_count=pop_gen_dict[keep_locus_pop].count("0")
            one_freq=one_count/(one_count+zero_count)
            pop_SFS_freq_dict[keep_locus_pop].append(one_freq)
            if one_freq>0.5:
                pop_SFS_freq_dict_folded[keep_locus_pop].append(1-one_freq)
            else:
                pop_SFS_freq_dict_folded[keep_locus_pop].append(one_freq)


            #make multidimensional SFS
with open(output, "w") as out_file:
    #determine number of columns
    #out_file.write("\t"+"\t".join([str(i) for i in range(0, down_projection+1)])+"\n")
    out_file.write("\t".join([str(i) for i in range(0, down_projection+1)])+"\n")
    for out_pop in pop_SFS_dict:
        out_string=out_pop
        for out_count_cat in range(0, down_projection+1):
            if out_count_cat in pop_SFS_dict[out_pop]:
                out_string+="\t"+str(pop_SFS_dict[out_pop][out_count_cat])
            else:
                out_string+="\t0"
        out_file.write(out_string+"\n")

with open(folded_output, "w") as folded_out_file:
    if down_projection%2==0:
        n_columns=int(down_projection/2)
    else:
        n_columns=int((down_projection-1)/2)
    #folded_out_file.write("\t"+"\t".join([str(i) for i in range(0, n_columns+1)])+"\n")
    folded_out_file.write("\t".join([str(i) for i in range(0, n_columns+1)])+"\n")
    for out_pop_folded in pop_SFS_dict_folded:
        out_string_folded=out_pop_folded
        for out_count_cat_folded in range(0, n_columns+1):
            if out_count_cat_folded in pop_SFS_dict_folded[out_pop_folded]:
                out_string_folded+="\t"+str(pop_SFS_dict_folded[out_pop_folded][out_count_cat_folded])
            else:
                out_string_folded+="\t0"
        folded_out_file.write(out_string_folded+"\n")

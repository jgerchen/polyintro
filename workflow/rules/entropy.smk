configfile: "../config/entropy.yaml"

rule prune_vcf:
	input:
		vcf_input=config["input_vcf"],
		prune_individuals=config["resource_dir"]+"/{species}_prune_vcf_individuals.lst",
		pop_map=config["resource_dir"]+"/{species}_pop_map.tsv"
	output:
		vcf_pruned=config["temp_output"]+"/{species}_ld_pruned.vcf.gz",
		structure_pruned=config["temp_output"]+"/{species}_ld_pruned.struct"
	threads: 1
	resources: 
		mem_mb=8000,
		disk_mb=8000,
		runtime="2:00:00"
	#conda:
	#	"../envs/prune_vcf.yaml"
	log:	config["log_dir"]+"/{species}_prune_vcf.log"
	shell:
		"""
		
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/prune_vcf.sh
		fi
		temp_folder={config[temp_dir]}/prune_vcf_folder
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input.vcf_input} ${{temp_folder}}/vcf_input.vcf.gz
		cp {input.prune_individuals} ${{temp_folder}}/prune_individuals.lst
		cp {input.pop_map} ${{temp_folder}}/pop_map.tsv
		python3 scripts/make_plink_input.py -v ${{temp_folder}}/vcf_input.vcf.gz -i ${{temp_folder}}/prune_individuals.lst -o ${{temp_folder}}/vcf_pruned_plink.vcf.gz >> {log} 2>> {log} 
		plink2 --bad-ld --allow-extra-chr --vcf ${{temp_folder}}/vcf_pruned_plink.vcf.gz --indep-pairwise {config[prune_vcf_window_size]}kb 1 {config[prune_vcf_max_ld]} --out ${{temp_folder}}/ld_pruned >> {log} 2>> {log}
		python3 scripts/subset_vcf.py -v ${{temp_folder}}/vcf_input.vcf.gz -i ${{temp_folder}}/ld_pruned.prune.in -o ${{temp_folder}}/vcf_pruned_final.vcf.gz >> {log} 2>> {log} 
		cp ${{temp_folder}}/vcf_pruned_final.vcf.gz {output.vcf_pruned}
		python3 scripts/vcf_to_structure.py -v ${{temp_folder}}/vcf_pruned_final.vcf.gz -p ${{temp_folder}}/pop_map.tsv -o ${{temp_folder}}/vcf_pruned_final.struct >> {log} 2>>  {log} 
		cp ${{temp_folder}}/vcf_pruned_final.struct {output.structure_pruned}
		"""
#make entropy input first, read ploidy, num inds and loci from here
rule make_entropy_input:
	input:
		vcf_pruned=config["temp_output"]+"/{species}_ld_pruned.vcf.gz",
		ploidy_list="../resources/{species}_ploidies.lst",
		PCA_PCs=config["temp_output"]+"/{species}_pca_pcs.tsv"
	output:
		entropy_input=config["temp_output"]+"/{species}_ld_pruned.mpgl",
		entropy_qk_files=expand("{TEMP_OUT}/{{species}}_qk{K}",TEMP_OUT=config["temp_output"], K=range(2,int(config["max_K"])+1))
	#conda:
	#	"../envs/entropy_input.yaml"
	log:	config["log_dir"]+"/{species}_entropy_input.log"
	threads: 1
	resources: 
		mem_mb=2000,
		disk_mb=4000,
		runtime="1:00:00"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/entropy_input.sh
		fi
		temp_folder={config[temp_dir]}/entropy_input
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input.vcf_pruned} ${{temp_folder}}/input_VCF 
		cp {input.ploidy_list} ${{temp_folder}}/input_ploidy_list
		python3 scripts/make_mpgl.py -v ${{temp_folder}}/input_VCF -p ${{temp_folder}}/input_ploidy_list -o ${{temp_folder}}/entropy.mpgl >> {log} 2>> {log}
		mv ${{temp_folder}}/entropy.mpgl {output.entropy_input}
		cp {input.PCA_PCs} ${{temp_folder}}/PCA_PCs.tsv		
		#cd ${{temp_folder}}
		Rscript scripts/k_means_ldak.r ${{temp_folder}}/PCA_PCs.tsv {config[max_K]} ${{temp_folder}}/{wildcards.species} >> {log} 2>> {log}
		mv ${{temp_folder}}/{wildcards.species}_qk* {config[temp_output]}
		"""
rule run_structure:
	input:
		entropy_input=config["temp_output"]+"/{species}_ld_pruned.mpgl",
		structure_pruned=config["temp_output"]+"/{species}_ld_pruned.struct"
	output:
		structure_output=config["temp_output"]+"/{species}_{struct_k}_{struct_rep}_structure.out"
	threads: 1
	resources: 
		mem_mb=2000,
		disk_mb=4000,
		runtime="4:00:00"
	#conda:
	#	"../envs/structure.yaml"
	log:	config["log_dir"]+"/{species}_k{struct_k}_rep{struct_rep}_structure.log"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/structure.sh
		fi
		temp_folder={config[temp_dir]}/structure_{wildcards.species}_{wildcards.struct_k}_{wildcards.struct_rep}
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		NINDS=$(head -n 1 {input.entropy_input} | cut -d \" \" -f 1)
		NLOCI=$(head -n 1 {input.entropy_input} | cut -d \" \" -f 2)
		cp {input.structure_pruned} ${{temp_folder}}/structure_input		
		cp ../resources/structure_mainparams ${{temp_folder}}/mainparams
		echo \"#define MAXPOPS {wildcards.struct_k}\" >> ${{temp_folder}}/mainparams
		echo \"#define BURNIN {config[structure_burnin]}\" >> ${{temp_folder}}/mainparams
		echo \"#define NUMREPS {config[structure_numreps]}\" >> ${{temp_folder}}/mainparams
		echo \"#define INFILE ${{temp_folder}}/structure_input\" >> ${{temp_folder}}/mainparams
		echo \"#define OUTFILE ${{temp_folder}}/structure_output\" >> ${{temp_folder}}/mainparams
		echo \"#define PLOIDY {config[max_ploidy]}\" >> ${{temp_folder}}/mainparams
		echo \"#define NUMINDS ${{NINDS}}\" >> ${{temp_folder}}/mainparams
		echo \"#define NUMLOCI ${{NLOCI}}\" >> ${{temp_folder}}/mainparams
		cp ../resources/structure_extraparams ${{temp_folder}}/extraparams
		structure -i ${{temp_folder}}/structure_input -o ${{temp_folder}}/structure_output -m ${{temp_folder}}/mainparams -e ${{temp_folder}}/extraparams >> {log} 2>>{log}	
		cp ${{temp_folder}}/structure_output_f {output} 
		"""
rule clumpp_structure_results:
	input:
		structure_output=expand("{TEMP_OUT}/{{species}}_{{struct_k}}_{struct_rep}_structure.out", TEMP_OUT=config["temp_output"], struct_rep=range(1, int(config["structure_replicate_runs"])+1)),
		entropy_input=config["temp_output"]+"/{species}_ld_pruned.mpgl"
	output:
		clumpp_output=config["temp_output"]+"/{species}_{struct_k}.clumpp",
		clumpp_misc_out=config["temp_output"]+"/{species}_{struct_k}_misc.clumpp"
	threads: 1
	resources: 
		mem_mb=2000,
		disk_mb=4000,
		runtime="1:00:00"
	log:	config["log_dir"]+"/{species}_k{struct_k}_clumpp.log"
	#conda:
	#	"../envs/structure.yaml"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/structure.sh
		fi
		NINDS=$(head -n 1 {input.entropy_input} | cut -d \" \" -f 1)
		temp_folder={config[temp_dir]}/clumpp_{wildcards.species}_{wildcards.struct_k}
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input} $temp_folder
		echo "DATATYPE 0" > ${{temp_folder}}/mainparams	
		#make indfile from structure outputs
		for struct_out in $temp_folder/*_structure.out
		do
			grep -A $NINDS \"(%Miss)\" $struct_out | tail -n $NINDS | awk '{{s=\"\"; for (i=6; i<=NF; i++) s=s $i \"\\t\"; print $1\"\\t\"$1\"\\t(0)\\t\"$4\"\\t:\\t\"s}}' >> $temp_folder/clumpp_indfile
			echo -e \"\\n\" >> $temp_folder/clumpp_indfile

		done
		
		echo \"DATATYPE 0\" > ${{temp_folder}}/mainparams	
		echo \"INDFILE ${{temp_folder}}/clumpp_indfile\" >> ${{temp_folder}}/mainparams	
		echo \"OUTFILE ${{temp_folder}}/clumpp_out\">> ${{temp_folder}}/mainparams
		echo \"MISCFILE ${{temp_folder}}/clumpp_misc\" >> ${{temp_folder}}/mainparams
		echo \"K \"{wildcards.struct_k} >> ${{temp_folder}}/mainparams
		echo \"C ${{NINDS}}\" >> ${{temp_folder}}/mainparams
		echo \"R {config[structure_replicate_runs]}\" >> ${{temp_folder}}/mainparams
		echo \"M 1\" >> ${{temp_folder}}/mainparams
		echo \"W 1\" >> ${{temp_folder}}/mainparams
		echo \"S 1\" >> ${{temp_folder}}/mainparams
		echo \"OVERRIDE_WARNINGS 1\" >> ${{temp_folder}}/mainparams
		echo \"ORDER_BY_RUN 0\" >> ${{temp_folder}}/mainparams
		echo \"PRINT_EVERY_PERM 0\" >> ${{temp_folder}}/mainparams
		echo \"PRINT_PERMUTED_DATA 0\" >> ${{temp_folder}}/mainparams
		CLUMPP ${{temp_folder}}/mainparams >> {log} 2>> {log} 
		mv  ${{temp_folder}}/clumpp_out {output.clumpp_output} 
		mv  ${{temp_folder}}/clumpp_misc {output.clumpp_misc_out}
		"""
rule run_entropy:
	input:
		entropy_input=config["temp_output"]+"/{species}_ld_pruned.mpgl",
		entropy_qk_file=config["temp_output"]+"/{species}_qk{entropy_k}",
		ploidy_list="../resources/{species}_ploidies.lst"
	output:
		entropy_out_hdf5=config["temp_output"]+"/{species}_{entropy_k}_{entropy_rep}_entropy.hdf5",
		entropy_out_mcmc_q=config["temp_output"]+"/{species}_{entropy_k}_{entropy_rep}_entropyq.mcmc",
		entropy_out_mcmc_p=config["temp_output"]+"/{species}_{entropy_k}_{entropy_rep}_entropyp.mcmc",
		entropy_out_posterior_q=config["temp_output"]+"/{species}_{entropy_k}_{entropy_rep}_entropypostq.tsv"
	threads: 1
	resources: 
		mem_mb=8000,
		disk_mb=4000,
		runtime="24:00:00"
	log:	config["log_dir"]+"/{species}_k{entropy_k}_rep{entropy_rep}_entropy.log"
	#conda:
	#	"../envs/entropy.yaml"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/entropy.sh
		fi
		temp_folder={config[temp_dir]}/entropy_{wildcards.species}_{wildcards.entropy_k}_{wildcards.entropy_rep}
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input.entropy_input} ${{temp_folder}}/entropy_input.mpgl
		cp {input.entropy_qk_file} ${{temp_folder}}/entropy_qk
		cp {input.ploidy_list} ${{temp_folder}}/ploidy.lst
		entropy -i ${{temp_folder}}/entropy_input.mpgl -w 1 -m 1 -D 1 -n ${{temp_folder}}/ploidy.lst -k {wildcards.entropy_k} -q ${{temp_folder}}/entropy_qk -l {config[entropy_numreps_l]} -b {config[entropy_burnin_b]} -t {config[entropy_store_n_step_t]} -o ${{temp_folder}}/entropy_out.hdf5 >> {log} 2>> {log}
		estpost.entropy -p q -s 2 -o ${{temp_folder}}/entropy_q.mcmc ${{temp_folder}}/entropy_out.hdf5 >> {log} 2>> {log}
		estpost.entropy -p p -s 2 -o ${{temp_folder}}/entropy_p.mcmc ${{temp_folder}}/entropy_out.hdf5 >> {log} 2>> {log}
		estpost.entropy -p q -s 0 -o ${{temp_folder}}/entropy_q.post ${{temp_folder}}/entropy_out.hdf5 >> {log} 2>> {log}
		mv ${{temp_folder}}/entropy_out.hdf5 {output.entropy_out_hdf5}
		mv ${{temp_folder}}/entropy_q.mcmc {output.entropy_out_mcmc_q}
		mv ${{temp_folder}}/entropy_p.mcmc {output.entropy_out_mcmc_p}
		mv ${{temp_folder}}/entropy_q.post {output.entropy_out_posterior_q}
		"""
rule entropy_assess_param_convergence:
	input:
		entropy_hdf5s=expand("{TEMP_OUT}/{{species}}_{{entropy_k}}_{entropy_rep}_entropy.hdf5",TEMP_OUT=config["temp_output"],entropy_rep=range(1,int(config["entropy_replicate_runs"])+1))
	output:
		entropy_rhat_p=config["temp_output"]+"/{species}_{entropy_k}_entropyp.rhat",
		entropy_rhat_q=config["temp_output"]+"/{species}_{entropy_k}_entropyq.rhat",
		entropy_out_posterior_q_all=config["temp_output"]+"/{species}_{entropy_k}_allentropypostq.tsv"
	threads: 1 
	resources: 
		mem_mb=2000,
		disk_mb=4000,
		runtime="1:00:00"
	log:	config["log_dir"]+"/{species}_k{entropy_k}_entropy_convergence.log"
	#conda:
	#	"../envs/entropy.yaml"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/entropy.sh
		fi
		temp_folder={config[temp_dir]}/entropy_param_convergence_{wildcards.species}_{wildcards.entropy_k}
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input.entropy_hdf5s} ${{temp_folder}}
		estpost.entropy -p p -s 4 ${{temp_folder}}/*.hdf5 -o ${{temp_folder}}/entropy_p.rhat
		estpost.entropy -p q -s 4 ${{temp_folder}}/*.hdf5 -o ${{temp_folder}}/entropy_q.rhat
		estpost.entropy -p q -s 0 ${{temp_folder}}/*.hdf5 -o ${{temp_folder}}/entropy_q_all_rep.post
		mv ${{temp_folder}}/entropy_p.rhat {output.entropy_rhat_p}
		mv ${{temp_folder}}/entropy_q.rhat {output.entropy_rhat_q}
		mv ${{temp_folder}}/entropy_q_all_rep.post {output.entropy_out_posterior_q_all}
		"""

rule plot_entropy:
	input:
		entropy_out_posterior_q_all=expand("{TEMP_OUT}/{{species}}_{entropy_k}_allentropypostq.tsv",TEMP_OUT=config["temp_output"], entropy_k=range(2,config["max_K"]+1)),
		ploidy_list="../resources/{species}_ploidies.lst",
		pop_map="../resources/{species}_pop_map.tsv"
	output:
		entropy_out_posterior_q_all_pdf=expand("../results/entropy/{{species}}_{entropy_k}_allentropypostq.pdf",TEMP_OUT=config["temp_output"], entropy_k=range(2,config["max_K"]+1)),
		entropy_rdata="../results/entropy/{species}_entropy.Rdata"
	threads: 1 
	resources: 
		mem_mb=8000,
		disk_mb=4000,
		runtime="1:00:00"
	log:	config["log_dir"]+"/{species}_plot_entropy.log"
	#conda:
	#	"../envs/plot_rmd.yaml"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/plot_rmd.sh
		fi
		temp_folder={config[temp_dir]}/plot_entropy_{wildcards.species}
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input} $temp_folder
		Rscript scripts/plot_entropy.R $temp_folder/{wildcards.species}_ploidies.lst $temp_folder/{wildcards.species}_pop_map.tsv $temp_folder/entropy_plot.Rdata $temp_folder/{wildcards.species}_{{2..{config[max_K]}}}_allentropypostq.tsv $temp_folder/{wildcards.species}_{{2..{config[max_K]}}}_allentropypostq.pdf >> {log} 2>> {log}
		mkdir -p ../results/entropy
		mv $temp_folder/entropy_plot.Rdata {output.entropy_rdata}
		mv $temp_folder/*.pdf ../results/entropy
		"""

rule plot_structure:
	input:
		structure_out_clumpp_all=expand("{TEMP_OUT}/{{species}}_{struct_k}.clumpp",TEMP_OUT=config["temp_output"], struct_k=range(2,config["max_K"]+1)),
		ploidy_list="../resources/{species}_ploidies.lst",
		pop_map="../resources/{species}_pop_map.tsv"
	output:
		structure_plot_clumpp_all=expand("../results/structure/{{species}}_{struct_k}_allstructureq.pdf", struct_k=range(2,config["max_K"]+1)),
		structure_rdata="../results/structure/{species}_structure.Rdata"
	threads: 1 
	resources: 
		mem_mb=8000,
		disk_mb=4000,
		runtime="1:00:00"
	log:	config["log_dir"]+"/{species}_plot_structure.log"
	#conda:
	#	"../envs/plot_rmd.yaml"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/plot_rmd.sh
		fi
		temp_folder={config[temp_dir]}/plot_structure_{wildcards.species}
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input} $temp_folder
		Rscript scripts/plot_structure.R $temp_folder/{wildcards.species}_ploidies.lst $temp_folder/{wildcards.species}_pop_map.tsv $temp_folder/structure_plot.Rdata $temp_folder/{wildcards.species}_{{2..{config[max_K]}}}.clumpp $temp_folder/{wildcards.species}_{{2..{config[max_K]}}}_allstructureq.pdf >> {log} 2>> {log}
		mkdir -p ../results/structure
		mv $temp_folder/structure_plot.Rdata {output.structure_rdata}
		mv $temp_folder/*.pdf ../results/structure
		"""


rule make_PCA:
	input:
		vcf_input=config["temp_output"]+"/{species}_ld_pruned.vcf.gz"
		#vcf_input=config["input_vcf"]
	output:
		PCA_PCs=config["temp_output"]+"/{species}_pca_pcs.tsv",
		PCA_eigenvalues=config["temp_output"]+"/{species}_pca_eigen.tsv",
		PCA_plot=config["result_dir"]+"/{species}_PCA.pdf"
	threads: 1
	resources:
		mem_mb=20000,
		disk_mb=20000,
		runtime="2:00:00"
	#conda:
	#	"../envs/r_popgen.yaml"
	log:	config["log_dir"]+"/{species}_PCA.log"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/r_popgen.sh
		fi
		temp_folder={config[temp_dir]}/{wildcards.species}_SFS
		trap 'rm -rf $temp_folder' TERM EXIT
		mkdir -p $temp_folder
		cp {input.vcf_input} $temp_folder/{wildcards.species}.vcf.gz
		Rscript scripts/make_PCA.r $temp_folder/{wildcards.species}.vcf.gz $temp_folder/{wildcards.species}_pca_pcs.tsv $temp_folder/{wildcards.species}_pca_eigen.tsv $temp_folder/{wildcards.species}_PCA.pdf >> {log} 2>> {log}
		cp $temp_folder/{wildcards.species}_pca_pcs.tsv {output.PCA_PCs}
		cp $temp_folder/{wildcards.species}_pca_eigen.tsv {output.PCA_eigenvalues}
		cp $temp_folder/{wildcards.species}_PCA.pdf {output.PCA_plot}
		"""

rule make_SFS:
	input:
		vcf_input=config["input_vcf"],
		pop_map=config["resource_dir"]+"/{species}_pop_map.tsv"
	output:
		SFS_table=config["temp_output"]+"/{species}_SFS.tsv",
		Folded_SFS_table=config["temp_output"]+"/{species}_SFS_folded.tsv"
	threads: 1
	resources:
		mem_mb=8000,
		disk_mb=10000,
		runtime="2:00:00"
	#conda:
	#	"../envs/r_popgen.yaml"
	log:	config["log_dir"]+"/{species}_SFS.log"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/r_popgen.sh
		fi
		temp_folder={config[temp_dir]}/{wildcards.species}_SFS
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input.vcf_input} $temp_folder
		cp {input.pop_map} $temp_folder
		i_vcf=$(ls $temp_folder/*.vcf*)
		python3 scripts/SFS.py -v $i_vcf -p $temp_folder/{wildcards.species}_pop_map.tsv -d {config[SFS_projection]} -o $temp_folder/{wildcards.species}_SFS.tsv -f $temp_folder/{wildcards.species}_SFS_folded.tsv &>>{log}
		cp $temp_folder/{wildcards.species}_SFS.tsv {output.SFS_table} 
		cp $temp_folder/{wildcards.species}_SFS_folded.tsv {output.Folded_SFS_table}
		"""


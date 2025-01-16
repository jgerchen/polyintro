rule make_PCA:
	input:
		vcf_input=config["temp_output"]+"/{species}_ld.pruned.vcf.gz"
		#vcf_input=config["input_vcf"]
	output:
		PCA_PCs=config["temp_output"]+"/{species}_pca_pcs.tsv",
		PCA_eigenvalues=config["temp_output"]+"/{species}_pca_eigen.tsv",
		PCA_plot=report(config["result_dir"]+"/{species}_PCA.pdf", category="PCA", labels={"Analysis":"PCA"})
	threads: 1
	resources:
		mem_mb=20000,
		disk_mb=20000,
		runtime="2h"
	params:
		pca_script=workflow.source_path("../scripts/make_PCA.r")
	log:	config["log_dir"]+"/{species}_PCA.log"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/r_popgen.sh
		fi
		temp_folder={config[temp_dir]}/{wildcards.species}_PCA
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {params.pca_script} $temp_folder
		cp {input.vcf_input} $temp_folder/{wildcards.species}.vcf.gz
		Rscript ${{temp_folder}}/make_PCA.r $temp_folder/{wildcards.species}.vcf.gz $temp_folder/{wildcards.species}_pca_pcs.tsv $temp_folder/{wildcards.species}_pca_eigen.tsv $temp_folder/{wildcards.species}_PCA.pdf
		cp $temp_folder/{wildcards.species}_pca_pcs.tsv {output.PCA_PCs}
		cp $temp_folder/{wildcards.species}_pca_eigen.tsv {output.PCA_eigenvalues}
		cp $temp_folder/{wildcards.species}_PCA.pdf {output.PCA_plot}
		"""

rule stampp_neighbornet:
	input:
		vcf_input=config["temp_output"]+"/{species}_ld.pruned.vcf.gz",
		pop_map=config["resource_dir"]+"/{species}_pop_map.tsv",
		ploidy_map=config["resource_dir"]+"/{species}_ploidies.tsv"
		#vcf_input=config["input_vcf"]
	output:
		heatmap_pop=report(config["result_dir"]+"/{species}_heatmap_pop.pdf", category="stampp", labels={"Analysis": "Heatmap pop"}),
		heatmap_ind=report(config["result_dir"]+"/{species}_heatmap_ind.pdf", category="stampp", labels={"Analysis": "Heatmap ind"}),
		nj_pop=report(config["result_dir"]+"/{species}_nj_tree_pop.pdf", category="stampp", labels={"Analysis":"NJ tree pop"}),
		nj_ind=report(config["result_dir"]+"/{species}_nj_tree_ind.pdf", category="stampp", labels={"Analysis":"NJ tree ind"}),
		neighnet_pop=report(config["result_dir"]+"/{species}_neighbornet_pop.pdf", category="stampp", labels={"Analysis":"Neighbornet pop"}),
		neighnet_ind=report(config["result_dir"]+"/{species}_neighbornet_ind.pdf", category="stampp", labels={"Analysis":"Neighbornet ind"}),
		nei_dist_phy=config["result_dir"]+"/{species}_nei_dist_ind.phy"
	threads: 1
	resources:
		mem_mb=20000,
		disk_mb=20000,
		runtime="2h"
	params:
		stampp_script=workflow.source_path("../scripts/neis_distance.r")
	log:	config["log_dir"]+"/{species}_stampp_neighbornet.log"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/r_popgen.sh
		fi
		temp_folder={config[temp_dir]}/{wildcards.species}_stampp_neighnet
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {params.stampp_script} $temp_folder
		cp {input.vcf_input} $temp_folder/{wildcards.species}.vcf.gz
		cp {input.pop_map} $temp_folder/{wildcards.species}_pop_map.tsv
		cp {input.ploidy_map} $temp_folder/{wildcards.species}_ploidy_map.tsv
		Rscript ${{temp_folder}}/neis_distance.r $temp_folder/{wildcards.species}.vcf.gz $temp_folder/{wildcards.species}_ploidy_map.tsv $temp_folder/{wildcards.species}_pop_map.tsv $temp_folder/{wildcards.species}_heatmap_pop.pdf $temp_folder/{wildcards.species}_heatmap_ind.pdf $temp_folder/{wildcards.species}_nj_tree_pop.pdf $temp_folder/{wildcards.species}_nj_tree_ind.pdf $temp_folder/{wildcards.species}_neighbornet_pop.pdf $temp_folder/{wildcards.species}_neighbornet_ind.pdf $temp_folder/{wildcards.species}_nei_dist_ind.phy
		cp $temp_folder/{wildcards.species}_heatmap_pop.pdf {output.heatmap_pop}
		cp $temp_folder/{wildcards.species}_heatmap_ind.pdf {output.heatmap_ind}
		cp $temp_folder/{wildcards.species}_nj_tree_pop.pdf {output.nj_pop}
		cp $temp_folder/{wildcards.species}_nj_tree_ind.pdf {output.nj_ind}
		cp $temp_folder/{wildcards.species}_neighbornet_pop.pdf {output.neighnet_pop}
		cp $temp_folder/{wildcards.species}_neighbornet_ind.pdf {output.neighnet_ind}
		cp $temp_folder/{wildcards.species}_nei_dist_ind.phy {output.nei_dist_phy}
		"""

rule make_SFS:
	input:
		vcf_input=config["input_vcf"],
		pop_map=config["resource_dir"]+"/{species}_pop_map_SFS.tsv"
	output:
		SFS_table=config["temp_output"]+"/{species}_SFS.tsv",
		Folded_SFS_table=config["temp_output"]+"/{species}_SFS_folded.tsv",
		SFS_table_2d=config["temp_output"]+"/{species}_SFS_2d.tsv",
		SFS_table_2d_folded=config["temp_output"]+"/{species}_SFS_2d_folded.tsv"
	threads: 1
	resources:
		mem_mb=8000,
		disk_mb=10000,
		runtime="2h"
	params:
		sfs_script=workflow.source_path("../scripts/SFS.py")
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
		cp {params.sfs_script} $temp_folder
		#make SFS input file
		i_vcf=$(ls $temp_folder/*.vcf*)
		python3 ${{temp_folder}}/SFS.py -v $i_vcf -p $temp_folder/{wildcards.species}_pop_map_SFS.tsv -o $temp_folder/{wildcards.species}_SFS.tsv -f $temp_folder/{wildcards.species}_SFS_folded.tsv -d $temp_folder/{wildcards.species}_SFS_2d.tsv -e $temp_folder/{wildcards.species}_SFS_2d_folded.tsv
		cp $temp_folder/{wildcards.species}_SFS.tsv {output.SFS_table} 
		cp $temp_folder/{wildcards.species}_SFS_2d.tsv {output.SFS_table_2d} 
		cp $temp_folder/{wildcards.species}_SFS_2d_folded.tsv {output.SFS_table_2d_folded} 
		cp $temp_folder/{wildcards.species}_SFS_folded.tsv {output.Folded_SFS_table}
		"""


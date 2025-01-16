configfile: "../config/treemix.yaml"
	
rule make_treemix_input:
	input:
		vcf_input=config["temp_output"]+"/{species}_ld.pruned.vcf.gz",
		pop_map=config["resource_dir"]+"/{species}_pop_map_treemix.tsv"

	output:
		treemix_input=config["temp_output"]+"/{species}_treemix.gz"
	threads: 1
	resources: 
		mem_mb=16000,
		disk_mb=50000,
		runtime="6h"
	params:
		make_treemix_script=workflow.source_path("../scripts/make_treemix.py")
	log:	config["log_dir"]+"/{species}_make_treemix_input.log"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/treemix.sh
		fi
		temp_folder={config[temp_dir]}/treemix_input
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {params.make_treemix_script} $temp_folder
		cp {input.vcf_input} ${{temp_folder}}/vcf_input.vcf.gz
		cp {input.pop_map} ${{temp_folder}}/pop_map.tsv
		python3 ${{temp_folder}}/make_treemix.py -v ${{temp_folder}}/vcf_input.vcf.gz -p ${{temp_folder}}/pop_map.tsv -o ${{temp_folder}}/pruned.gz
		cp ${{temp_folder}}/pruned.gz {output.treemix_input}
		"""

rule run_treemix:
	input:
		treemix_pruned=config["temp_output"]+"/{species}_treemix.gz",
		pop_map=config["resource_dir"]+"/{species}_pop_map_treemix.tsv"
	output:
		treemix_out_plot=report(config["result_dir"]+"/{species}_{migedge}_treemix.pdf", category="treemix", labels={"output":"tree","migration edges": "{migedge}"}),
		treemix_out_dir=directory(config["temp_output"]+"/{species}_{migedge}_treemix"),
		treemix_out_residuals=report(config["result_dir"]+"/{species}_{migedge}_residuals.pdf", category="treemix", labels={"output":"residuals","migration edges": "{migedge}"})
	threads: 1
	resources: 
		mem_mb=8000,
		disk_mb=20000,
		runtime="6h"
	log:	config["log_dir"]+"/{species}_run_treemix_{migedge}.log"
	shell:
		"""
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[prerun_scripts_dir]}/treemix.sh
		fi
		temp_folder={config[temp_dir]}/treemix_run_{wildcards.species}_{wildcards.migedge}
		mkdir -p $temp_folder 		
		trap 'rm -rf $temp_folder' TERM EXIT
		cp {input.treemix_pruned} $temp_folder
		cp {input.pop_map} $temp_folder
		cd $temp_folder
		mkdir treemix_out
		#run treemix
		if [ {wildcards.migedge} -eq 0 ]
		then
			treemix -i $temp_folder/{wildcards.species}_treemix.gz -root {config[treemix_outgroup]} -o $temp_folder/treemix_out/{wildcards.species}_{wildcards.migedge}
		else
			treemix -i $temp_folder/{wildcards.species}_treemix.gz -m {wildcards.migedge} -root {config[treemix_outgroup]} -o $temp_folder/treemix_out/{wildcards.species}_{wildcards.migedge}
		fi
		cp -r $temp_folder/treemix_out {output.treemix_out_dir}
		R_PLOT_FUNC=$(whereis plotting_funcs.R | cut -d \" \" -f 2)
		cp $R_PLOT_FUNC $temp_folder
		awk '{{if ($2!=\"xxx\") print $2}}' $temp_folder/{wildcards.species}_pop_map_treemix.tsv | sort | uniq > $temp_folder/poplist.txt
		Rscript -e \"source('\"$temp_folder\"/plotting_funcs.R'); plot_tree('\"$temp_folder\"/treemix_out/{wildcards.species}_{wildcards.migedge}')\"
		cp Rplots.pdf {output.treemix_out_plot}
		Rscript -e \"source('\"$temp_folder\"/plotting_funcs.R'); plot_resid('\"$temp_folder\"/treemix_out/{wildcards.species}_{wildcards.migedge}', '\"$temp_folder\"/poplist.txt')\"
		cp Rplots.pdf {output.treemix_out_residuals}
		"""



# establish the microhaplot directory with the Shiny components in it.
# (This way we can do this operation only once).  This also ensures
# that microhaplot is installed
rule make_microhap_folder:
	params:
		rd = "{run_dir}/{species_dir}"
	output:
		"{run_dir}/{species_dir}/microhaplot/ui.R",
		"{run_dir}/{species_dir}/microhaplot/server.R"
	log:
		"{run_dir}/{species_dir}/logs/microhap_shiny_stuff/microhap-install-etc.log"
#	conda:
#		"../envs/microhap.yaml"
	envmodules:
		"R/4.0.3"  # this is for SEDNA. I can't get a conda-installed R to work
	script:
		"../script/create_microhaplot_folder.R"




# for each collection of SAM files, make an rds file witin an
# microhaplot directory like this
# {run_dir}/{species_dir}/microhaplot/{marker_set}--{type}--{additional}.rds
# The filename might be, for example:
# LFAR--fullgex_remapped_to_thinned--Otsh_v1.0.rds
#
# or
#
# ROSA--target_fastas--rosawr.rds
#
# The input here is all of the sam files needed to create
# the single output. (This is an aggregation step and
# we want all the dependencies to propagate appropriately.)
# rule: micohaplot
rule microhap_extract_fullgex_remapped:
	input:
		sams = lambda wc: bam_tree_equivalent_files_from_marker_sets(wc, type = "fullgex_remapped", trunk = "sams", ext = ".sam"),
		input_vcf = fullgex_remapped_mh_vcf_from_marker_set_genome_microhap_vcf,
		ui = "{run_dir}/{species_dir}/microhaplot/ui.R",
		shiny = "{run_dir}/{species_dir}/microhaplot/server.R" 
	log:
		"{run_dir}/{species_dir}/logs/microhap_extract_fullgex_remapped/{marker_set}--{genome}--{microhap_variants}.log"
	envmodules:
		"R/4.0.3"  # this is for SEDNA. I can't get a conda-installed R to work
	threads: 10
	output:
		rds="{run_dir}/{species_dir}/microhaplot/{marker_set}--fullgex_remapped_to_thinned--{genome}--{microhap_variants}.rds"
	script:
		"../script/extract_microhaps.R"


# here is the rule for target fastas
rule microhap_extract_target_fastas:
	input:
		sams = lambda wc: bam_tree_equivalent_files_from_marker_sets(wc, type = "target_fasta", trunk = "sams", ext = ".sam"),
		input_vcf = target_fasta_mh_vcf_from_marker_set_genome_microhap_vcf,
		ui = "{run_dir}/{species_dir}/microhaplot/ui.R",
		shiny = "{run_dir}/{species_dir}/microhaplot/server.R" 
	log:
		"{run_dir}/{species_dir}/logs/microhap_extract_target_fasta/{marker_set}--{target_fasta}--{microhap_variants}.log"
	envmodules:
		"R/4.0.3"  # this is for SEDNA. I can't get a conda-installed R to work
	threads: 10
	output:
		rds="{run_dir}/{species_dir}/microhaplot/{marker_set}--target_fastas--{target_fasta}--{microhap_variants}.rds"
	script:
		"../script/extract_microhaps.R"
	#shell:
	#	"for i in {input.sams}; do echo $i; done  > {output.rds}"


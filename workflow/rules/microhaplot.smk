
# establish the microhaplot directory with the Shiny components in it.
# (This way we can do this operation only once)
rule microhap_shiny_stuff:
	params:
		rd = "{run_dir}"
	output:
		"{run_dir}/microhaplot/ui.R",
		"{run_dir}/microhaplot/server.R"
	log:
		"{run_dir}/logs/microhap_shiny_stuff/microhap-install-etc.log"
#	conda:
#		"../envs/microhap.yaml"
	envmodules:
		"R/4.0.3"  # this is for SEDNA. I can't get a conda-installed R to work
	script:
		"../script/create_microhaplot_folder.R"

# for each collection of SAM files, make an rds file witin an
# microhaplot directory like this
# {run_dir}/microhaplot/{marker_set}--{type}--{additional}.rds
# The filename might be, for example:
# LFAR--fullgex_remapped_to_thin--Otsh_v1.0.rds
#
# or
#
# ROSA--target_fastas--rosawr.rds
#
# The input here is all of the sam files needed to create
# the single output. (This is an aggregation step and
# we want all the dependencies to propagate appropriately.)
# rule: micohaplot
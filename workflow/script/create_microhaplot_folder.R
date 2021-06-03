# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# get the run_dir,  and then create a microhaplot directory
# there with the Shiny elements in it.

run_dir = snakemake@params$rd


inst_packs <- rownames(installed.packages())
if(!("dplyr" %in% inst_packs)) {
	message("Installing the remotes package")
	install.packages("dplyr", repos = "http://cran.rstudio.com")
}
if(!("remotes" %in% inst_packs)) {
	message("Installing the remotes package")
	install.packages("remotes", repos = "http://cran.rstudio.com")
}
message("Installing microhaplotextract")
remotes::install_github(
	"eriqande/microhaplot",
	ref = "just-for-extracting",
	build_vignettes = FALSE,
	build_opts = c("--no-resave-data", "--no-manual"),
	upgrade = "never" 
)

microhaplotextract::mvShinyHaplot(run_dir)

# get rid of the microhaplot example files:
file.remove(
	file.path(
		run_dir, 
		c(
			"microhaplot/fish1.rds",
			"microhaplot/fish1_posinfo.rds",
			"microhaplot/fish2.rds",
			"microhaplot/fish2_posinfo.rds"
		)
	)
)



# here is a little piece to load the config into python
# and to get some other snakemake functions needed
import pandas as pd
from snakemake.utils import validate
from snakemake.io import load_configfile
from snakemake.io import expand
from snakemake.io import unpack

config = load_configfile("config/config.yaml")

config["run_dir"] = "data/test-210129"


# here are some lines to definee a spoofed wildcards variable for testing
class Foo(object):
    pass

wildcards = Foo()


# then after sourcing the input function definitions from common.smk we
# can test them here for correctness:
wildcards.genome = "Otsh_v1.0"
genome_url_from_genome(wildcards)

wildcards.marker_set = "LFAR"



samples_file = r"{run_dir}/{sample_file}".format(
    run_dir=config["run_dir"],
    sample_file = config["samples"]
    )
# read and validate
samples = pd.read_csv(samples_file).set_index("sample", drop=False)



# get path to units file within the run_dir
units_file = r"{run_dir}/{units_file}".format(
    run_dir=config["run_dir"],
    units_file = config["units"]
    )
units = pd.read_csv(units_file, dtype={"sample": str, "Markers": str}).set_index(["sample", "Markers"], drop=False).sort_index()

#### Filter units into two versions: genome-focused and target-fasta-focused
# first get all the markers into lists of genome-focused and target-fasta focused ones
MS = list(config["marker_sets"].keys())
tf_markers = [k for k in MS if "target_fasta" in config["marker_sets"][k]]
gf_markers = [k for k in MS if "genome" in config["marker_sets"][k]]

# figure out which rows in units correspond to each of those tf and gf things.
gflist = units["Markers"].isin(gf_markers)
tflist = units["Markers"].isin(tf_markers)

# break units into genome-focused and target-fasta-focused ones
gf_units = units[gflist == True]
tf_units = units[tflist == True]



bam_tree_equivalent_files_from_marker_sets(wildcards, type = "fullgex_remapped", trunk = "sams", ext = ".sam")





region_files_from_marker_set_and_genome(wildcards)

fna_from_genome(wildcards)


wildcards.run_dir = ".test/data"
wildcards.sample = "CH_003"
fq1_from_sample_and_run(wildcards)
fq2_from_sample_and_run(wildcards)



wildcards.marker_set = "ROSA"
wildcards.target_fasta = "rosa_seqs_1"
fna_from_marker_set_and_target_fasta(wildcards)


wildcards.marker_set = "ROSA"
fullg_bam_inputs_for_calling_from_marker_set_and_genome(wildcards)

target_fasta_bam_inputs_for_calling_from_marker_set_and_fasta(wildcards)


# here we read the units file and make sure it is working right
#units  = pd.read_csv("data/201019_M02749_0096_000000000-CV34D/units.csv", dtype={"sample": str, "Markers": str}).set_index(["sample", "Markers"], drop=False).sort_index()

requested_vcfs_from_units_and_config()


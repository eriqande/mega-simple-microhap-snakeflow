

# here is a little piece to load the config into python
# and to get some other snakemake functions needed
import pandas as pd
from snakemake.utils import validate
from snakemake.io import load_configfile
from snakemake.io import expand

config = load_configfile("config/config.yaml")


# here are some lines to definee a spoofed wildcards variable for testing
class Foo(object):
    pass

wildcards = Foo()


# then after sourcing the input function definitions from common.smk we
# can test them here for correctness:
wildcards.genome = "Otsh_v1.0"
genome_url_from_genome(wildcards)

wildcards.marker_set = "WRAP"
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
units  = pd.read_csv("data/201019_M02749_0096_000000000-CV34D/units.csv", dtype={"sample": str, "Markers": str}).set_index(["sample", "Markers"], drop=False).sort_index()

requested_vcfs_from_units_and_config()


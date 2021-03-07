import pandas as pd
from snakemake.utils import validate


#### Config file and sample spreadsheets ####
configfile: "config/config.yaml"

# get path to samples file within the run_dir
samples_file = r"{run_dir}/{sample_file}".format(
    run_dir=config["run_dir"],
    sample_file = config["samples"]
    )
# read and validate
samples = pd.read_csv(samples_file).set_index("Sample_ID", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


# get path to units file within the run_dir
units_file = r"{run_dir}/{units_file}".format(
    run_dir=config["run_dir"],
    units_file = config["units"]
    )
units = pd.read_csv(units_file, dtype={"Sample_ID": str, "Markers": str}).set_index(["Sample_ID", "Markers"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")



#### Functions for turning wildcards into input values ####

def fna_from_genome(wildcards):
    """Get path to genome fasta from a given genome"""
    return r"resources/genomes/{genome}/{genome}.fna".format(genome=wildcards.genome)

def fq1_from_sample_and_run(wildcards):
    """Get path to a sample's read1 fastq file"""
    return r"{run_dir}/raw/{fq}".format(
        run_dir=wildcards.run_dir,
        fq=samples.loc[wildcards.sample, "fq1"]
    )

def fq2_from_sample_and_run(wildcards):
    """Get path to a sample's read1 fastq file"""
    return r"{run_dir}/raw/{fq}".format(
        run_dir=wildcards.run_dir,
        fq=samples.loc[wildcards.sample, "fq2"]
    )


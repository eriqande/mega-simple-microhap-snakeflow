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

#### Filter units into two versions: genome-focused and target-fasta-focused

gflist = units["Markers"].isin(config["genome_focused_marker_sets"]["name"])
tflist = units["Markers"].isin(config["target_fasta_focused_marker_sets"]["name"])

gf_units = units[gflist == True]
tf_units = units[tflist == True]


#### Functions for turning wildcards into input values ####

def region_files_from_marker_set_and_genome(wildcards):
    """Get path to the regions file given the genome and marker set"""
    return config["marker_sets"][wildcards.marker_set]["genome"][wildcards.genome]["regions"]

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


def fna_from_marker_set_and_target_fasta(wildcards):
    """get path to a target fasta"""
    return config["marker_sets"][wildcards.marker_set]["target_fasta"][wildcards.target_fasta]


def fullg_bam_inputs_for_calling_from_marker_set_and_genome(wildcards):
    """get list of input bams for every sample for a given marker_set and genome"""
    # first, get the pandas data frame of samples for the particular marker set
    DF = gf_units[gf_units["Markers"].isin([wildcards.marker_set])]
    # then cycle over those rows and make a list of paths
    ret = list()
    for index, row in DF.iterrows():
        S = row['Sample_ID']
        ret = ret + [r"{R}/bams/fullg-extracted/{M}/{G}/{S}.bam".format(
        R = wildcards.run_dir,
        M = wildcards.marker_set,
        G = wildcards.genome,
        S = S)]
    return ret


def target_fasta_bam_inputs_for_calling_from_marker_set_and_fasta(wildcards):
    """get list of input bams for every sample for a given marker_set and genome"""
    # first, get the pandas data frame of samples for the particular marker set
    DF = tf_units[tf_units["Markers"].isin([wildcards.marker_set])]
    # then cycle over those rows and make a list of paths
    ret = list()
    for index, row in DF.iterrows():
        S = row['Sample_ID']
        ret = ret + [r"{R}/bams/target_fastas/{M}/{T}/{S}.bam".format(
        R = wildcards.run_dir,
        M = wildcards.marker_set,
        T = wildcards.target_fasta,
        S = S)]
    return ret








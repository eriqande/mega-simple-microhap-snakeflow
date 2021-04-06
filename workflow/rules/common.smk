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
samples = pd.read_csv(samples_file).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


# get path to units file within the run_dir
units_file = r"{run_dir}/{units_file}".format(
    run_dir=config["run_dir"],
    units_file = config["units"]
    )
units = pd.read_csv(units_file, dtype={"sample": str, "Markers": str}).set_index(["sample", "Markers"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

#### Filter units into two versions: genome-focused and target-fasta-focused
#  eventually I should be able to pull these different sets from the marker_sets
# tree in config, but until then I use the pre-made lists in the config
gflist = units["Markers"].isin(config["genome_focused_marker_sets"]["name"])
tflist = units["Markers"].isin(config["target_fasta_focused_marker_sets"]["name"])

gf_units = units[gflist == True]
tf_units = units[tflist == True]


#### Functions for turning wildcards into input values ####

def genome_url_from_genome(wildcards):
    """Get the URL for downloading the genome"""
    return config["genomes"][wildcards.genome]["url"]

def region_files_from_marker_set_and_genome(wildcards):
    """Get path to the regions file given the genome and marker set"""
    return config["marker_sets"][wildcards.marker_set]["genome"][wildcards.genome]["regions"]

def fna_from_genome(wildcards):
    """Get path to genome fasta from a given genome"""
    return r"resources/genomes/{genome}/{genome}.fna".format(genome=wildcards.genome)

def fai_from_genome(wildcards):
    """Get path to genome fasta from a given genome"""
    return r"resources/genomes/{genome}/{genome}.fna.fai".format(genome=wildcards.genome)


def fna_bwt_from_genome(wildcards):
    """Get path to genome fasta from a given genome"""
    return r"resources/genomes/{genome}/{genome}.fna.bwt".format(genome=wildcards.genome)

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
    return config["marker_sets"][wildcards.marker_set]["target_fasta"]["fasta"][wildcards.target_fasta]


def fullg_bam_inputs_for_calling_from_marker_set_and_genome(wildcards):
    """get list of input bams for every sample for a given marker_set and genome"""
    # first, get the pandas data frame of samples for the particular marker set
    DF = gf_units[gf_units["Markers"].isin([wildcards.marker_set])]
    # then cycle over those rows and make a list of paths
    ret = list()
    for index, row in DF.iterrows():
        S = row['sample']
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
        S = row['sample']
        ret = ret + [r"{R}/bams/target_fastas/{M}/{T}/{S}.bam".format(
        R = wildcards.run_dir,
        M = wildcards.marker_set,
        T = wildcards.target_fasta,
        S = S)]
    return ret



#### Functions for defining output files from units and config ####


# from unit.csv (now available in gf_units and tf_units) get all the different
# vcf outputs that we should be expecting.  This means cycling over all the
# different genomes and target_fasta seqs in config, as well.
# BIG NOTE: This is explicitly for creating VCFs from the bams.  These VCFs will
# then be available to merge and define new "canonical" variation. But they
# are not directly used to feed into the microhaplot or SNP-yanking sections
# of the workflow.
def requested_vcfs_from_units_and_config():
    """get list of vcf outputs we expect from the units and the config file"""
    # start with the genome-focused ones
    # here are the unique marker sets called for:
    MS = list(set(list(gf_units["Markers"])))
    # now, expand each of those by the genomes they might be associated with
    gf = list()
    for m in MS:
        # list of full genomes they are associated with
        g = [str(k) for k in config["marker_sets"][m]["genome"].keys()]
        gf = gf + expand("{rd}/vcfs/{ms}/fullg/{g}/variants-bcftools.vcf", rd = config["run_dir"], ms = m, g = g)
    # then do the target-fasta focused ones
    MS = list(set(list(tf_units["Markers"])))
    # now, expand each of those by the specific target-fasta seqs they might be associated with
    tf = list()
    for m in MS:
        # list of target-fasta fastas they are associated with
        t = [str(k) for k in config["marker_sets"][m]["target_fasta"]["fasta"].keys()]
        tf = tf + expand("{rd}/vcfs/{ms}/target_fasta/{t}/variants-bcftools.vcf", rd = config["run_dir"], ms = m, t = t)
    return gf + tf




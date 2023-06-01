import pandas as pd
from snakemake.utils import validate
from glob import glob


#### Global wildcard constraints ####
#wildcard_constraints:
#    run_dir = "^(?!.*MULTI_RUN_RESULTS.*)"


#### Config file and sample spreadsheets ####

# Chinook is the default, but the configfile should
# really be set on the command line with --configfile config/XXXXX/config.yaml
#configfile: "config/Chinook/config.yaml"  

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


#### Functions for turning wildcards into input values ####

def genome_url_from_genome(wildcards):
    """Get the URL for downloading the genome"""
    return config["genomes"][wildcards.genome]["url"]

def region_files_from_marker_set_and_genome(wildcards):
    """Get path to the regions file given the genome and marker set"""
    return config["marker_sets"][wildcards.marker_set]["genome"][wildcards.genome]["regions"]

def fq1_or_trim1_from_sample_and_run(wildcards):
    """Get path to a sample's read1 fastq file"""
    if not config["use_trimmomatic"]:
        return r"{run_dir}/raw/{fq}".format(
            run_dir=wildcards.run_dir,
            fq=samples.loc[wildcards.sample, "fq1"]
        )
    else:
        return r"{run_dir}/{species_dir}/trimmomatic/{sm}.1.fastq.gz".format(
            run_dir=wildcards.run_dir,
            species_dir=config["species"],
            sm=wildcards.sample)

def fq2_or_trim2_from_sample_and_run(wildcards):
    """Get path to a sample's read1 fastq file"""
    if not config["use_trimmomatic"]:
        return r"{run_dir}/raw/{fq}".format(
            run_dir=wildcards.run_dir,
            fq=samples.loc[wildcards.sample, "fq2"]
        )
    else:
        return r"{run_dir}/{species_dir}/trimmomatic/{sm}.2.fastq.gz".format(
            run_dir=wildcards.run_dir,
            species_dir=config["species"],
            sm=wildcards.sample)


def fq1_from_sample_and_run(wildcards):
    """Get path to a sample's read1 fastq file"""
    return r"{run_dir}/raw/{fq}".format(
        run_dir=wildcards.run_dir,
        fq=samples.loc[wildcards.sample, "fq1"]
    )

def fq2_from_sample_and_run(wildcards):
    """Get path to a sample's read2 fastq file"""
    return r"{run_dir}/raw/{fq}".format(
        run_dir=wildcards.run_dir,
        fq=samples.loc[wildcards.sample, "fq2"]
    )


def rg_from_sample(wildcards):
    """get the read group to put in a sample's BAM file"""
    # key thing here is the set ID to {sample} and SM to NMFS_DNA_ID
    return r"@RG\tID:{sample}\tSM:{nmfs}\tPL:ILLUMINA".format(
        sample = wildcards.sample,
        nmfs = samples.loc[wildcards.sample, "NMFS_DNA_ID"]
    )

def fna_from_marker_set_and_target_fasta(wildcards):
    """get path to a target fasta stored in the Config."""
    return config["marker_sets"][wildcards.marker_set]["target_fasta"][wildcards.target_fasta]["fasta"]


# here is a more general version  that I can call
# to also do idxstats, sams, and different extenstions with
# it by calling it in a lambda function.
# type is either fullg or target_fasta or fullgex_remapped
# trunk is whether it is bam, sam, or idxstats, etc.
# ext is the extensions the file should have.  IT MUST INCLUDE THE PERIOD (i.e. ".bam", not "bam") 
def bam_tree_equivalent_files_from_marker_sets(wildcards, type, trunk, ext):
    if(type == "fullg"):
        # first, get the pandas data frame of samples for the particular marker set
        DF = gf_units[gf_units["Markers"].isin([wildcards.marker_set])]
        # then cycle over those rows and make a list of paths
        ret = list()
        for index, row in DF.iterrows():
            S = row['sample']
            ret = ret + [r"{R}/{species_dir}/{TRUNK}/fullg-extracted/{M}/{G}/{S}{EXT}".format(
            R = wildcards.run_dir,
            species_dir = wildcards.species_dir,
            TRUNK = trunk,
            M = wildcards.marker_set,
            G = wildcards.genome,
            S = S,
            EXT = ext)]
        return ret
    else:
        if(type == "target_fasta"):
            # first, get the pandas data frame of samples for the particular marker set
            DF = tf_units[tf_units["Markers"].isin([wildcards.marker_set])]
            # then cycle over those rows and make a list of paths
            ret = list()
            for index, row in DF.iterrows():
                S = row['sample']
                ret = ret + [r"{R}/{species_dir}/{TRUNK}/target_fastas/{M}/{T}/{S}{EXT}".format(
                R = wildcards.run_dir,
                species_dir = wildcards.species_dir,
                TRUNK = trunk,
                M = wildcards.marker_set,
                T = wildcards.target_fasta,
                S = S,
                EXT = ext)]
            return ret
        else:
            if(type == "fullgex_remapped"):
                # first, get the pandas data frame of samples for the particular marker set
                DF = gf_units[gf_units["Markers"].isin([wildcards.marker_set])]
                # then cycle over those rows and make a list of paths
                ret = list()
                for index, row in DF.iterrows():
                    S = row['sample']
                    ret = ret + [r"{R}/{species_dir}/{TRUNK}/fullgex_remapped_to_thinned/{M}/{G}/{S}{EXT}".format(
                    R = wildcards.run_dir,
                    species_dir = wildcards.species_dir,
                    TRUNK = trunk,
                    M = wildcards.marker_set,
                    G = wildcards.genome,
                    S = S,
                    EXT = ext)]
                return ret
            else:
                raise ValueError("type must be fullg or target_fasta or fullgex_remapped")



# get the canonical variants VCF for microhaplot for fullgex_remapped
def fullgex_remapped_mh_vcf_from_marker_set_genome_microhap_vcf(wildcards):
    return config["marker_sets"][wildcards.marker_set]["genome"][wildcards.genome]["microhap_variants"][wildcards.microhap_variants]


# get the canonical variants VCF for microhaplot for target fastas
def target_fasta_mh_vcf_from_marker_set_genome_microhap_vcf(wildcards):
    return config["marker_sets"][wildcards.marker_set]["target_fasta"][wildcards.target_fasta]["microhap_variants"][wildcards.microhap_variants]



# this is for returning the BAMs of a particular species, marker set, 
# mapping type, genome, etc. from multiple completed runs
def existing_bams_for_multirun_fullgex_remapped(wildcards):
    # first get name of files 
    dirs_file = "MULTI_RUN_RESULTS/{md_dir}/dirs.txt".format(md_dir=wildcards.multi_run_dir)
    with open(dirs_file) as f:
        dirs_list = f.readlines()
    dirs_list = [x.strip() for x in dirs_list if x.strip()] # strip whitespace and empty lines
    globs_list = ["{d}/{s}/bams/fullgex_remapped_to_thinned/{m}/{g}/*.bam".format(d=dir, s=wildcards.species_dir, m=wildcards.marker_set, g=wildcards.genome) for dir in dirs_list]
    bam_list = list()
    for G in globs_list:
        bam_list = bam_list + glob(G)
    return bam_list



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
        gf = gf + expand("{rd}/{sd}/vcfs/{ms}/fullg/{g}/variants-bcftools.vcf", rd = config["run_dir"], sd = config["species"], ms = m, g = g)
    # then do the target-fasta focused ones
    MS = list(set(list(tf_units["Markers"])))
    # now, expand each of those by the specific target-fasta seqs they might be associated with
    tf = list()
    for m in MS:
        # list of target-fasta fastas they are associated with
        t = [str(k) for k in config["marker_sets"][m]["target_fasta"].keys()]
        tf = tf + expand("{rd}/{sd}/vcfs/{ms}/target_fasta/{t}/variants-bcftools.vcf", rd = config["run_dir"], sd = config["species"], ms = m, t = t)
    # then, also do the genome-focused ones that have been extracted and remapped to the
    # thinned genomes.    
    MS = list(set(list(gf_units["Markers"])))
    # now, expand each of those by the genomes they might be associated with
    gfex = list()
    for m in MS:
        # list of full genomes they are associated with
        g = [str(k) for k in config["marker_sets"][m]["genome"].keys()]
        gfex = gfex + expand("{rd}/{sd}/vcfs/{ms}/fullgex_remapped/{g}/variants-bcftools.vcf", rd = config["run_dir"], sd = config["species"], ms = m, g = g)

    return gf + gfex + tf




# this is incomplete.  Currently it just gets the fullgex_and_mapped_to_extracted ones
# because the others need some revamping in the config tree.
def requested_microhap_outputs_from_units_and_config():
    # do the genome-focused ones that have been extracted and remapped to the
    # thinned genomes.    
    MS = list(set(list(gf_units["Markers"])))
    # now, expand each of those by the genomes they might be associated with
    gfex = list()
    for m in MS:
        for g in [str(k) for k in config["marker_sets"][m]["genome"].keys()]:
            for v in  [str(k) for k in config["marker_sets"][m]["genome"][g]["microhap_variants"].keys()]:
                gfex = gfex + ["{run_dir}/{species_dir}/microhaplot/{marker_set}--fullgex_remapped_to_thinned--{genome}--{microhap_variants}.rds".format(
                    run_dir = config["run_dir"],
                    species_dir = config["species"],
                    marker_set = m,
                    genome = g,
                    microhap_variants = v
                )]
    # now we will do the target fasta ones
    MS = list(set(list(tf_units["Markers"])))
    tf = list()
    for m in MS:
        for t in [str(k) for k in config["marker_sets"][m]["target_fasta"].keys()]:
            for v in [str(k) for k in config["marker_sets"][m]["target_fasta"][t]["microhap_variants"].keys()]:
                tf = tf + ["{run_dir}/{species_dir}/microhaplot/{marker_set}--target_fastas--{targ_fast}--{microhap_variants}.rds".format(
                    run_dir = config["run_dir"],
                    species_dir = config["species"],
                    marker_set = m,
                    targ_fast = t,
                    microhap_variants = v
                )]
    return gfex + tf



# generalized version: trunk should be like, vcf of idxstats, etc.
# filename should be the filename..

# from unit.csv (now available in gf_units and tf_units) get all the different
# vcf outputs that we should be expecting.  This means cycling over all the
# different genomes and target_fasta seqs in config, as well.
# BIG NOTE: This is explicitly for creating VCFs from the bams.  These VCFs will
# then be available to merge and define new "canonical" variation. But they
# are not directly used to feed into the microhaplot or SNP-yanking sections
# of the workflow.
def requested_outputs_from_units_and_config(trunk, filename):
    """get list of vcf outputs we expect from the units and the config file"""
    # start with the genome-focused ones
    # here are the unique marker sets called for:
    MS = list(set(list(gf_units["Markers"])))
    # now, expand each of those by the genomes they might be associated with
    gf = list()
    for m in MS:
        # list of full genomes they are associated with
        g = [str(k) for k in config["marker_sets"][m]["genome"].keys()]
        gf = gf + expand("{rd}/{sd}/{tr}/fullg-extracted/{ms}/{g}/{fn}", rd = config["run_dir"], sd = config["species"], tr = trunk, ms = m, g = g, fn = filename)
    # then do the target-fasta focused ones
    MS = list(set(list(tf_units["Markers"])))
    # now, expand each of those by the specific target-fasta seqs they might be associated with
    tf = list()
    for m in MS:
        # list of target-fasta fastas they are associated with
        t = [str(k) for k in config["marker_sets"][m]["target_fasta"].keys()]
        tf = tf + expand("{rd}/{sd}/{tr}/target_fastas/{ms}/{t}/{fn}", rd = config["run_dir"], sd = config["species"], tr = trunk, ms = m, t = t, fn = filename)
    # then, also do the genome-focused ones that have been extracted and remapped to the
    # thinned genomes.    
    MS = list(set(list(gf_units["Markers"])))
    # now, expand each of those by the genomes they might be associated with
    gfex = list()
    for m in MS:
        # list of full genomes they are associated with
        g = [str(k) for k in config["marker_sets"][m]["genome"].keys()]
        gfex = gfex + expand("{rd}/{sd}/{tr}/fullgex_remapped_to_thinned/{ms}/{g}/{fn}", rd = config["run_dir"], sd = config["species"], tr = trunk, ms = m, g = g, fn = filename)

    return gf + gfex + tf



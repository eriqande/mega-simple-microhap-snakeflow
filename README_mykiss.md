README\_mykiss
================

-   [Running A Workflow with Trimmomatic (for Jeff Rodzen) for
    Mykiss](#running-a-workflow-with-trimmomatic-for-jeff-rodzen-for-mykiss)
-   [Preparing samples.csv and
    units.csv](#preparing-samplescsv-and-unitscsv)
-   [The config tree](#the-config-tree)
-   [Getting the `starter_variation`](#getting-the-starter_variation)
-   [What things looked like](#what-things-looked-like)

## Running A Workflow with Trimmomatic (for Jeff Rodzen) for Mykiss

The whole pipeline can be run like this:

``` sh
snakemake --config run_dir=data/Rodzen-data/0045_FRH_steelhead use_trimmomatic=True --configfile config/Mykiss/config.yaml --use-conda --use-envmodules --cores 20
```

Note that if you give the `use_trimmomatic=True` config on the command
line, you have to use Python formatted logicals—i.e., either `True` or
`False` written exactly as shown. (`true` or `TRUE` will not work for
you!)

Also note that `--use-envmodules` is specific to my cluster. It allows
it to use R, so it probably should be dropped if running on one’s own
Linux box with R available by default.

This is a start for Mykiss, but it still needs a lot of cleaning,
filtering, and curating, because some of those amplicons may be
amplifying multiple regions.

## Preparing samples.csv and units.csv

At NMFS we do things a little non-standardly in that we store the
NMFS\_DNA\_ID as a piece of the Sample\_Plate column of the
`SampleSheet.csv`. Jeff has his sample ID’s in the Sample\_ID column of
`SampleSheet.csv`, so an additional parameter is required when preparing
the `samples.csv` and `units.csv` files with R:

``` r
# do this within R from the top level of the
# mega-simple-microhap-snakeflow directory
source("preprocess/sample-sheet-processing-functions.R")

create_samples_and_units(
    # point it to the SampleSheet in whatever data directory you are working with
    "data/Rodzen-data/0045_FRH_steelhead/SampleSheet.csv", 
    NMFS_DNA_ID_from_Sample_ID = TRUE  # This is the extra required option
)
```

For this to work, Eric renamed Jeff’s original
`SampleSheet_FRH_Steelhead_GTseq_110720.csv` to `SampleSheet.csv` and
then added `OMY_RODZEN` to every sample row in the Description column,
and added `Test_for_mykiss` to every sample row in the Sample\_Project
column. The first column added tells it what marker set these were typed
at and the second is a required column.

The result looks like this:

``` csv
[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Investigator Name,Jeff,,,,,,,,
Experiment Name,FRH_Steelhead_GTseqtest,,,,,,,,
Date,11/7/20,,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,FRH_Steelhead_GTseqtest,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
76,,,,,,,,,
76,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,,,,,,,,,
AdapterRead2,,,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
M190760FRH,RT3497,FRH20190407-01,C12,NMFS_GTseq_i7_i001,CGTGAT,i5-CGTCTA-36,CGTCTA,Test_for_mykiss,OMY_RODZEN
M190762FRH,RT3498,FRH20190407-01,D12,NMFS_GTseq_i7_i001,CGTGAT,i5-TGGGGA-48,TGGGGA,Test_for_mykiss,OMY_RODZEN
M190819FRH,RT3500,FRH20190407-01,F12,NMFS_GTseq_i7_i001,CGTGAT,i5-TTCTAG-72,TTCTAG,Test_for_mykiss,OMY_RODZEN

...and so forth...
```

When running the workflow, you should start with a data directory that
includes `SampleSheet.csv` and the directory `raw` that includes all the
fastq.gz files. Then:

1.  create `samples.csv` and `units.csv` as shown above.
2.  Run the workflow with the snakemake command shown above.

## The config tree

The config file `config/Mykiss/config.yaml` at this juncture, looks
something like: Anyway, we have some starter VCFs for variation now, and
they are listed in the config:

``` yaml
marker_sets:
  OMY_RODZEN:
    genome:
      Omyk_1.0:
        regions: config/Mykiss/regions/OMY_RODZEN-Omyk_1.0.txt
        microhap_variants:
          starter_variation: config/Mykiss/canonical_variation/myk-rodzen-fullgex-484-sites-from-0045_FRH.vcf
    target_fasta:
      omy_rodzen_tf:
        fasta: config/Mykiss/target_fastas/rodzen_mykiss_GTseq.fasta
        microhap_variants:
          starter_variation: config/Mykiss/canonical_variation/myk-rodzen-target-fasta-545-sites-from-0045_FRH.vcf
```

The `starter_variation` VCF sets are just produced from some of the
variation that I saw in the initial runs, as described below. Note that
this will likely need further filtering and cleaning, and certainly more
variants will need to be added to it as more populations are typed at
these markers. (See, for example, [these
steps](https://github.com/eriqande/mega-simple-microhap-snakeflow#multi-run-variant-calling))

## Getting the `starter_variation`

I first ran the workflow with empty VCF files for the variation:

``` sh
snakemake --config run_dir=data/Rodzen-data/0045_FRH_steelhead use_trimmomatic=True --configfile config/Mykiss/config.yaml --use-conda --use-envmodules --cores 20 --keep-going
```

I added the `--keep-going` up there because I knew that it was going to
fail at the microhap step, because I don’t have any variants for these
guys. But after running through it once, I have the VCF files produced
by bcftools, so we can use those for preliminary variants to pass to
microhaplot.

After that run completed (and I fixed a few things with the pipeline) I
used the resulting VCF file as a starting VCF for micohaplot, like this:

``` sh
# filter out any sites that are missing in more than
# 50% of the indivs, and then just keep one indiv

#: in: /home/eanderson/Documents/git-repos/mega-simple-microhap-snakeflow/data/Rodzen-data/0045_FRH_steelhead/Mykiss/vcfs/OMY_RODZEN/fullgex_remapped/Omyk_1.0

bcftools view -i 'F_MISSING < 0.5' variants-bcftools.vcf | bcftools view -s M160023FRH > ../../../../../../../../config/Mykiss/canonical_variation/myk-rodzen-fullgex-484-sites-from-0045_FRH.vcf

# for the target fasta ones, we also had to drop one site.
bcftools view -t ^Omy_RAD7016-31:208 -i 'F_MISSING < 0.5' variants-bcftools.vcf |  bcftools view -s M160023FRH > ../../../../../../../../config/Mykiss/canonical_variation/myk-rodzen-target-fasta-545-sites-from-0045_FRH.vcf
```

## What things looked like

A cursory look at things in microhaplot suggested that there were at
least a few SNPs that are fully heterozygous, suggesting that perhaps
they are amplifying multiple regions.

Jeff will need to go through all his markers and figure out which ones
to use.

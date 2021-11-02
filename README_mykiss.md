README\_mykiss
================

Eric has added some marker sets for Mykiss into the flow in order to get
one of Jeff Rodzen’s marker sets in as as example.

To run this, Eric took the `SampleSheet_FRH_Steelhead_GTseq_110720.csv`
and added `OMY_RODZEN` to every sample row in the Description column,
and added `Test_for_mykiss` to every sample row in the Sample\_Project
column. The first column added tells it what marker set these were typed
at and the second is a required column. I then renamed the result to
SampleSheet.csv\`.

Then, I copied the `0045_FRH_steelhead` to
`/home/eanderson/Documents/git-repos/mega-simple-microhap-snakeflow/data/Rodzen-data/0045_FRH_steelhead`
on our cluster and removed only `SampleSheet.csv` and the fastq.gz files
in `raw`. So, a listing looks like this:

``` sh
(snakemake) [node08: mega-simple-microhap-snakeflow]--% pwd
/home/eanderson/Documents/git-repos/mega-simple-microhap-snakeflow
(snakemake) [node08: mega-simple-microhap-snakeflow]--% ls data/Rodzen-data/0045_FRH_steelhead/
raw  SampleSheet.csv
```

Then, we have to create the samples and units files:

``` r
# working within R
> getwd()
[1] "/home/eanderson/Documents/git-repos/mega-simple-microhap-snakeflow"

> source("preprocess/sample-sheet-processing-functions.R")

> create_samples_and_units("data/Rodzen-data/0045_FRH_steelhead/SampleSheet.csv", NMFS_DNA_ID_from_Sample_ID = TRUE)
```

Note that we add the option `NMFS_DNA_ID_from_Sample_ID = TRUE` because
that is where Jeff’s main IDs are. (Rather than our rather nonstandard
incorporation of those ID’s into the Sample\_Plate field).

And when that is done, we run it!

``` sh
snakemake --config run_dir=data/Rodzen-data/0045_FRH_steelhead use_trimmomatic=True --configfile config/Mykiss/config.yaml --use-conda --use-envmodules --cores 20 --keep-going
```

Note that if you give the `use_trimmomatic=True` config on the command
line, you have to use Python formatted logicals—i.e., either `True` or
`False` written exactly as shown. (`true` or `TRUE` will not work for
you!)

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

So, what I saw from that was a few sites that are fully heterozygous.
Jeff needs to go through all his markers and figure out which ones to
use.

Anyway, we have some starter VCFs for variation now, and they are listed
in the config:

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

And now, the whole pipeline can be run like this:

``` sh
snakemake --config run_dir=data/Rodzen-data/0045_FRH_steelhead use_trimmomatic=True --configfile config/Mykiss/config.yaml --use-conda --use-envmodules --cores 20
```

Note that `--use-envmodules` is specific to my cluster. It allows it to
use R.

Anyway, that is a start for Mykiss, but it needs a lot of cleaning,
filtering, and curating, because some of those amplicons are probably
amplifying multiple regions.

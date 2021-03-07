mega-simple-microhap-snakeflow
================

This repository holds the Snakefile and associated files for Eric’s
first attempt (still under construction\!) at making a SnakeMake-based
workflow for our microhaplotypes.

# The way forward

## samples.csv and units.csv

After a lot of learning and then figuring out what we need for our
microhap workflow, this is how we are going to do it:

1.  Each run will come off the machine, and it will be associated with a
    `samples.csv` file that gives the run-names (i.e. like
    CH\_117\_LFAR\_001) of the individual samples. These are apparently
    based on the individually barcoded units in the run, so each
    run-name gets broken into its own pair of fastq files. There might
    be multiple species on a run.
2.  The first step, which is not part of the snakemake workflow is to
    use an R script to break the run into different folders, one for
    each species. Each folder has a `samples.csv` file and a `raw`
    directory that houses the fastq files. This step also creates a
    `units.csv` file that includes one line for each “unit” which in
    this case is the combination of a run-name *and* a marker set (which
    would be based on the PCR primers, like LFAR, WRAP, ROSA,
    TRANSITION, etc.)

A little explanation will be helpful here. In `.test/data/samples.csv` I
have a simple example:

``` r
library(tidyverse)

samples <- read_csv(".test/data/samples.csv")

samples
```

    ## # A tibble: 8 x 6
    ##   Sample_ID NMFS_DNA_ID Plate Marker_Sets     fq1               fq2             
    ##   <chr>     <chr>       <chr> <chr>           <chr>             <chr>           
    ## 1 CH_001    T0123       W56   LFAR            .test/data/raw/C… .test/data/raw/…
    ## 2 CH_002    T0124       W56   ROSA            .test/data/raw/C… .test/data/raw/…
    ## 3 CH_003    T0125       W56   WRAP            .test/data/raw/C… .test/data/raw/…
    ## 4 CH_004    T0126       W56   TRANSITION      .test/data/raw/C… .test/data/raw/…
    ## 5 CH_005    T0127       W56   LFAR,WRAP       .test/data/raw/C… .test/data/raw/…
    ## 6 CH_006    T0128       W56   TRANSITION,ROSA .test/data/raw/C… .test/data/raw/…
    ## 7 CH_007    T0129       W56   ROSA,TRANSITIO… .test/data/raw/C… .test/data/raw/…
    ## 8 CH_008    T0130       W56   WRAP,ROSA,LFAR  .test/data/raw/C… .test/data/raw/…

From this, the units.tsv file would have been made, something like this:

``` r
units <- samples %>%
  mutate(Markers = str_split(Marker_Sets, "\\s*,\\s*")) %>%
  unnest(cols = Markers) %>%
  select(-Marker_Sets) %>%
  select(Sample_ID, Markers, everything())
```

Those have been stored in `.test/data/units.csv` and look like this:

``` r
read_csv(".test/data/units.csv")
```

    ## # A tibble: 14 x 6
    ##    Sample_ID Markers   NMFS_DNA_ID Plate fq1                 fq2                
    ##    <chr>     <chr>     <chr>       <chr> <chr>               <chr>              
    ##  1 CH_001    LFAR      T0123       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ##  2 CH_002    ROSA      T0124       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ##  3 CH_003    WRAP      T0125       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ##  4 CH_004    TRANSITI… T0126       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ##  5 CH_005    LFAR      T0127       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ##  6 CH_005    WRAP      T0127       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ##  7 CH_006    TRANSITI… T0128       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ##  8 CH_006    ROSA      T0128       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ##  9 CH_007    ROSA      T0129       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ## 10 CH_007    TRANSITI… T0129       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ## 11 CH_007    WRAP      T0129       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ## 12 CH_008    WRAP      T0130       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ## 13 CH_008    ROSA      T0130       W56   .test/data/raw/CH-… .test/data/raw/CH-…
    ## 14 CH_008    LFAR      T0130       W56   .test/data/raw/CH-… .test/data/raw/CH-…

## The config file

This file is in `config/newconfig.yaml`. So far it looks like this:

``` yaml
species: Chinook

run: gtseq01

# the following are names of files assumed to be
# within the run directory
samples: .test/data/samples.csv
units: .test/data/units.csv

# set genomes up to be downloaded from url.  They will
# be placed in resources/genomes/name_of_genome/name_of_genome.fna
# where name_of_genome is like Otsh_v1.0
genomes:
  Otsh_v1.0:
    url: put_here_the_download_url


genome_focused_marker_sets:
  name: ["LFAR", "WRAP", "ROSA"]

target_fasta_focused_marker_sets:
  name: TRANSITION

marker_sets:
  LFAR:
    genome:
      Otsh_v1.0:
        regions: config/regions/LFAR-Otsh_v1.0.txt
  WRAP:
    genome:
      Otsh_v1.0:
        regions: config/regions/WRAP-Otsh_v1.0.txt
  ROSA:
    genome:
      Otsh_v1.0:
        regions: config/regions/ROSA-Otsh_v1.0.txt
  TRANSITION:
    target_fasta: config/target_fastas/transition-panel.fna

```

# Notes on Development

For someone like me who is not super familiar with python, it is nice to
be able to test various constructions in the python interpreter. In
order to have access to variables like the `config` dictionary when
doing that, we can read it in with functions from snakemake’s io module
like this:

``` python
from snakemake.io import load_configfile
config = load_configfile("config/config.yaml")

# then test to see if we have things where we expect them:
config["run"]
config["marker_sets"]["LFAR"]["genome"].keys()
config["marker_sets"]["LFAR"]["genome"]["Otsh_v1.0"]["regions"]

# or this, which returns a list of the regions files for LFAR over
# all the genomes (although there is only one)
[config["marker_sets"]["LFAR"]["genome"][x]["regions"] for x in config["marker_sets"]["LFAR"]["genome"].keys()]

# or get a list of all the marker sets that are genome-focused
config["genome_focused_marker_sets"]["names"]
```

# Genomes and target\_fastas

The target fastas are going to be stored in `config/target_fastas`. The
genomes will have to be downloaded and stored in
`resources/genomes/Otsh_v1.0`, for example.

For testing, I just touch such a genome file and a target\_fastas file:

``` sh
mkdir -p resources/genomes/Otsh_v1.0/
touch resources/genomes/Otsh_v1.0/Otsh_v1.0.fna

mkdir -p config/target_fastas
touch config/target_fastas/transition-panel.fna
```

# Old stuff

To use it, place the flashed fastqs in the `data/samples` directory.
(Note\! We will put the flashing step in here too, soon\!). Then list
the names of the samples and the paths to those fastqs in the
`config.yaml` file like you see for the simple 5 fastq case:

``` yaml
samples:
  CHLFAR1: data/fastq/CHLFAR1.extendedFrags.fastq.gz
  CHLFAR2: data/fastq/CHLFAR2.extendedFrags.fastq.gz
  CHLFAR3: data/fastq/CHLFAR3.extendedFrags.fastq.gz
  CHLFAR4: data/fastq/CHLFAR4.extendedFrags.fastq.gz
  CHLFAR5: data/fastq/CHLFAR5.extendedFrags.fastq.gz
```

And give the path to the genome (assumed to be indexed by bwa-mem
already…, but we could incorporate that too…) there as well, like:

``` yaml
genome: genome/Otsh_v1.0_genomic.fna
```

And finally, list the genomic coordinates of the target microhaplotypes
in the full genome:

``` yaml
regions:
  "NC_037130.1:1847885-1848185":
  "NC_037130.1:1062935-1063235":
  "NC_037130.1:786919-787219":
  "NC_037130.1:929833-930133":
  "NC_037130.1:536847-537147":
  "NC_037130.1:1845977-1846277":
  "NC_037130.1:828619-828919":
  "NC_037130.1:864908-865208":
```

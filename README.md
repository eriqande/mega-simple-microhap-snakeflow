mega-simple-microhap-snakeflow
================

  - [The way forward](#the-way-forward)
      - [samples.csv and units.csv](#samples.csv-and-units.csv)
      - [The config file](#the-config-file)
  - [Notes on Development](#notes-on-development)
      - [Testing python snippets](#testing-python-snippets)
          - [Interpreter access to `config`
            variable](#interpreter-access-to-config-variable)
          - [Interactive python-session testing of input
            functions](#interactive-python-session-testing-of-input-functions)
      - [Getting filepaths and wildcards
        right](#getting-filepaths-and-wildcards-right)
          - [Put fast commands in to check the wildcard
            framework](#put-fast-commands-in-to-check-the-wildcard-framework)
          - [Definitely get the `tree`
            command](#definitely-get-the-tree-command)
          - [Remember to use Unix file iteration/alternation to specify
            requested output
            files](#remember-to-use-unix-file-iterationalternation-to-specify-requested-output-files)
  - [Genomes and target\_fastas](#genomes-and-target_fastas)
  - [Old stuff](#old-stuff)

This repository holds the Snakefile and associated files for Eric’s
first attempt (still under construction\!) at making a SnakeMake-based
workflow for our microhaplotypes.

# The way forward

OK, after a lot of thinking and experimenting, this is how I want to do
it:

The input data will all come in a single directory with the following
items in it:

  - `samples.csv`: the list of samples
  - `units.csv`: the units file.
  - `raw/`: a directory full of the paired end fastq file.

In the config file (or more likely on the command line) we specify the
path (relative or absolute) of this directory, calling it the `run_dir`.

Actually, we might allow that the `run_dir` could be a list, if we
wanted it to be. In which case the workflow would operate over multiple
directories. Actually, that would be a huge PITA, so I will not do that…

But, I think that I will specify the run\_dir as a wildcard in the
workflow, because otherwise I have to repeatedly use the format method
on a bunch of paths…

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
    ##   Sample_ID NMFS_DNA_ID Plate Marker_Sets        fq1             fq2            
    ##   <chr>     <chr>       <chr> <chr>              <chr>           <chr>          
    ## 1 CH_001    T0123       W56   LFAR               CH-001-aaqr-R1… CH-001-aaqr-R2…
    ## 2 CH_002    T0124       W56   ROSA               CH-002-aaqr-R1… CH-002-aaqr-R2…
    ## 3 CH_003    T0125       W56   WRAP               CH-003-aaqr-R1… CH-003-aaqr-R2…
    ## 4 CH_004    T0126       W56   TRANSITION         CH-004-aaqr-R1… CH-004-aaqr-R2…
    ## 5 CH_005    T0127       W56   LFAR,WRAP          CH-005-aaqr-R1… CH-005-aaqr-R2…
    ## 6 CH_006    T0128       W56   TRANSITION,ROSA    CH-006-aaqr-R1… CH-006-aaqr-R2…
    ## 7 CH_007    T0129       W56   ROSA,TRANSITION,W… CH-007-aaqr-R1… CH-007-aaqr-R2…
    ## 8 CH_008    T0130       W56   WRAP,ROSA,LFAR     CH-008-aaqr-R1… CH-008-aaqr-R2…

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

run_dir: .test/data

# the following are names of files that assumed to be
# within the run_dir
samples: samples.csv
units: units.csv

# set genomes up to be downloaded from url.  They will
# be placed in resources/genomes/name_of_genome/name_of_genome.fna
# where name_of_genome is like Otsh_v1.0
genomes:
  Otsh_v1.0:
    url: put_here_the_download_url
  nookie2:
    url: this_would_need_a_download_url_too


# The name element must be specified as a list here!
genome_focused_marker_sets:
  name: ["LFAR", "WRAP", "ROSA"]

# The name element must be specified as a list here!
target_fasta_focused_marker_sets:
  name: ["TRANSITION", "ROSA"]


# here we give details about our marker sets.  Note that a marker
# set can be both genome-focused and target_fasta-focused, as
# I have done with ROSA here, just for an example.

# Note also that we have allowed for two different target_fastas
# for ROSA which will give us two different possible ROSA sets, possibly.
# It seemed like that would be good in order to test different
# target_fastas side-by-side, etc.

# Also, each target_fasta is assumed to belong to only one marker_set,
# but a particular full genome can belong to multiple marker sets.

# And we can map a marker set to more than one full genome...as the
# LFAR example with a genome called "nookie2" shows.
marker_sets:
  LFAR:
    genome:
      Otsh_v1.0:
        regions: config/regions/LFAR-Otsh_v1.0.txt
      nookie2:
        regions: config/regions/LFAR-nookie2.txt
  WRAP:
    genome:
      Otsh_v1.0:
        regions: config/regions/WRAP-Otsh_v1.0.txt
  ROSA:
    genome:
      Otsh_v1.0:
        regions: config/regions/ROSA-Otsh_v1.0.txt
    target_fasta:
      rosa_seqs_1: config/target_fastas/rosa-seqs-1-example.fasta
      rosa_seqs_2: config/target_fastas/rosa-seqs-2-example.fasta
  TRANSITION:
    target_fasta:
      transition-panel: config/target_fastas/transition-panel.fna
```

# Notes on Development

## Testing python snippets

For someone like me who is not super familiar with python, it is nice to
be able to test various constructions in the python interpreter. Here
are some fun tricks to allows that to happen.

### Interpreter access to `config` variable

In order to have access to variables like the `config` dictionary when
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

### Interactive python-session testing of input functions

Often we need to write functions that return possible file names and
paths to the `input` field of a rule. It is nice to be able to test
those just in an interactive python session. We can do that by making
our own `wildcards` object, and then assigning attributes to it and
passing it to our functions.

That goes like this:

``` python

# create a generic object class (called Foo here)
class Foo(object):
    pass

# make a variable, wildcards, which is an object of that class
wildcards = Foo()

# now, if you want to test specific values you can do like:
wildcards.marker_set = "LFAR"
wildcards.genome = "Otsh_v1.0"

# then, if you have a function that uses the marker_set and the genome
# attributes of the wildcards, like:

def region_files_from_marker_set_and_genome(wildcards):
    """Get path to the regions file given the genome and marker set"""
    return config["marker_sets"][wildcards.marker_set]["genome"][wildcards.genome]["regions"]

# then you can test it (assuming you have defined the config variable into your
# interactive session of python using the trick from above), like this:

region_files_from_marker_set_and_genome(wildcards)

# and the result is:
'config/regions/LFAR-Otsh_v1.0.txt'

# which is just what we want it to be.
```

## Getting filepaths and wildcards right

Just understanding the wildcarding framework of snakemake can take a
while to get your head around. This microhap workflow is pretty beastly
because there are many different options (genomes, target\_fastas,
regions, etc.) One cool thing to realize was that all those different,
possible, analysis pathways could be handled by using wildcards. And
once all those wildcards were in place, the actual analyses done would
be determined by what the desired end result is. (More on that later: a
key to all this is figuring out how to elegantly use the config file and
the samples.csv and units.csv to properly expand needed inputs in the
aggregation steps (like VCF creation)).

### Put fast commands in to check the wildcard framework

While developing the workflow, I wanted to start out focusing on making
sure that I had the wildcarding and paths correct, without worrying too
much, at first, about making sure the exact command (i.e., for mapping,
etc.). So, when making new rules, I would just echo what I thought the
command should be into the output file. That way I could check
everything.

Note that snakemake has the `--touch` option, which might do something
similar. But, while developing stuff and making sure the wildcarding
framework is correct, just populating the `shell` blocks with simple
commands like:

``` python
shell:
    "echo bwa mem -R {params.rg} {input.g} {input.EF} > {output.bam}; "
    "touch {output.bai}"
```

seems to be a pretty good way to go.

I will just point out here that by using this approach I have managed to
verify the logic of the workflow, while developing it, without actually
having to get the real data from anyone. I just created some sample
names and genome names, and touched “fastq” files. This has let me
quickly figure out if the logic works, for example, in cases where one
marker set might be typed via a full-genome approach and also via a
target fasta approach.

### Definitely get the `tree` command

`tree` is a Unix utility that plots diagrammatic trees of your directory
structure. You can get it on a Mac via Homebrew with:

``` sh
brew install tree
```

This turns out to be **absolutely indispensable** when it comes to
grokking out the file hierarchies that get created when you run a
snakemake workflow. It is especially helpful during development if your
workflow is particularly large and it is hard to keep track of all the
files that are getting created. Here is a screenshot of my terminal
after a `tree` command.

![images\_for\_readme/tree-grap.png](images_for_readme/tree-grab.png)

### Remember to use Unix file iteration/alternation to specify requested output files

When developing a rule and wondering if the wildcards are set up
properly, it is really helpful to request the exact target file that
should be created by the rule. Remember, as pointed out in the snakemake
manual, that you can use the Unix `{this,that,the_other}` construct to
put those files on the command line. The same goes for putting series of
numbers using `{1..8}`, for example. Thus, while developing this, I can
request certain output files from snakemake like this:

``` sh
snakemake -np \
    .test/data/bams/fullg/{Otsh_v1.0,nookie2}/CH_00{1..8}.bam \
    .test/data/bams/target_fastas/ROSA/{rosa_seqs_1,rosa_seqs_2}/CH_00{1..8}.bam \
    .test/data/bams/target_fastas/TRANSITION/transition-panel/CH_00{1..8}.bam \
    resources/bedfiles/fullg/{LFAR,WRAP,ROSA}-Otsh_v1.0.bed  \
    resources/thinned_genomes/Otsh_v1.0/{LFAR,WRAP,ROSA}/thinned.fa
```

This is good for testing your grammar in those early steps without
having to specifically deal with the aggregation that might come later
(i.e., the VCF-making step).

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

mega-simple-microhap-snakeflow
================

  - [Production runs](#production-runs)
  - [Multi-run variant calling](#multi-run-variant-calling)
      - [Making a VCF for microhaplot after multi-run variant
        calling](#making-a-vcf-for-microhaplot-after-multi-run-variant-calling)
      - [A word on getting those files off a remote
        server](#a-word-on-getting-those-files-off-a-remote-server)
  - [“Post-hoc” Rules](#post-hoc-rules)
      - [multi\_dir\_variant\_calling.smk](#multi_dir_variant_calling.smk)
  - [Development related stuff. Will be cleaned up
    later.](#development-related-stuff.-will-be-cleaned-up-later.)
      - [samples.csv and units.csv](#samples.csv-and-units.csv)
      - [The config file](#the-config-file)
  - [Notes on Development](#notes-on-development)
      - [Some python that is good to
        know](#some-python-that-is-good-to-know)
          - [Flattening lists](#flattening-lists)
      - [Testing python snippets](#testing-python-snippets)
          - [Unpacking dictionary keys and
            values](#unpacking-dictionary-keys-and-values)
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
      - [Handling units and samples](#handling-units-and-samples)
  - [Genomes and target\_fastas](#genomes-and-target_fastas)

This repository holds the Snakefile and associated files for Eric’s
first attempt (still under construction\!) at making a SnakeMake-based
workflow for our microhaplotypes.

It is near completion. To try it out on the test data:

1.  Get this repository with something like `git clone
    https://github.com/eriqande/mega-simple-microhap-snakeflow.git`

2.  Make sure you have Miniconda installed.

3.  Get a full install of a Snakemake conda environment by following the
    directions
    [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

4.  Activate that snakemake conda environment

5.  From within the `mega-simple-microhap-snakeflow` directory give this
    command:
    
    ``` sh
    snakemake --config run_dir=.test/data --configfile config/Chinook/config.yaml  --use-conda -np
    ```
    
    It should spit out a bunch of stuff that ends with something like
    this:
    
    ``` sh
    Job counts:
        count   jobs
        1   all
        2   call_fullg_marker_sets_with_bcftools
        2   call_fullgex_remapped_markers_with_bcftools
        2   call_target_fasta_marker_sets_with_bcftools
        20  extract_reads_from_full_genomes
        33  flash_paired_ends
        20  flashed_fastqs_from_fullg_extracted_bams
        1   make_microhap_folder
        20  map_fullg_extracted_to_thinned_genomes
        20  map_to_full_genome
        26  map_to_target_fastas
        2   microhap_extract_fullgex_remapped
        2   microhap_extract_target_fastas
        151
    This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
    ```

6.  If that worked, you can do a full run with this:
    
    ``` sh
    snakemake --config run_dir=.test/data --configfile config/Chinook/config.yaml  --use-conda --cores 1
    ```
    
    The test data set is small enough that a single core is enough.
    However, the first time you run it, you can expect to spend quite a
    bit of time (30 minute to a couple hours, depending on the speed of
    your computer) downloading the Otsh\_v1.0 genome and indexing it
    with bwa.

## Production runs

Here is an idea of how a production-run command line might look for five
directories. After checking out a node with 20 cores…

``` sh
for i in data/200715_M02749_0092_000000000-CV34F \
    data/201012_M02749_0095_000000000-CWHDK \
    data/210129_M02749_0102_000000000-J33JD \
    data/200828_M02749_0093_000000000-CV2FP \
    data/201019_M02749_0096_000000000-CV34D; do
    snakemake --config run_dir=$i --configfile config/Chinook/config.yaml --use-conda --use-envmodules --cores 20
done
    
```

## Multi-run variant calling

This is our current setup to simply call variants after running multiple
directories. It just using globbing to figure out which BAMs were
produced. It does not parse the sample sheets of each run to figure out
which individuals should have BAMs and request those.

This is only imlemented for the fullgex-remapped-to-thinned stuff, at
the moment, but I will get it going for the target\_fastas soon, too.

Here is how to do it:

1.  Establish a directory in the top level called
    `MULTI_RUN_RESULTS/some_directory_name` where `some_directory_name`
    is some informative name.

2.  Inside that directory make a file called `dirs.txt` which lists the
    run paths to the directories you want to include individuals from
    (if any for the marker sets you will request). The paths should be
    relative to the top-level of the `mega-simple-microhap-snakeflow`
    directory. For example, after the above to steps you might have:
    
    ``` sh
    # A directory
    MULTI_RUN_RESULTS/5_early_runs
    
    # a file dirs.txt within that with the following contents shown
    (snakemake) [node11: mega-simple-microhap-snakeflow]--% cat MULTI_RUN_RESULTS/5_early_runs/dirs.txt
    data/200715_M02749_0092_000000000-CV34F
    data/200828_M02749_0093_000000000-CV2FP
    data/201012_M02749_0095_000000000-CWHDK
    data/201019_M02749_0096_000000000-CV34D
    data/210129_M02749_0102_000000000-J33JD
    
    # that is just 5 lines, each with a path to a previously run run directory!
    ```

3.  Request the output file you want by replacing the wildcards, as
    appropriate, in the following output path:
    
    ``` sh
    MULTI_RUN_RESULTS/{multi_run_dir}/{species_dir}/vcfs/{marker_set}/fullgex_remapped/{genome}/variants-from-multi-runs-bcftools.vcf
    ```
    
    In our example, if we wanted to do this for both WRAP and LFAR, the
    output targets we would request would be:
    
    ``` sh
    MULTI_RUN_RESULTS/5_early_runs/Chinook/vcfs/{LFAR,WRAP}/fullgex_remapped/Otsh_v1.0/variants-from-multi-runs-bcftools.vcf
    ```
    
    Note the use of the Unix brace expansion to get two paths out of
    that—one for LFAR and one for WRAP.

4.  Request that target with snakemake. The full command would look
    like:
    
    ``` sh
    snakemake --use-conda --cores 20  \
    MULTI_RUN_RESULTS/5_early_runs/Chinook/vcfs/{LFAR,WRAP}/fullgex_remapped/Otsh_v1.0/variants-from-multi-runs-bcftools.vcf
    ```

### Making a VCF for microhaplot after multi-run variant calling

The obvious reason for doing multi-run variant calling is to have
sufficient numbers of individuals from enough populations to ensure that
you have most of the interesting variants for forming microhaplotypes.
Here is how you create a new VCF for microhaplot, after doing the above
multi-run variant calling. Here I show some simple steps for the WRAP
markers:

First, filter out sites with a lot of missing data:

``` sh
WRAP=MULTI_RUN_RESULTS/5_early_runs/Chinook/vcfs/WRAP/fullgex_remapped/Otsh_v1.0/variants-from-multi-runs-bcftools.vcf
(snakemake) [node11: mega-simple-microhap-snakeflow]--% bcftools view -i 'F_MISSING < 0.30' $WRAP | awk '!/^#/' | wc
     55   37455  867160
(snakemake) [node11: mega-simple-microhap-snakeflow]--% bcftools view -i 'F_MISSING < 0.10' $WRAP | awk '!/^#/' | wc
     53   36093  824120
(snakemake) [node11: mega-simple-microhap-snakeflow]--% bcftools view -i 'F_MISSING < 0.03' $WRAP | awk '!/^#/' | wc
     53   36093  824120
# OK! Looks like 53 variants in the 24 regions/microhaps

# How many alternate alleles out of how many total alleles at each site?
(snakemake) [node11: mega-simple-microhap-snakeflow]--% bcftools view -i 'F_MISSING < 0.03' $WRAP | bcftools query -f '%CHROM\t%POS\t%AC/%AN\n'
NC_037104.1:55923357-55923657   151 583/1322
NC_037104.1:55966251-55966551   149 632/1322
NC_037104.1:55966251-55966551   151 641/1322
NC_037104.1:56061938-56062238   151 604/1322
NC_037104.1:56061938-56062238   182 1322/1322
NC_037104.1:56088878-56089178   151 452/1322
NC_037108.1:73538966-73539266   151 524/1318
NC_037108.1:73538966-73539266   160 57/1320
NC_037108.1:73538966-73539266   214 56/1320
NC_037108.1:73540716-73541016   151 526/1322
NC_037108.1:73543706-73544006   130 148/1322
NC_037108.1:73543706-73544006   133 17/1322
NC_037108.1:73543706-73544006   151 475/1322
NC_037108.1:73553140-73553440   151 475/1324
NC_037112.1:24500367-24500667   151 468/1320
NC_037112.1:24542569-24542869   151 455/1324
NC_037112.1:24542569-24542869   165 63/1324
NC_037112.1:24542569-24542869   172 891/1324
NC_037112.1:24593758-24594058   151 454/1320
NC_037112.1:24593758-24594058   181 141/1320
NC_037112.1:24609163-24609463   151 874/1322
NC_037112.1:24609163-24609463   188 447/1320
NC_037112.1:24618993-24619293   109 9/1322
NC_037112.1:24618993-24619293   131 455/1322
NC_037112.1:24618993-24619293   151 867/1322
NC_037112.1:24704405-24704705   151 454/1322
NC_037112.1:24704405-24704705   207 29/1322
NC_037112.1:24721041-24721341   141 450/1320
NC_037112.1:24721041-24721341   162 448/1320
NC_037112.1:24999768-25000068   125 455/1320
NC_037112.1:24999768-25000068   147 456/1320
NC_037112.1:24999768-25000068   155 455/1320
NC_037112.1:25015012-25015312   118 532/1322
NC_037112.1:25015012-25015312   151 453/1322
NC_037112.1:25015012-25015312   161 2/1322
NC_037112.1:25015012-25015312   163 801,6/1322
NC_037112.1:28278997-28279297   125 256/1326
NC_037112.1:28278997-28279297   151 508/1324
NC_037112.1:28296342-28296642   130 11/1324
NC_037112.1:28296342-28296642   149 467/1324
NC_037112.1:28296342-28296642   150 77/1308
NC_037112.1:28296342-28296642   151 482/1324
NC_037112.1:28296342-28296642   165 15/1324
NC_037112.1:28320487-28320787   126 53/1322
NC_037112.1:28320487-28320787   128 39/1322
NC_037112.1:28320487-28320787   137 831/1318
NC_037112.1:28320487-28320787   151 540/1322
NC_037112.1:28350510-28350810   116 22/1322
NC_037112.1:28350510-28350810   141 1170/1322
NC_037112.1:28350510-28350810   151 530/1322
NC_037112.1:28350510-28350810   163 363/1322
NC_037121.1:6243491-6243791 151 474/1324
NC_037121.1:6268440-6268740 151 448/1318

# So, some of them are quite rare, but let's go ahead and keep
# those in there...That is still fewer than two variants per
# amplicon, on average.

# Finally, let's make a VCF file with only one individual in it
# to use for our microhaplot VCF:
(snakemake) [node11: mega-simple-microhap-snakeflow]--% bcftools view -i 'F_MISSING < 0.03' $WRAP | bcftools view -s T170774  > config/Chinook/canonical_variation/WRAP-24-amplicons-53-variants.vcf
```

Now, to add that into the workflow, let’s request another set of
canonical variation in the Chinook config:

``` yaml
WRAP:
    genome:
      Otsh_v1.0:
        regions: config/Chinook/regions/WRAP-Otsh_v1.0.txt
        microhap_variants:
          all_variants: config/Chinook/canonical_variation/WRAP-all-snps-round-1.vcf
          after_5_runs: config/Chinook/canonical_variation/WRAP-24-amplicons-53-variants.vcf  <--- this line added
  
```

After that, invoking Snakemake will see that there are new microhap
variants for things to be run at, if there are any WRAP fish in the run.
Cool\!

### A word on getting those files off a remote server

If we just did 5 runs worth of microhaplot for LFAR and WRAP, the
results would be found in directories like this:

``` sh
(snakemake) [node11: mega-simple-microhap-snakeflow]--% pwd
/home/eanderson/Documents/git-repos/mega-simple-microhap-snakeflow
(snakemake) [node11: mega-simple-microhap-snakeflow]--% ls data/*/Chinook/microhaplot/*after_5*
data/200715_M02749_0092_000000000-CV34F/Chinook/microhaplot/LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
data/200715_M02749_0092_000000000-CV34F/Chinook/microhaplot/LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
data/200828_M02749_0093_000000000-CV2FP/Chinook/microhaplot/LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
data/200828_M02749_0093_000000000-CV2FP/Chinook/microhaplot/LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
data/200828_M02749_0093_000000000-CV2FP/Chinook/microhaplot/WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
data/200828_M02749_0093_000000000-CV2FP/Chinook/microhaplot/WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
data/201012_M02749_0095_000000000-CWHDK/Chinook/microhaplot/LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
data/201012_M02749_0095_000000000-CWHDK/Chinook/microhaplot/LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
data/201019_M02749_0096_000000000-CV34D/Chinook/microhaplot/WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
data/201019_M02749_0096_000000000-CV34D/Chinook/microhaplot/WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
data/210129_M02749_0102_000000000-J33JD/Chinook/microhaplot/LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
data/210129_M02749_0102_000000000-J33JD/Chinook/microhaplot/LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
data/210129_M02749_0102_000000000-J33JD/Chinook/microhaplot/WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
data/210129_M02749_0102_000000000-J33JD/Chinook/microhaplot/WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
```

They all have the same name, so it is the directory structure that
distinguishes them. To get these all to your laptop, rsync with the R
option is your friend (and you can test with the `-n` dry-run option
before doing it):

``` sh
 rsync -avR eanderson@sedna.nwfsc2.noaa.gov:'/home/eanderson/Documents/git-repos/mega-simple-microhap-snakeflow/./data/*/Chinook/microhaplot/*-after_5_*' ./
```

Note the use of the single quotes to avoid expanding the wildcards on
the local machine, and instead send them to the remote machine. And
*also* note the `.` in: `mega-simple-microhap-snakeflow/./data`. That
tells rsync to only include the part of the path to the left of the dot.
Cool\! After running the above command, we have:

``` sh
(base) /from_cluster/--% (master)  tree .
.
└── data
    ├── 200715_M02749_0092_000000000-CV34F
    │   └── Chinook
    │       └── microhaplot
    │           ├── LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
    │           └── LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
    ├── 200828_M02749_0093_000000000-CV2FP
    │   └── Chinook
    │       └── microhaplot
    │           ├── LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
    │           ├── LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
    │           ├── WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
    │           └── WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
    ├── 201012_M02749_0095_000000000-CWHDK
    │   └── Chinook
    │       └── microhaplot
    │           ├── LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
    │           └── LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
    ├── 201019_M02749_0096_000000000-CV34D
    │   └── Chinook
    │       └── microhaplot
    │           ├── WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
    │           └── WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
    └── 210129_M02749_0102_000000000-J33JD
        └── Chinook
            └── microhaplot
                ├── LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
                ├── LFAR--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds
                ├── WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs.rds
                └── WRAP--fullgex_remapped_to_thinned--Otsh_v1.0--after_5_runs_posinfo.rds

16 directories, 14 files
```

# “Post-hoc” Rules

## multi\_dir\_variant\_calling.smk

This is a recently added rule that lets one quickly call variants from
muliple runs together. The inspiration for this came from the
possibility of a situtation where run A might include individuals from
only one population of a species and run B might include individuals
from only a different population. If these pops happen to be fixed for
different variants, then those will never show up in the VCFs, unless
all them are run together. So, this rule is there to allow the user to
easily call variation from previously created BAMs. It *does not* use
any sample sheet information to figure out which individuals should be
included. Rather it just globs the existing BAM files for each marker
set that have already been run. Not only that, but it currently requires
the user to set up the directory for the outputs to be run in, and also
a file that lists the input run directories. And the user has to request
the output target with the correct file path, etc.

# Development related stuff. Will be cleaned up later.

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

    ## # A tibble: 33 x 8
    ##    sample Marker_Sets NMFS_DNA_ID Sample_ID Sample_Name Sample_Project fq1  
    ##    <chr>  <chr>       <chr>       <chr>     <chr>       <chr>          <chr>
    ##  1 s0001  TRANSITION… T194879     CH_33547  CH-33547    CH_microhaps_… CH-3…
    ##  2 s0002  TRANSITION… T194887     CH_33548  CH-33548    CH_microhaps_… CH-3…
    ##  3 s0003  TRANSITION… T194895     CH_33549  CH-33549    CH_microhaps_… CH-3…
    ##  4 s0004  TRANSITION… T194903     CH_33550  CH-33550    CH_microhaps_… CH-3…
    ##  5 s0005  TRANSITION… T194911     CH_33551  CH-33551    CH_microhaps_… CH-3…
    ##  6 s0006  TRANSITION… T194919     CH_33552  CH-33552    CH_microhaps_… CH-3…
    ##  7 s0007  TRANSITION… T194927     CH_33553  CH-33553    CH_microhaps_… CH-3…
    ##  8 s0008  TRANSITION… T194935     CH_33554  CH-33554    CH_microhaps_… CH-3…
    ##  9 s0009  TRANSITION… T194943     CH_33555  CH-33555    CH_microhaps_… CH-3…
    ## 10 s0010  TRANSITION… T194951     CH_33556  CH-33556    CH_microhaps_… CH-3…
    ## # … with 23 more rows, and 1 more variable: fq2 <chr>

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

    ## # A tibble: 46 x 2
    ##    sample Markers   
    ##    <chr>  <chr>     
    ##  1 s0001  TRANSITION
    ##  2 s0001  ROSA      
    ##  3 s0002  TRANSITION
    ##  4 s0002  ROSA      
    ##  5 s0003  TRANSITION
    ##  6 s0003  ROSA      
    ##  7 s0004  TRANSITION
    ##  8 s0004  ROSA      
    ##  9 s0005  TRANSITION
    ## 10 s0005  ROSA      
    ## # … with 36 more rows

## The config file

This file is in `config/newconfig.yaml`. So far it looks like this:

``` yaml
species: Chinook

run_dir: .test/data  # set it on the command line to the directory in which raw, samples.csv, and units.csv reside

# the following are names of files that are assumed to be
# within the run_dir.  NOTE! These names might be hardwired
# in extract_microhaplotypes.R.  Don't change these...
samples: samples.csv
units: units.csv

# set genomes up to be downloaded from url.  They will
# be placed in resources/{species_dir}/genomes/name_of_genome/name_of_genome.fna
# where name_of_genome is like Otsh_v1.0
genomes:
  Otsh_v1.0:
    url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/872/995/GCF_002872995.1_Otsh_v1.0/GCF_002872995.1_Otsh_v1.0_genomic.fna.gz


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

# But, also note that we can have multiple variant sets for each genome
# within a genome-focused marker set.

# The "variants" field is reserved only for the canonical SNP-sets that we would
# use to call specific variants.  Not for finding new variants, etc.
marker_sets:
  LFAR:
    genome:
      Otsh_v1.0:
        regions: config/Chinook/regions/LFAR-Otsh_v1.0.txt
        microhap_variants:
          all_variants: config/Chinook/canonical_variation/LFAR-all-snps-round-1.vcf
          after_5_runs: config/Chinook/canonical_variation/LFAR-8-amplicons-22-variants.vcf
  WRAP:
    genome:
      Otsh_v1.0:
        regions: config/Chinook/regions/WRAP-Otsh_v1.0.txt
        microhap_variants:
          all_variants: config/Chinook/canonical_variation/WRAP-all-snps-round-1.vcf
          after_5_runs: config/Chinook/canonical_variation/WRAP-24-amplicons-53-variants.vcf
  ROSA:
    target_fasta:
      rosawr:
        fasta: config/Chinook/target_fastas/greb1_rosa_with_wrdiag_reference.fasta
        microhap_variants:
          all_variants: config/Chinook/canonical_variation/ROSA-rosawr-all-snps-round-1.vcf
  TRANSITION:
    target_fasta:
      transition_204:
        fasta: config/Chinook/target_fastas/chinook_204amps_transition.fasta
        microhap_variants:
          snplicons: config/Chinook/canonical_variation/chinook_204amps_ultima4_transition_multi_targets.vcf
  SNPS84:
    target_fasta:
      snps84:
        fasta: config/Chinook/target_fastas/chinook_84amps.fasta
        microhap_variants:
          single_snps: config/Chinook/canonical_variation/chinook_84amps_single_targets.vcf











```

# Notes on Development

## Some python that is good to know

### Flattening lists

Since the `+` operator catenates elements of lists together:

``` python
>>> a = [1]; b = [2]; a + b
[1, 2]
```

You can also use sum to flatten lists:

``` python
>>> lumpy = [[1, 2, 3], [8, 9, 10], [4], [5]]
>>> lumpy
[[1, 2, 3], [8, 9, 10], [4], [5]]

# the second argument here is the starting list to start catenating to.
# So we pass it an empty list...
>>> sum(lumpy, [])
[1, 2, 3, 8, 9, 10, 4, 5]
```

## Testing python snippets

For someone like me who is not super familiar with python, it is nice to
be able to test various constructions in the python interpreter. Here
are some fun tricks to allows that to happen.

### Unpacking dictionary keys and values

If you have several genomes, for example, in your config:

``` python
[*config["genome].keys()]
```

The `*` “unpacks” the `odict_keys` object into the enclosing list
environment.

Note that you might also use `list(config["genome"].keys())`.

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

Note that if you want to see what all the current values are in the
wildcards object that you have defined for testing, you can do that
with:

``` python
 wildcards.__dict__
```

which will produce output like:

    {'marker_set': 'TRANSITION', 'genome': 'Otsh_v1.0', 'target_fasta': 'transition-panel', 'run_dir': '.test/data'}

#### Using values that caused things to fail

Sometimes you get a message like this when there is a failure:

``` sh
InputFunctionException in line 24 of /Users/eriq/Documents/git-repos/mega-simple-microhap-snakeflow/workflow/rules/full_genome_map_and_extract.smk:
Error:
  KeyError: 'nookie2'
Wildcards:
  run_dir=.test/data
  marker_set=ROSA
  genome=nookie2
  sample=CH_002
Traceback:
  File "/Users/eriq/Documents/git-repos/mega-simple-microhap-snakeflow/workflow/rules/common.smk", line 40, in region_files_from_marker_set_and_genome
```

It would be nice to set a local wildcards variable to have those values.
Here is a script that does that, if you first copy out these lines and
have them in your clipboard:

``` 
  run_dir=.test/data
  marker_set=ROSA
  genome=nookie2
  sample=CH_002
```

Then run this script:

``` sh
# if the error values of snakemake wildcards are on your clipboard like this:
#
#  run_dir=.test/data
#  marker_set=ROSA
#  genome=nookie2
#  sample=CH_002
#
# Then, run this script and it will put the python code to assign such
# values to a local wildcards variable onto your clipboard

pbpaste | sed 's/^ *//g; s/ *$//g;' |  awk -F"=" '{printf("wildcards.%s = \"%s\";\n", $1,$2)}' | pbcopy
```

i.e., put that into a script called `wilcard_em.sh` in your PATH, and
then invoke it and your clipboard will then have:

``` python
wildcards.run_dir = ".test/data";
wildcards.marker_set = "ROSA";
wildcards.genome = "nookie2";
wildcards.sample = "CH_002";
```

which you can paste into python for testing.

## Getting filepaths and wildcards right

Just understanding the wildcarding framework of snakemake can take a
while to get your head around. This microha workflow is pretty beastly
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
    resources/thinned_genomes/Otsh_v1.0/{LFAR,WRAP,ROSA}/thinned.fa \
    .test/data/bams/fullg-extracted/{LFAR,WRAP,ROSA}/Otsh_v1.0/CH_00{1..8}.bam \
    .test/data/bams/fullg-extracted/LFAR/nookie2/CH_00{1..8}.bam
```

This is good for testing your grammar/logic in those early steps without
having to specifically deal with the aggregation that might come later
(i.e., the VCF-making step).

The above block was great for seeing if all my bam creation and
extraction steps were solid.

Now, to see if my rules that involve expanding over units to create the
vcfs are working, I might try this:

``` sh
snakemake -np \
    .test/data/vcfs/LFAR/fullg/{Otsh_v1.0,nookie2}/variants-bcftools.vcf
```

## Handling units and samples

I could get whatever I need from these in a couple lines in R, but
snakemake is done in python, so, I need to figure that out. I am sure
there are muliple ways to do all these things, I am just going to find
one that works.

I decided it would be easiest to break units into two pandas data
frames: one for the Markers that had target fastas, and another for the
Markers that were genome-focused.

``` python
# code like this goes into common.smk
tf_units = units[tflist == True]
tflist = units["Markers"].isin(config["target_fasta_focused_marker_sets"]["name"])
```

I need to cycle over the rows in units, and for each one return a path
for combination of:

  - sample
  - marker\_set
  - genome of marker\_set (if any)
  - target\_fasta of marker set (if any)

Note that for testing things in a python session, it might be useful to
do:

``` python
from snakemake.io import expand
```

Then we should be able to iterate over the rows in units like this:

``` python
for index, row in tf_units.iterrows():
  print(row['Sample_ID'], row['Markers'])

for index, row in gf_units.iterrows():
  print(row['Sample_ID'], row['Markers'])

```

So, we just have to figure out how to expand those over the possible
genomes and target fastas…

``` python
ret = list()
for index, row in tf_units.iterrows():
  S = row['Sample_ID']
  M = row['Markers']
  ret = ret + expand("{M}--{k}--{S}", M = M, S = S, k =  [str(k) for k in config["marker_sets"][M]["target_fasta"].keys()] )

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

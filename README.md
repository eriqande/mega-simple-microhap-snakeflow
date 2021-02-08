mega-simple-microhap-snakeflow
================

This repository holds the Snakefile and associated files for Eric’s
first attempt (still under construction\!) at making a SnakeMake-based
workflow for our microhaplotypes.

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

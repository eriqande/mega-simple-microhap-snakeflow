configfile: "config.yaml"

# make a bedfile for our target regions from the full genome
rule region_bedfiles_full:
  params:
    reg=expand("{region}", region=config["regions"])
  output:
    "bedfiles/regions_full.bed"
  shell:
    "echo {params.reg} | "
    "awk '{{for(i=1;i<=NF;i++) print $i}}' | "
    "sed 's/[-:]/ /g' | "
    "awk 'BEGIN {{OFS=\"\t\"}} {{print $1, $2-1, $3, $1\":\"$2\"-\"$3}}' > {output}"



# we make a bedfile that we can use for checking the depth
# on the reads mapped to the thin reference. We pull this from
# the very first sample's bamfile
rule region_bedfiles_thin:
  input:
    one=expand("mapped_reads_thin/{sample}.bam", sample=config["samples"])[0]
  output:
    "bedfiles/regions_thin.bed"
  shell:
    "samtools view -H {input.one} | awk '/^@SQ/ {{print $2, $3}}' | "
    "sed 's/SN://g; s/LN://g;' | "
    "awk 'BEGIN {{OFS=\"\t\"}} {{print $1, 0, $2, $1}}' > {output} "



# pull the target regions out of the full genome and make a "thinned reference"
# from them.  We can see how many amplicons are getting mapped to non-target
# regions by comparing the number of reads mapped to each fragment in this thinned
# genome to the number of reads mapped to the target region when mapping to
# the whole genome.
rule make_thinned_genome:
  input:
    fa=config["genome"]
  params:
    regs=expand("{region}", region=config["regions"])
  output:
    fa="thinned_genome/thinned.fa",
    fai="thinned_genome/thinned.fa.fai",
    idx="thinned_genome/thinned.fa.bwt"
  shell:
    "samtools faidx {input.fa} {params.regs} > {output.fa}; "
    "samtools faidx {output.fa}; "
    "bwa index {output.fa}"



# map all the reads to the full genome
rule bwa_map:
  input:
    config["genome"],
    lambda wildcards: config["samples"][wildcards.sample]
  output:
    bam="mapped_reads/{sample}.bam",
    bai="mapped_reads/{sample}.bam.bai"
  params:
    rg=r"@RG\tID:{sample}\tSM:{sample}"
  log:
    "logs/bwa_map/{sample}.log"
  benchmark:
    "benchmarks/bwa_map/{sample}.benchmark.txt"
  shell:
    "(bwa mem -R '{params.rg}' {input} | samtools view -Sb - | "
    "samtools sort -T mapped_reads/{wildcards.sample} -O bam - > {output.bam}; "
    "samtools index {output.bam}) 2> {log}"



# map all the reads to the thinned genome
rule bwa_map_thinned:
  input:
    "thinned_genome/thinned.fa",
    lambda wildcards: config["samples"][wildcards.sample]
  output:
    "mapped_reads_thin/{sample}.bam"
  params:
    rg=r"@RG\tID:{sample}\tSM:{sample}"
  log:
    "logs/bwa_map_thinned/{sample}.log"
  benchmark:
    "benchmarks/bwa_map_thinned/{sample}.benchmark.txt"
  shell:
    "(bwa mem -R '{params.rg}' {input} | samtools view -Sb - | "
    "samtools sort -T mapped_reads_thin/{wildcards.sample} -O bam - > {output}) 2> {log}"



# extract the reads overlapping our target regions when mapped
# to the full genome
rule extract_regions:
  input:
    "mapped_reads/{sample}.bam"
  output:
    bam="extracted_reads/{sample}.bam",
    bai="extracted_reads/{sample}.bam.bai"
  params:
    regs=expand("{region}", region=config["regions"])
  shell:
    "samtools view -u {input} {params.regs} | "
    "samtools sort -T extracted_reads/{wildcards.sample} -O bam - > {output.bam}; "
    "samtools index {output.bam}"



# compute the coverage at the target regions for each sample
# when mapped to the full genome
rule coverage_full:
  input:
    a="bedfiles/regions_full.bed",
    b="extracted_reads/{sample}.bam"
  output:
    "coverage/samples/{sample}_cov_full.tsv"
  params:
    SM="{sample}"
  shell:
    "bedtools coverage -a {input.a} -b {input.b} | "
    "awk -v SM={params.SM} '{{printf(\"full\t%s\t%s\\n\", SM, $0);}}' > {output}"


# compute the coverage at the target regions for each sample
# when mapped to the full genome
rule coverage_thin:
  input:
    a="bedfiles/regions_thin.bed",
    b="mapped_reads_thin/{sample}.bam"
  output:
    "coverage/samples/{sample}_cov_thin.tsv"
  params:
    SM="{sample}"
  shell:
    "bedtools coverage -a {input.a} -b {input.b} | "
    "awk -v SM={params.SM} '{{printf(\"thin\t%s\t%s\\n\", SM, $0);}}' > {output}"
 


# Now catenate all the full and thin coverges into a single
# tsv file that we can analyze easily in the tidyverse
rule catenate_coverages:
  input:
    expand("coverage/samples/{sample}_cov_{cond}.tsv", sample=config["samples"], cond=["full", "thin"])
  output:
    "coverage/all_coverages.tsv"
  shell:
    "(echo \"genome_condition\tsample\tchrom\tstart\tstop\ttarget\tnum_reads\tnum_bases\tfeature_length\tfract_bases\"; "
    "cat {input}) > {output}"

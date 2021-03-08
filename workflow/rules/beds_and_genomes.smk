

# we make a bedfile that we can use for checking the depth
# on the reads mapped to the thin reference. We pull this from
# the very first sample's bamfile
rule region_bedfiles_thin:
  input:
    one=expand("results/mapped_reads_thin/{sample}.bam", sample=config["samples"])[0]
  output:
    "results/bedfiles/regions_thin.bed"
  envmodules:
    "bio/samtools"
  conda:
    "../envs/bwasam.yaml"
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
    fa="results/thinned_genome/thinned.fa",
    fai="results/thinned_genome/thinned.fa.fai",
    idx="results/thinned_genome/thinned.fa.bwt"
  envmodules:
    "aligners/bwa",
    "bio/samtools"
  conda:
    "../envs/bwasam.yaml"
  shell:
    "samtools faidx {input.fa} {params.regs} > {output.fa}; "
    "samtools faidx {output.fa}; "
    "bwa index {output.fa}"


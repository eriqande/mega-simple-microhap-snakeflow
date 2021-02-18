# map all the reads to the full genome
rule bwa_map:
  input:
    config["genome"],
    lambda wildcards: config["samples"][wildcards.sample]
  output:
    bam="results/mapped_reads/{sample}.bam",
    bai="results/mapped_reads/{sample}.bam.bai"
  params:
    rg=r"@RG\tID:{sample}\tSM:{sample}"
  log:
    "results/logs/bwa_map/{sample}.log"
  benchmark:
    "results/benchmarks/bwa_map/{sample}.benchmark.txt"
  envmodules:
    "aligners/bwa",
    "bio/samtools"
  shell:
    "(bwa mem -R '{params.rg}' {input} | samtools view -Sb - | "
    "samtools sort -T mapped_reads/{wildcards.sample} -O bam - > {output.bam}; "
    "samtools index {output.bam}) 2> {log}"



# map all the reads to the thinned genome
rule bwa_map_thinned:
  input:
    "results/thinned_genome/thinned.fa",
    lambda wildcards: config["samples"][wildcards.sample]
  output:
    "results/mapped_reads_thin/{sample}.bam"
  params:
    rg=r"@RG\tID:{sample}\tSM:{sample}"
  log:
    "results/logs/bwa_map_thinned/{sample}.log"
  benchmark:
    "results/benchmarks/bwa_map_thinned/{sample}.benchmark.txt"
  envmodules:
    "aligners/bwa",
    "bio/samtools"
  shell:
    "(bwa mem -R '{params.rg}' {input} | samtools view -Sb - | "
    "samtools sort -T mapped_reads_thin/{wildcards.sample} -O bam - > {output}) 2> {log}"



# extract the reads overlapping our target regions when mapped
# to the full genome
rule extract_regions:
  input:
    "results/mapped_reads/{sample}.bam"
  output:
    bam="results/extracted_reads/{sample}.bam",
    bai="results/extracted_reads/{sample}.bam.bai"
  params:
    regs=expand("{region}", region=config["regions"])
  envmodules:
    "bio/samtools"
  shell:
    "samtools view -u {input} {params.regs} | "
    "samtools sort -T extracted_reads/{wildcards.sample} -O bam - > {output.bam}; "
    "samtools index {output.bam}"


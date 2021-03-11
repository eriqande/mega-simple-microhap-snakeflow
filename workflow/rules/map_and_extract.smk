
#### THIS SECTION IS MAPPING TO FULL GENOMES, THEN EXTRACTING  ####

# map reads to the requested full genome
rule map_to_full_genome:
  input:
    g=fna_from_genome,
    bwt=fna_bwt_from_genome,
    EF="{run_dir}/flash/{sample}-extendedFrags.fna.gz"
  output:
    bam="{run_dir}/bams/fullg/{genome}/{sample}.bam",
    bai="{run_dir}/bams/fullg/{genome}/{sample}.bam.bai"
  params:
    rg=r"@RG\tID:{sample}\tSM:{sample}"
  log:
    "{run_dir}/logs/map_to_full_genome/{genome}/{sample}.log"
  conda:
    "../envs/bwasam.yaml"
  shell:
    "echo bwa mem -R {params.rg} {input.g} {input.EF} > {output.bam}; "
    "touch {output.bai}"


# extract the reads overlapping our target regions when mapped
# to the full genome
rule extract_reads_from_full_genomes:
  input:
    bam="{run_dir}/bams/fullg/{genome}/{sample}.bam",
    regfile=region_files_from_marker_set_and_genome
  output:
    bam="{run_dir}/bams/fullg-extracted/{marker_set}/{genome}/{sample}.bam",
    bai="{run_dir}/bams/fullg-extracted/{marker_set}/{genome}/{sample}.bam.bai"
  envmodules:
    "bio/samtools"
  log:
    "{run_dir}/logs/extract_reads_from_full_genomes/{marker_set}/{genome}/{sample}.log"
  conda:
    "../envs/bwasam.yaml"
  shell:
    "echo samtools view -u {input.bam} $(cat {input.regfile}) > {output.bam} 2> {log}; touch {output.bai} "
#   "samtools sort -T extracted_reads/{wildcards.sample} -O bam - > {output.bam}; "
#   "samtools index {output.bam}"




#### THIS SECTION IS MAPPING TO TARGET FASTAS WITH NO EXTRACTING NECESSARY ####


# map reads to the requested full genome
rule map_to_target_fastas:
  input:
    g=fna_from_marker_set_and_target_fasta,
    EF="{run_dir}/flash/{sample}-extendedFrags.fna.gz"
  output:
    bam="{run_dir}/bams/target_fastas/{marker_set}/{target_fasta}/{sample}.bam",
    bai="{run_dir}/bams/target_fastas/{marker_set}/{target_fasta}/{sample}.bam.bai"
  params:
    rg=r"@RG\tID:{sample}\tSM:{sample}"
  log:
    "{run_dir}/logs/map_to_target_fastas/{marker_set}/{target_fasta}/{sample}.log"
  conda:
    "../envs/bwasam.yaml"
  shell:
    "echo bwa mem -R {params.rg} {input.g} {input.EF} > {output.bam} 2> {log}; "
    "touch {output.bai}"

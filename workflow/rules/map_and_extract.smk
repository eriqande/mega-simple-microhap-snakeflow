
# To be expedient, after making BAMs, we then, as part of the same
# block, convert them to SAMs (for microhaplot) and also run idxstats
# on them.  This is not super fine-grained the way a Snakemake workflow
# would typically be, but it doesn't take much time for these other steps
# and we are going to want them for microhaplot in almost all cases.


#### THIS SECTION IS MAPPING TO FULL GENOMES, THEN EXTRACTING  ####

# map reads to the requested full genome
rule map_to_full_genome:
  input:
    g=fna_from_genome,
    bwt=fna_bwt_from_genome,
    EF="{run_dir}/flash/{sample}.extendedFrags.fastq.gz"
  params:
    rg=rg_from_sample
  log:
    bwa="{run_dir}/logs/map_to_full_genome/{genome}/{sample}.bwa.log",
    samtools="{run_dir}/logs/map_to_full_genome/{genome}/{sample}.samtools.log",
  conda:
    "../envs/bwasam.yaml"
  output:
    bam="{run_dir}/bams/fullg/{genome}/{sample}.bam",
    bai="{run_dir}/bams/fullg/{genome}/{sample}.bam.bai"
  shell:
    " bwa mem -R '{params.rg}' {input.g} {input.EF} 2> {log.bwa} | "
    " samtools view -u -  2> {log.samtools} |  "
    " samtools sort -T {wildcards.run_dir}/bams/fullg/{wildcards.genome}/{wildcards.sample} "
    "   -O bam -o {output.bam} - 2>> {log.samtools}; "
    " samtools index {output.bam} 2>> {log.samtools}"


# extract the reads overlapping our target regions when mapped
# to the full genome
rule extract_reads_from_full_genomes:
  input:
    bam="{run_dir}/bams/fullg/{genome}/{sample}.bam",
    regfile=region_files_from_marker_set_and_genome
  log:
    "{run_dir}/logs/extract_reads_from_full_genomes/{marker_set}/{genome}/{sample}.log"
  envmodules:
    "bio/samtools"
  conda:
    "../envs/bwasam.yaml"
  output:
    bam="{run_dir}/bams/fullg-extracted/{marker_set}/{genome}/{sample}.bam",
    bai="{run_dir}/bams/fullg-extracted/{marker_set}/{genome}/{sample}.bam.bai",
    sam="{run_dir}/sams/fullg-extracted/{marker_set}/{genome}/{sample}.sam",
    idx="{run_dir}/idxstats/fullg-extracted/{marker_set}/{genome}/{sample}_idxstats.txt",
  shell:
    " samtools view -u {input.bam} $(cat {input.regfile})  2> {log} |  "
    " samtools sort -T {wildcards.run_dir}/bams/fullg-extracted/{wildcards.marker_set}/{wildcards.genome}/{wildcards.sample} "
    "   -O bam -o {output.bam} -  2>> {log}; "
    " samtools index {output.bam}; "
    " samtools view {output.bam} > {output.sam} 2>> {log}; "
    " samtools idxstats {output.bam} > {output.idx} 2>> {log}; "




#### THIS SECTION IS MAPPING TO TARGET FASTAS WITH NO EXTRACTING NECESSARY ####


# map reads to the requested target_fastas
rule map_to_target_fastas:
  input:
    g="resources/target_fastas/{marker_set}/{target_fasta}/ref.fna",
    EF="{run_dir}/flash/{sample}.extendedFrags.fastq.gz"
  params:
    rg=rg_from_sample
  log:
    bwa="{run_dir}/logs/map_to_target_fastas/{marker_set}/{target_fasta}/{sample}.bwa.log",
    samtools="{run_dir}/logs/map_to_target_fastas/{marker_set}/{target_fasta}/{sample}.samtools.log"
  conda:
    "../envs/bwasam.yaml"
  output:
    bam="{run_dir}/bams/target_fastas/{marker_set}/{target_fasta}/{sample}.bam",
    bai="{run_dir}/bams/target_fastas/{marker_set}/{target_fasta}/{sample}.bam.bai",
    sam="{run_dir}/sams/target_fastas/{marker_set}/{target_fasta}/{sample}.sam",
    idx="{run_dir}/idxstats/target_fastas/{marker_set}/{target_fasta}/{sample}_idxstats.txt",
  shell:
    "bwa mem -R '{params.rg}' {input.g} {input.EF}  2> {log.bwa} | "
    " samtools view -u -  2> {log.samtools} |  "
    " samtools sort -T {wildcards.run_dir}/bams/target_fastas/{wildcards.marker_set}/{wildcards.target_fasta}/{wildcards.sample} "
    "   -O bam -o {output.bam} - 2>> {log.samtools}; "
    " samtools index {output.bam} 2>> {log.samtools}"
    " samtools view {output.bam} > {output.sam} 2>> {log}; "
    " samtools idxstats {output.bam} > {output.idx} 2>> {log}; "


#### THIS SECTION IS MAPPING TO THINNED GENOMES 
# this is mostly used to see how much extra stuff is being amplified that
# is off-target.

# map reads to the requested thinned genomes
rule map_to_thinned_genomes:
  input:
    g="resources/thinned_genomes/{genome}/{marker_set}/thinned.fna",
    EF="{run_dir}/flash/{sample}.extendedFrags.fastq.gz"
  params:
    rg=rg_from_sample
  log:
    bwa="{run_dir}/logs/map_to_thinned_genomes/{genome}/{marker_set}/{sample}.bwa.log",
    samtools="{run_dir}/logs/map_to_thinned_genomes/{genome}/{marker_set}/{sample}.samtools.log"
  conda:
    "../envs/bwasam.yaml"
  output:
    bam="{run_dir}/bams/thinned_genomes/{genome}/{marker_set}/{sample}.bam",
    bai="{run_dir}/bams/thinned_genomes/{genome}/{marker_set}/{sample}.bam.bai"
  shell:
    "bwa mem -R '{params.rg}' {input.g} {input.EF}  2> {log.bwa} | "
    " samtools view -u -  2> {log.samtools} |  "
    " samtools sort -T {wildcards.run_dir}/bams/thinned_genomes/{wildcards.genome}/{wildcards.marker_set}/{wildcards.sample} "
    "   -O bam -o {output.bam} - 2>> {log.samtools}; "
    " samtools index {output.bam} 2>> {log.samtools}"

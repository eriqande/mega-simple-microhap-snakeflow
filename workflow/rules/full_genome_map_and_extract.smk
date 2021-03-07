

# map reads to the requested full genome
rule map_to_full_genome:
  input:
    g=fna_from_genome,
    fq1=fq1_from_sample_and_run,
    fq2=fq2_from_sample_and_run
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
    "echo bwa mem {input.g} {input.fq1} {input.fq2} > {output.bam}; "
    "touch {output.bai}"



# map reads to the requested full genome
rule map_to_full_genome:
  input:
    g=fna_from_genome,
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




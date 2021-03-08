



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

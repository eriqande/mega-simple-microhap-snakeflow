

# every one of the read-pair files gets flashed as a first step.
# We might want to change this up for the things that get mapped
# to a full genome.  It might be better for them to get mapped
# as paired reads.
rule flash_paired_ends:
  input:
    fq1=fq1_from_sample_and_run,
    fq2=fq2_from_sample_and_run
  output:
    EF="{run_dir}/flash/{sample}-extendedFrags.fna.gz",
  log:
    "{run_dir}/logs/flash/{sample}.log"
  conda:
    "../envs/bwasam.yaml"
  shell:
    "echo flash {input.fq1} {input.fq2} > {output.EF} 2> {log}"


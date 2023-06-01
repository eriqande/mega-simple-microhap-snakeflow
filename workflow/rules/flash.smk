

# every one of the read-pair files gets flashed as a first step.
# We might want to change this up for the things that get mapped
# to a full genome.  It might be better for them to get mapped
# as paired reads.
rule flash_paired_ends:
  input:
    fq1=fq1_or_trim1_from_sample_and_run,
    fq2=fq2_or_trim2_from_sample_and_run
  output:
    EF="{run_dir}/{species_dir}/flash/{sample}.extendedFrags.fastq.gz",
  log:
    stdout="{run_dir}/{species_dir}/logs/flash/{sample}.stdout.log",
    stderr="{run_dir}/{species_dir}/logs/flash/{sample}.stderr.log"
  conda:
    "../envs/flash.yaml"
  shell:
    " if [ {input.fq1} = {input.fq2} ]; then  "
    "    cp {input.fq1} {output.EF}; "
    " else          "  
    "   flash -m 4 -M 100 -z -O  "
    "     --output-prefix={wildcards.run_dir}/{wildcards.species_dir}/flash/{wildcards.sample} "
    "     {input.fq1} {input.fq2} > {log.stdout}  2> {log.stderr}; "
    " fi "


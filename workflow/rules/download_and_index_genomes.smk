


# download the requested genome and rename it as such
rule download_genome:
  params:
    url=genome_url_from_genome
  log:
    wget="resources/{species_dir}/logs/download_genome/wget-{genome}.log",
    gunzip="resources/{species_dir}/logs/download_genome/gunzip-{genome}.log"
  conda:
    "../envs/wget.yaml"
  output:
    fna="resources/{species_dir}/genomes/{genome}/{genome}.fna"
  shell:
    "wget -O {output.fna}.gz {params.url} 2> {log.wget}; "
    " gunzip {output.fna}.gz  2> {log.gunzip}"

rule faidx_genome:
  input:
    fna="resources/{species_dir}/genomes/{genome}/{genome}.fna"
  log:
    "resources/{species_dir}/logs/faidx_genome/{genome}.log"
  conda:
    "../envs/samtools.yaml"
  output:
    "resources/{species_dir}/genomes/{genome}/{genome}.fna.fai"
  shell:
    "samtools faidx {input.fna} 2> {log}"


rule bwt_index_genome:
  input:
    fna="resources/{species_dir}/genomes/{genome}/{genome}.fna"
  log:
    "resources/{species_dir}/logs/bwt_index_genome/{genome}.log"
  conda:
    "../envs/bwa.yaml"
  output:
    multiext("resources/{species_dir}/genomes/{genome}/{genome}.fna", ".amb", ".ann", ".bwt", ".pac", ".sa")
  shell:
    "bwa index {input.fna}  2> {log} "



# we also need to bwa index the target fastas
rule copy_and_index_target_fasta:
  input:
    fna=fna_from_marker_set_and_target_fasta
  log:
    "resources/{species_dir}/logs/bwt_index_target_fasta/{marker_set}/{target_fasta}.log"
  conda:
    "../envs/bwasam.yaml"
  output:
    fna="resources/{species_dir}/target_fastas/{marker_set}/{target_fasta}/ref.fna",
    exts=multiext("resources/{species_dir}/target_fastas/{marker_set}/{target_fasta}/ref.fna", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai")
  shell:
    " cp {input.fna} {output.fna}; "
    " samtools faidx {output.fna}; "
    " bwa index {output.fna}  2> {log} "



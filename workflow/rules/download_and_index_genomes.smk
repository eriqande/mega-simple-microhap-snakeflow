


# download the requested genome and rename it as such
rule download_genome:
  params:
    url=genome_url_from_genome
  output:
    fna="resources/genomes/{genome}/{genome}.fna"
  shell:
    "echo 'wget -O {output.fna}.gz {params.url}; "
    " gunzip {output.fna}.gz' > {output.fna}"


rule bwt_index_genome:
  input:
    fna=fna_from_genome
  output:
    bwt="resources/genomes/{genome}/{genome}.fna.bwt"
  shell:
    "echo bwa index {input.fna} > {output.bwt}"

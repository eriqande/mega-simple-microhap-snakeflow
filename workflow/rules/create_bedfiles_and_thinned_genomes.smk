# make a bedfile for our target regions from the full genome
rule region_bedfiles_full:
  input:
    regfile=region_files_from_marker_set_and_genome
  output:
    bed="resources/{species_dir}/bedfiles/fullg/{marker_set}-{genome}.bed"
  log:
    "resources/{species_dir}/logs/region_bedfiles_full/{marker_set}-{genome}.log"
  conda:
    "../envs/sedgawk.yaml"
  shell:
    "cat {input.regfile} 2> {log} | "
    "gawk '{{for(i=1;i<=NF;i++) print $i}}'  2>> {log} | "
    "sed 's/[-:]/ /g'  2>> {log} | "
    "gawk 'BEGIN {{OFS=\"\t\"}} {{print $1, $2-1, $3, $1\":\"$2\"-\"$3}}' > {output.bed} 2>> {log}"


# pull the target regions out of the full genome and make a "thinned reference"
# from them.  We can see how many amplicons are getting mapped to non-target
# regions by comparing the number of reads mapped to each fragment in this thinned
# genome to the number of reads mapped to the target region when mapping to
# the whole genome.  We will put these thinned_genomes into "resources"
rule make_thinned_genomes:
  input:
    regfile=region_files_from_marker_set_and_genome,
    fa="resources/{species_dir}/genomes/{genome}/{genome}.fna"
  log:
    faidx="resources/{species_dir}/logs/make_thinned_genomes/faidx-{marker_set}-{genome}.log",
    bwaidx="resources/{species_dir}/logs/make_thinned_genomes/bwa_index-{marker_set}-{genome}.log",
    faidx_output="resources/{species_dir}/logs/make_thinned_genomes/faidx_output-{marker_set}-{genome}.log",
  envmodules:
    "aligners/bwa",
    "bio/samtools"
  conda:
    "../envs/bwasam.yaml"
  output:
    idx_files=multiext("resources/{species_dir}/thinned_genomes/{genome}/{marker_set}/thinned.fna", ".bwt", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
    fa="resources/{species_dir}/thinned_genomes/{genome}/{marker_set}/thinned.fna"
  shell:
    "samtools faidx {input.fa} $(cat {input.regfile}) > {output.fa} 2> {log.faidx}; "
    "bwa index {output.fa} 2> {log.bwaidx}; "
    "samtools faidx {output.fa} 2> {log.faidx_output}"


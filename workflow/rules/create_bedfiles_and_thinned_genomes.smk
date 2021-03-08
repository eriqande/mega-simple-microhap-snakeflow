# make a bedfile for our target regions from the full genome
rule region_bedfiles_full:
  input:
    regfile=region_files_from_marker_set_and_genome
  output:
    bed="resources/bedfiles/fullg/{marker_set}-{genome}.bed"
  shell:
    "cat {input.regfile} | "
    "awk '{{for(i=1;i<=NF;i++) print $i}}' | "
    "sed 's/[-:]/ /g' | "
    "awk 'BEGIN {{OFS=\"\t\"}} {{print $1, $2-1, $3, $1\":\"$2\"-\"$3}}' > {output.bed}"


# pull the target regions out of the full genome and make a "thinned reference"
# from them.  We can see how many amplicons are getting mapped to non-target
# regions by comparing the number of reads mapped to each fragment in this thinned
# genome to the number of reads mapped to the target region when mapping to
# the whole genome.  We will put these thinned_genomes into "resources"
rule make_thinned_genomes:
  input:
    regfile=region_files_from_marker_set_and_genome,
    fa=fna_from_genome
  output:
    fa="resources/thinned_genomes/{genome}/{marker_set}/thinned.fa",
    fai="resources/thinned_genomes/{genome}/{marker_set}/thinned.fai",
    idx="resources/thinned_genomes/{genome}/{marker_set}/thinned.fa.bwt"
  envmodules:
    "aligners/bwa",
    "bio/samtools"
  conda:
    "../envs/bwasam.yaml"
  shell:
    "echo samtools faidx {input.fa} $(cat {input.regfile}) > {output.fa}; "
    "echo samtools faidx {output.fa} > {output.fai}; "
    "echo bwa index {output.fa}; touch {output.idx}"


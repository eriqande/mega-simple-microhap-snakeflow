

# Note that we can make more rules to get variants from freebayes, if desired.

rule call_fullg_marker_sets_with_bcftools:
  input:
    bams= lambda wc: bam_tree_equivalent_files_from_marker_sets(wc, type = "fullg", trunk = "bams", ext = ".bam"), #fullg_bam_inputs_for_calling_from_marker_set_and_genome,
    fna=fna_from_genome,
    bed="resources/bedfiles/fullg/{marker_set}-{genome}.bed"
  log:
    mpileup="{run_dir}/logs/call_fullg_marker_sets_with_bcftools/bcftools_mpileup-{marker_set}-{genome}.log",
    call="{run_dir}/logs/call_fullg_marker_sets_with_bcftools/bcftools_call-{marker_set}-{genome}.log",
    sort="{run_dir}/logs/call_fullg_marker_sets_with_bcftools/bcftools_sort-{marker_set}-{genome}.log"
  conda:
    "../envs/bcftools.yaml"
  output:
    vcf="{run_dir}/vcfs/{marker_set}/fullg/{genome}/variants-bcftools.vcf"
  shell:
    "bcftools mpileup -f {input.fna} -R {input.bed} -a AD,DP,INFO/AD "
    "-B  -q 20 -Q 20 {input.bams} 2> {log.mpileup} | bcftools call -v -m  2> {log.call} | bcftools sort > {output.vcf} 2> {log.sort}"


rule call_target_fasta_marker_sets_with_bcftools:
  input:
    bams=  lambda wc: bam_tree_equivalent_files_from_marker_sets(wc, type = "target_fasta", trunk = "bams", ext = ".bam"), #target_fasta_bam_inputs_for_calling_from_marker_set_and_fasta,
    fna=fna_from_marker_set_and_target_fasta
  log:
    mpileup="{run_dir}/logs/call_target_fasta_marker_sets_with_bcftools/bcftools_mpileup-{marker_set}-{target_fasta}.log",
    call="{run_dir}/logs/call_target_fasta_marker_sets_with_bcftools/bcftools_call-{marker_set}-{target_fasta}.log",
    sort="{run_dir}/logs/call_target_fasta_marker_sets_with_bcftools/bcftools_sort-{marker_set}-{target_fasta}.log"
  conda:
    "../envs/bcftools.yaml"
  output:
    vcf="{run_dir}/vcfs/{marker_set}/target_fasta/{target_fasta}/variants-bcftools.vcf"
  shell:
    "bcftools mpileup -f {input.fna} -a AD,DP,INFO/AD "
    " -B  -q 20 -Q 20 {input.bams}  2> {log.mpileup} | bcftools call -v -m  2> {log.call} | bcftools sort > {output.vcf}  2> {log.sort}"


rule call_fullgex_remapped_markers_with_bcftools:
  input:
    bams= lambda wc: bam_tree_equivalent_files_from_marker_sets(wc, type = "fullgex_remapped", trunk = "bams", ext = ".bam"), #fullg_bam_inputs_for_calling_from_marker_set_and_genome,
    fna=thinned_fna_from_genome_and_marker_set,
  log:
    mpileup="{run_dir}/logs/call_fullgex_remapped_markers_with_bcftools/bcftools_mpileup-{marker_set}-{genome}.log",
    call="{run_dir}/logs/call_fullgex_remapped_markers_with_bcftools/bcftools_call-{marker_set}-{genome}.log",
    sort="{run_dir}/logs/call_fullgex_remapped_markers_with_bcftools/bcftools_sort-{marker_set}-{genome}.log"
  conda:
    "../envs/bcftools.yaml"
  output:
    vcf="{run_dir}/vcfs/{marker_set}/fullgex_remapped/{genome}/variants-bcftools.vcf"
  shell:
    "bcftools mpileup -f {input.fna} -a AD,DP,INFO/AD "
    "-B  -q 20 -Q 20 {input.bams} 2> {log.mpileup} | bcftools call -v -m  2> {log.call} | bcftools sort > {output.vcf} 2> {log.sort}"


# to have a look at them: bcftools query -f '%CHROM %POS %REF %ALT{0} %INFO/AD\n' variants-bcftools.vcf


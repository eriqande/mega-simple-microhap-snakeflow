

# Note that we can make more rules to get variants from freebayes, if desired.

rule call_fullg_marker_sets_with_bcftools:
  input:
    bams=fullg_bam_inputs_for_calling_from_marker_set_and_genome,
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
    bams=target_fasta_bam_inputs_for_calling_from_marker_set_and_fasta,
    fna=fna_from_marker_set_and_target_fasta
  log:
    mpileup="{run_dir}/logs/call_fullg_marker_sets_with_bcftools/bcftools_mpileup-{marker_set}-{target_fasta}.log",
    call="{run_dir}/logs/call_fullg_marker_sets_with_bcftools/bcftools_call-{marker_set}-{target_fasta}.log",
    sort="{run_dir}/logs/call_fullg_marker_sets_with_bcftools/bcftools_sort-{marker_set}-{target_fasta}.log"
  conda:
    "../envs/bcftools.yaml"
  output:
    vcf="{run_dir}/vcfs/{marker_set}/target_fasta/{target_fasta}/variants-bcftools.vcf"
  shell:
    "echo 'bcftools mpileup -f {input.fna} -a AD,DP,INFO/AD "
    "-B  -q 20 -Q 20 {input.bams}  2> {log.mpileup} | bcftools call -v -m  2> {log.call} | bcftools sort' > {output.vcf}  2> {log.sort}"



# to have a look at them: bcftools query -f '%CHROM %POS %REF %ALT{0} %INFO/AD\n' variants-bcftools.vcf


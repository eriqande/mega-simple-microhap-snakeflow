


# this calls variants from BAMS across multiple run directories.
# It just uses globbing to call from the BAMS that are already
# created for a particular marker set.
rule multi_dir_call_fullgex_remapped_markers_with_bcftools:
  input:
    bams= existing_bams_for_multirun_fullgex_remapped,
    bai= lambda wc: ["{file}.bai".format(file=f) for f in existing_bams_for_multirun_fullgex_remapped(wc)],
    fna="resources/{species_dir}/thinned_genomes/{genome}/{marker_set}/thinned.fna",
  log:
    mpileup="MULTI_RUN_RESULTS/{multi_run_dir}/{species_dir}/logs/call_fullgex_remapped_markers_with_bcftools/bcftools_mpileup-{marker_set}-{genome}.log",
    call="MULTI_RUN_RESULTS/{multi_run_dir}/{species_dir}/logs/call_fullgex_remapped_markers_with_bcftools/bcftools_call-{marker_set}-{genome}.log",
    sort="MULTI_RUN_RESULTS/{multi_run_dir}/{species_dir}/logs/call_fullgex_remapped_markers_with_bcftools/bcftools_sort-{marker_set}-{genome}.log",
    norm="MULTI_RUN_RESULTS/{multi_run_dir}/{species_dir}/logs/call_fullgex_remapped_markers_with_bcftools/bcftools_norm-{marker_set}-{genome}.log"
  threads: 6
  conda:
    "../envs/bcftools.yaml"
  output:
    vcf="MULTI_RUN_RESULTS/{multi_run_dir}/{species_dir}/vcfs/{marker_set}/fullgex_remapped/{genome}/variants-from-multi-runs-bcftools.vcf"
  shell:
    "bcftools mpileup -f {input.fna} -a AD,DP,INFO/AD "
    "-B  -q 20 -Q 20 {input.bams} 2> {log.mpileup} | bcftools call -v -m  2> {log.call} | bcftools sort 2> {log.sort} | bcftools norm -m +any -f {input.fna} > {output.vcf}  2> {log.norm} "

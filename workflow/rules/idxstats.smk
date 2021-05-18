



# Here we map the fullg-extracted reads to the thinned genomes.
rule combine_fullgex_remapped_idxstats_into_single_file:
  input:
  	idx=lambda wc: bam_tree_equivalent_files_from_marker_sets(wc, type = "fullgex_remapped", trunk = "idxstats", ext = "_idxstats.txt")
  log:
    "{run_dir}/{species_dir}/logs/combine_fullgex_remapped_idxstats_into_single_file/{marker_set}/{genome}/log.txt"
  output:
      "{run_dir}/{species_dir}/idxstats/fullgex_remapped_to_thinned/{marker_set}/{genome}/SamplesCatenated.txt",
  shell:
  	" for i in {input.idx}; do j=$(basename $i); k=${{j/_idxstats.txt/}}; "
  	" awk -v S=$k  'BEGIN {{OFS=\"\\t\"}} {{print S,$0}}'  $i;  done | "
  	" awk '$2 == \"*\" {{next}} {{print}}' > {output} 2> {log} "




# Here we map the fullg-extracted reads to the thinned genomes.
rule combine_fullg_idxstats_into_single_file:
  input:
  	idx=lambda wc: bam_tree_equivalent_files_from_marker_sets(wc, type = "fullg", trunk = "idxstats", ext = "_idxstats.txt")
  log:
    "{run_dir}/{species_dir}/logs/combine_fullg_idxstats_into_single_file/{marker_set}/{genome}/log.txt"
  output:
      "{run_dir}/{species_dir}/idxstats/fullg-extracted/{marker_set}/{genome}/SamplesCatenated.txt",
  shell:
  	" for i in {input.idx}; do j=$(basename $i); k=${{j/_idxstats.txt/}}; "
  	" awk -v S=$k  'BEGIN {{OFS=\"\\t\"}} {{print S,$0}}'  $i;  done | "
  	" awk '$2 == \"*\" {{next}} {{print}}' > {output} 2> {log} "




rule combine_target_fasta_idxstats_into_single_file:
  input:
  	idx=lambda wc: bam_tree_equivalent_files_from_marker_sets(wc, type = "target_fasta", trunk = "idxstats", ext = "_idxstats.txt")
  log:
    "{run_dir}/{species_dir}/logs/combine_fullg_idxstats_into_single_file/{marker_set}/{target_fasta}/log.txt"
  output:
      "{run_dir}/{species_dir}/idxstats/target_fastas/{marker_set}/{target_fasta}/SamplesCatenated.txt",
  shell:
  	" for i in {input.idx}; do j=$(basename $i); k=${{j/_idxstats.txt/}}; "
  	" awk -v S=$k  'BEGIN {{OFS=\"\\t\"}} {{print S,$0}}'  $i;  done | "
  	" awk '$2 == \"*\" {{next}} {{print}}' > {output} 2> {log} "


# a rule to create tabular and graphical summaries from the SamplesCatenated.txt
# I have a few new wildcards here that I should have used from the beginning:
# map_type = fullgex_remapped_to_thinned, fullg-extracted, target_fastas
# tf_or_gen = "the value that would go into the target_fasta or genome wildcards"
rule summarise_idxstats:
  input:
    "{run_dir}/{species_dir}/idxstats/{map_type}/{marker_set}/{tf_or_gen}/SamplesCatenated.txt"
  output:
    csv = "{run_dir}/{species_dir}/idxstats/{map_type}/{marker_set}/{tf_or_gen}/ordered-read-counts-table.csv",
    heatmap = "{run_dir}/{species_dir}/idxstats/{map_type}/{marker_set}/{tf_or_gen}/read-counts-heatmap.pdf",
    samp_bars = "{run_dir}/{species_dir}/idxstats/{map_type}/{marker_set}/{tf_or_gen}/read-counts-samples.pdf",
    marker_bars = "{run_dir}/{species_dir}/idxstats/{map_type}/{marker_set}/{tf_or_gen}/read-counts-markers.pdf",
  log:
    "{run_dir}/{species_dir}/logs/summarise_idxstats/{map_type}/{marker_set}/{tf_or_gen}/log.txt"
  envmodules:
    "R/4.0.3" 
  script:
    "../script/summarize-idxstats.R"

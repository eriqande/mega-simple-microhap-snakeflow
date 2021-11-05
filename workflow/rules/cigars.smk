


# count up the number of different kinds of CIGARs
# and their corresponding QUALs for each sample, chrom,
# mapping position.  This just uses gawk to do this in
# one fell swoop.  Note that I am finally dealing with
# the different mapping targets (fullg-extracted, target_fastas,
# fullgex_remapped) programmatically (with the wildcard map_mode).
# this is an aggregation step.
rule count_cigars:
	input:
		sams = lambda wc: bam_tree_equivalent_files2(wildcards=wc, trunk="sams", ext=".sam")
	output:
		cig_counts = "{run_dir}/{species_dir}/cigars/{map_mode}/{marker_set}/{genome_or_tf}/cigar_counts.tsv"
	log:
		"{run_dir}/{species_dir}/logs/count_cigars/{map_mode}---{marker_set}---{genome_or_tf}.log"
	conda:
		"../envs/gawk.yaml"
	shell:
		" gawk ' "
		"   BEGIN {{OFS=\"\t\"}} "
		"   {{samp = FILENAME; sub(/^.*\\//, \"\", samp); sub(/\\.sam/, \"\", samp);  print samp,$3,$4,$5,$6}} "
		"     '  {input.sams} 2> {log} | "
		"  gawk 'BEGIN {{SUBSEP=\"\t\"}} {{n[$0]++;}} END {{for(i in n) print i,n[i]}}' | "
		"  sort > {output.cig_counts} 2>> {log} "
		







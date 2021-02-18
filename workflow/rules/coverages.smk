
# compute the coverage at the target regions for each sample
# when mapped to the full genome
rule coverage_full:
  input:
    a="results/bedfiles/regions_full.bed",
    b="results/extracted_reads/{sample}.bam"
  output:
    "results/coverage/samples/{sample}_cov_full.tsv"
  params:
    SM="{sample}"
  envmodules:
    "bio/bedtools"
  conda:
    "../envs/bedtools.yaml"
  shell:
    "bedtools coverage -a {input.a} -b {input.b} | "
    "awk -v SM={params.SM} '{{printf(\"full\t%s\t%s\\n\", SM, $0);}}' > {output}"


# compute the coverage at the target regions for each sample
# when mapped to the full genome
rule coverage_thin:
  input:
    a="results/bedfiles/regions_thin.bed",
    b="results/mapped_reads_thin/{sample}.bam"
  output:
    "results/coverage/samples/{sample}_cov_thin.tsv"
  params:
    SM="{sample}"
  envmodules:
    "bio/bedtools"
  conda:
    "../envs/bedtools.yaml"
  shell:
    "bedtools coverage -a {input.a} -b {input.b} | "
    "awk -v SM={params.SM} '{{printf(\"thin\t%s\t%s\\n\", SM, $0);}}' > {output}"



# Now catenate all the full and thin coverages into a single
# tsv file that we can analyze easily in the tidyverse
rule catenate_coverages:
  input:
    expand("results/coverage/samples/{sample}_cov_{cond}.tsv", sample=config["samples"], cond=["full", "thin"])
  output:
    "results/coverage/all_coverages.tsv"
  shell:
    "(echo \"genome_condition\tsample\tchrom\tstart\tstop\ttarget\tnum_reads\tnum_bases\tfeature_length\tfract_bases\"; "
    "cat {input}) > {output}"



# Make a plot of the coverages
rule plot_coverages:
  input:
    "results/coverage/all_coverages.tsv"
  output:
    "results/plots/coverages-plot.pdf"
  envmodules:
    "R"
  script:
    "../script/summarize-coverages.R"


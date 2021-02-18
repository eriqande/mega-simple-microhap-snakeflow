
rule make_bam_list:
  input:
    bam=expand("results/extracted_reads/{sample}.bam", sample=config["samples"])
  output:
    "results/infiles/bamlist.list"
  shell:
    "echo {input.bam} | "
    " awk '{{for(i=1;i<=NF;i++) print $i}}' > {output} "



rule haplotype_caller:
    input:
        # single or list of bam files
        bamlist="results/infiles/bamlist.list",
        bams=expand("results/extracted_reads/{sample}.bam", sample=config["samples"]),
        ref=config["genome"],
        bed="results/bedfiles/regions_full.bed"
        # known="dbsnp.vcf"  # optional
    output:
        vcf="results/calls/everyone.vcf",
    log:
        "results/logs/gatk/haplotypecaller/everyone.log"
    params:
        extra=" --max-reads-per-alignment-start  8000 ",  # optional
        java_opts=" -Xmx4g ", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    conda:
      "../envs/gatk.yaml"
    envmodules:
      "bio/bedtools"
    shell:
      "gatk --java-options '{params.java_opts}' HaplotypeCaller {params.extra} "
      "-L {input.bed} "
      "-R {input.ref} -I {input.bamlist} "
      "-O {output.vcf} > {log} 2>&1 "

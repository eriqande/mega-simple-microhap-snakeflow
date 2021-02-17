
rule make_bam_list:
  input:
    bam=expand("extracted_reads/{sample}.bam", sample=config["samples"])
  output:
    "infiles/bamlist.list"
  shell:
    "echo {input.bam} | "
    " awk '{{for(i=1;i<=NF;i++) print $i}}' > {output} "



rule haplotype_caller:
    input:
        # single or list of bam files
        bamlist="infiles/bamlist.list",
        ref=config["genome"],
        bed="bedfiles/regions_full.bed"
        # known="dbsnp.vcf"  # optional
    output:
        vcf="calls/everyone.vcf",
    log:
        "logs/gatk/haplotypecaller/everyone.log"
    params:
        extra=" --max-reads-per-alignment-start  8000 ",  # optional
        java_opts=" -Xmx4g ", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    conda:
      "../envs/gatk.yaml"
    shell:
      "gatk --java-options '{params.java_opts}' HaplotypeCaller {params.extra} "
      "-L {input.bed} "
      "-R {input.ref} -I {input.bamlist} "
      "-O {output.vcf}"


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
        extra=" --max-reads-per-alignment-start  16000   --native-pair-hmm-threads 16 ",  # optional
        java_opts=" -Xmx64g ", # optional
    resources: cpus=16, mem_mb=75200
    conda:
      "../envs/gatk.yaml"
    shell:
      "gatk --java-options '{params.java_opts}' HaplotypeCaller {params.extra} "
      "-L {input.bed} "
      "-R {input.ref} -I {input.bamlist} "
      "-O {output.vcf} > {log} 2>&1 "
      
      
# GATK is hopeless.  Here is a better way forward:
#  bcftools mpileup -f $GENOME -R results/bedfiles/regions_full.bed -a AD,DP,INFO/AD  -b results/infiles/bamlist.list  -B  -q 20 -Q 20 | bcftools call -v -m  | bcftools sort > splud2.vcf
# 
# To have a quick look at it we can do: bcftools query -f '%CHROM  %POS  %REF  %ALT{0}   %INFO/AD\n' splud2.vcf

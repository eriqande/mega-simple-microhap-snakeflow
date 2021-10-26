


rule trimmomatic_paired_end:
    input:
        fq1=fq1_from_sample_and_run,
    	fq2=fq2_from_sample_and_run
    output:
        r1="{run_dir}/{species_dir}/trimmomatic/{sample}.1.fastq.gz",
        r2="{run_dir}/{species_dir}/trimmomatic/{sample}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="{run_dir}/{species_dir}/trimmomatic/{sample}.1.unpaired.fastq.gz",
        r2_unpaired="{run_dir}/{species_dir}/trimmomatic/{sample}.2.unpaired.fastq.gz"
    log:
    	stdout="{run_dir}/{species_dir}/logs/trimmomatic_paired_end/{sample}.stdout.log",
    	stderr="{run_dir}/{species_dir}/logs/trimmomatic_paired_end/{sample}.stderr.log"
    conda:
    	"../envs/trimmomatic.yaml"
    params:
        adapto=config["trimmo_adapter_opt"],
        opts=config["trimmo_opts"]
    shell:
    	"trimmomatic PE {input.fq1} {input.fq2} "
    	"   {output.r1} {output.r1_unpaired} "
    	"   {output.r2} {output.r2_unpaired} "
    	"   {params.adapto} {params.opts} > {log.stdout} 2> {log.stderr} "
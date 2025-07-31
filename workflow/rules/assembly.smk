"""
Assembly Module for AUTOVAR Pipeline

This module performs genome assembly:
1. SPAdes assembly
2. QUAST quality assessment
"""

# Perform de novo assembly using SPAdes
rule SPAdes:
    input:
        fastp_fq1 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R1.fq.gz",
        fastp_fq2 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R2.fq.gz"
    output:
        assembly = directory("results/Assembly/{species}/{sample}/SPAdes/"),
        link = "results/Assembly/{species}/{sample}/SPAdes/.done",
        scaffolds = protected("results/Assembly/{species}/{sample}/SPAdes/scaffolds.fasta")
    conda:
        "envs/assembly.yaml"
    params:
        threads = 12,
    log:
        "logs/Assembly/{species}/{sample}/SPAdes/SPAdes.log"
    shell:
        """
        spades.py --isolate -t {params.threads} \
        -1 {input.fastp_fq1} \
        -2 {input.fastp_fq2} \
        -o {output.assembly} \
        > {log} 2>&1 && \
        touch {output.link}
        """

# Assess assembly quality using QUAST
rule QUAST:
    input:
        scaffolds = "results/Assembly/{species}/{sample}/SPAdes/scaffolds.fasta",
        fastp_fq1 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R1.fq.gz",
        fastp_fq2 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R2.fq.gz",
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"],
        link = "results/Assembly/{species}/{sample}/SPAdes/.done",
    output:
        quast_dir = directory("results/Assembly/{species}/{sample}/QUAST/"),
        done_file = touch("results/Assembly/{species}/{sample}/QUAST/.done")
    conda:
        "envs/assembly.yaml"
    params:
        threads = 12
    log:
        "logs/Assembly/{species}/{sample}/QUAST/QUAST.log"
    shell:
        """
        mkdir -p {output.quast_dir} && \
        quast -r {input.ref} \
        -t {params.threads} \
        -1 {input.fastp_fq1} \
        -2 {input.fastp_fq2} \
        {input.scaffolds} \
        -o {output.quast_dir} \
        > {log} 2>&1
        """
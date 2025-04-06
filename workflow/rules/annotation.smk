"""
Annotation Module for AUTOVAR Pipeline

This module performs variant annotation:
1. SnpEff annotation
2. AMR detection (abricate)
3. Variant effect prediction
"""

# Download SnpEff database for species
rule snpeff_download_database:
    input:
        #none
    output:
        touch("data/reference/{species}/snpeff_database_done")
    log:    
        "logs/annotation/{species}/snpeff_database/snpeff.log"
    conda:
        "envs/annotation.yaml"
    params:
        species = lambda wildcards: references_config["reference"][wildcards.species]["name"]
    shell:
        "snpEff download {params.species} > {log} 2>&1"

# Annotate variants using SnpEff
rule snpeff_annotate:
    input:
        combined = "results/SNPIndel/{species}/Merge/variants_combined_{species}.vcf",
        database= "data/reference/{species}/snpeff_database_done"
    output:
        annotated_vcf = "results/annotation/{species}/snpeff/annotated_{species}.vcf",
        annotated_txt = "results/annotation/{species}/snpeff/annotated_{species}.txt",
        annotated_html = "results/annotation/{species}/snpeff/annotated_{species}.html"
    params:
        database_name = lambda wildcards: references_config["reference"][wildcards.species]["name"],
    log:
        "logs/annotation/{species}/snpeff_annotate/snpeff.log"
    conda:
        "envs/annotation.yaml"
    shell:
        """
        snpEff -v {params.database_name} {input.combined} > {output.annotated_vcf} 2> {output.annotated_txt} -htmlStats {output.annotated_html} 
        > {log} 2>&1
        """

# Extract fields from annotated VCF
rule snpsift_extract:
    input:
        annotated_vcf = "results/annotation/{species}/snpeff/annotated_{species}.vcf"
    output:
        tsv = "results/annotation/{species}/snpeff/extract_{species}.tsv"
    log:
        "logs/annotation/{species}/snpsift_extract/snpsift.log"
    conda:
        "envs/annotation.yaml"
    shell:
        """
        SnpSift extractFields {input.annotated_vcf} \
        CHROM POS REF ALT FILTER DP AF ANN[0].GENE ANN[0].EFFECT ANN[0].IMPACT \
        > {output.tsv} \
        2> {log}
        """

# Download Abricate databases
rule abricate_databanse_download:
    input:
        #none
    output:
        touch("results/annotation/{species}/{sample}/abricate_databases.done")
    log:    
        "logs/annotation/{species}/{sample}/abriacte_databanse_download/abricate.log"
    conda:
        "envs/annotation.yaml"
    shell:
        """
        abricate --setupdb || true \
        > {log} 2>&1
        """

# Detect AMR genes using Abricate
rule abricate:
    input:
        scaffolds = "results/Assembly/{species}/{sample}/SPAdes/scaffolds.fasta",
        link = "results/annotation/{species}/{sample}/abricate_databases.done"
    output:
        abricate = "results/annotation/{species}/{sample}/abricate/abricate.tab"
    log:
        "logs/annotation/{species}/{sample}/abricate/abricate.log"
    conda:
        "envs/annotation.yaml"
    shell:
        """
        abricate --db ncbi \
        {input.scaffolds} > {output.abricate} \
        2> {log}
        """
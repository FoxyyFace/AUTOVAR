"""
Quality Control Module for AUTOVAR Pipeline

This module performs quality control on sequencing data:
1. Raw data QC (FastQC)
2. Adapter trimming (trim_galore)
3. Quality filtering (FastP)
4. Final QC reports (MultiQC)
"""

# This Snakemake workflow performs quality control (QC) on sequencing data.
# It includes steps for raw data QC, adapter trimming, quality filtering, and final QC aggregation.
# Tools used include FastQC, MultiQC, and FastP, along with custom scripts for adapter trimming.

# Rule to perform quality control on raw sequencing data using FastQC.
rule rawdata_fastqc:
    input:
        fq1 = get_raw_fq1,
        fq2 = get_raw_fq2
    output:
        directory("results/QC/{species}/{sample}/fastqc/raw/")
    params:
        threads = 8,
        memory = 10000
    conda:
        "envs/qc_check.yaml"
    log:
        "logs/QC/{species}/{sample}/rawdata_fastqc/fastqc.log"
    shell:
        """
        mkdir -p {output} && 
        fastqc -t {params.threads} --memory {params.memory} \
        {input.fq1} {input.fq2} \
        -o {output} \
        > {log} 2>&1
        """

# Aggregate FastQC reports for raw data
rule raw_multiqc:
    input:
        lambda wildcards: expand(
            "results/QC/{species}/{sample}/fastqc/raw/",
            species=wildcards.species,
            sample=samples_config['samples'][wildcards.species].keys()
        )
    output:
        directory("results/QC/{species}/multiqc/raw/"),
        touch("results/QC/{species}/multiqc/raw/.done")
    conda:
        "envs/qc_check.yaml"
    log:
        "logs/QC/{species}/raw_multiqc/multiqc.log"
    shell:
        """
        mkdir -p {output[0]} && 
        multiqc {input} -o {output[0]} > {log} 2>&1
        """

# Trim adapter sequences using trim_galore
# Default max cycle is 15.It will stop when adapter is cleared or max cycle is reached.
# There are two modes: one is base on absolute count, the other is based on percentage.You can control it in the "config.yaml"
rule trim_adapter:
    input:
        fq1 = get_raw_fq1, 
        fq2 = get_raw_fq2
    output:
        output_dir = directory("results/QC/{species}/{sample}/trim_galore/"),
        fq1 = "results/QC/{species}/{sample}/trim_galore/{sample}_1_val_1.fq.gz",
        fq2 = "results/QC/{species}/{sample}/trim_galore/{sample}_2_val_2.fq.gz"
    params:
        threads = 8,
        quality = 20,
        max_n = 3,
        strict_qc = lambda wildcards: config["global"].get("strict_qc", False)
    log:
        "logs/QC/{species}/{sample}/trim_adapter/trim_galore.log"
    conda:
        "envs/qc_trim.yaml"
    script:
        "scripts/trim_galore_cycle.py"

# QC on trimmed data
rule trim_adapter_fastqc:
    input:
        trimmed_fq1 = "results/QC/{species}/{sample}/trim_galore/{sample}_1_val_1.fq.gz",
        trimmed_fq2 = "results/QC/{species}/{sample}/trim_galore/{sample}_2_val_2.fq.gz"
    output:
        directory("results/QC/{species}/{sample}/fastqc/trimmed_trim/"),
        r1_zip = "results/QC/{species}/{sample}/fastqc/trimmed_trim/{sample}_1_val_1_fastqc.zip",
        r2_zip = "results/QC/{species}/{sample}/fastqc/trimmed_trim/{sample}_2_val_2_fastqc.zip"
    conda:
        "envs/qc_check.yaml"
    params:
        threads = 8,
        memory = "10000"
    log:
        "logs/QC/{species}/{sample}/trim_adapter_fastqc/fastqc.log"
    shell:
        """
        mkdir -p {output[0]} && 
        fastqc -t {params.threads} --memory {params.memory} \
        {input.trimmed_fq1} {input.trimmed_fq2} \
        -o {output[0]} \
        > {log} 2>&1
        """

# Generate FastP parameters
# For setting the "-t/T and -f/F" parmaters.
rule fastp_params:
    input:
        r1 = "results/QC/{species}/{sample}/fastqc/trimmed_trim/{sample}_1_val_1_fastqc.zip",
        r2 = "results/QC/{species}/{sample}/fastqc/trimmed_trim/{sample}_2_val_2_fastqc.zip"
    output:
        "results/QC/{species}/{sample}/fastp_params/{sample}.json"
    log:
        "logs/QC/{species}/{sample}/fastp_params/fastp_params.log"
    script:
        "scripts/fastp_params.py"

# Quality filtering using FastP
rule fastp:
    input:
        trimmed_fq1 = "results/QC/{species}/{sample}/trim_galore/{sample}_1_val_1.fq.gz",
        trimmed_fq2 = "results/QC/{species}/{sample}/trim_galore/{sample}_2_val_2.fq.gz",
        params_file = "results/QC/{species}/{sample}/fastp_params/{sample}.json"
    output:
        directory("results/QC/{species}/{sample}/fastp/"),
        html = "results/QC/{species}/{sample}/fastp/report.html",
        json = "results/QC/{species}/{sample}/fastp/metrics.json",
        fastp_fq1 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R1.fq.gz",
        fastp_fq2 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R2.fq.gz"
    params:
        f_value = lambda wildcards, input: json.load(open(input.params_file))["f_value"],
        F_value = lambda wildcards, input: json.load(open(input.params_file))["F_value"],
        t_value = lambda wildcards, input: json.load(open(input.params_file))["t_value"],
        T_value = lambda wildcards, input: json.load(open(input.params_file))["T_value"]
    log:
        "logs/QC/{species}/{sample}/fastp/fastp.log"
    conda:
        "envs/qc_trim.yaml" 
    shell: 
        """
        fastp \
        -i {input.trimmed_fq1}  -I {input.trimmed_fq2}  \
        -o {output.fastp_fq1}  -O {output.fastp_fq2}  \
        -f {params.f_value} -F {params.F_value} \
        -t {params.t_value} -T {params.T_value} \
        -q 20 -l 50 \
        --json {output.json}  \
        --html {output.html}  \
        > {log} 2>&1 
        """

# QC on FastP-processed data
rule fastp_fastqc:
    input:
        fastp_fq1 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R1.fq.gz",
        fastp_fq2 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R2.fq.gz"
    output:
        directory("results/QC/{species}/{sample}/fastqc/final/")
    conda:
        "envs/qc_check.yaml"
    params:
        threads = 8,
        memory = "10000"
    log:
        "logs/QC/{species}/{sample}/fastp_fastqc/fastqc.log"
    shell:
        """
        mkdir -p {output} && 
        fastqc -t {params.threads} --memory {params.memory} \
        {input.fastp_fq1} {input.fastp_fq2} \
        -o {output} \
        > {log} 2>&1
        """

# Aggregate all QC reports
rule final_multiqc:
    input: 
        raw_fastqc = lambda wildcards: expand(
            "results/QC/{species}/{sample}/fastqc/raw/",
            species=wildcards.species,
            sample=samples_config['samples'][wildcards.species].keys()
        ),
        trim_fastqc = lambda wildcards: expand(
            "results/QC/{species}/{sample}/fastqc/trimmed_trim/", 
            species=wildcards.species,
            sample=samples_config['samples'][wildcards.species].keys()
        ),
        fastp_fastqc = lambda wildcards: expand(
            "results/QC/{species}/{sample}/fastqc/final/", 
            species=wildcards.species,
            sample=samples_config['samples'][wildcards.species].keys()
        ),
        fastp_metrics = lambda wildcards: expand(
            "results/QC/{species}/{sample}/fastp/metrics.json", 
            species=wildcards.species,
            sample=samples_config['samples'][wildcards.species].keys()
        ),
        link = "results/QC/{species}/multiqc/raw/.done"
    output: 
        directory("results/QC/{species}/multiqc/final/")
    log:
        "logs/QC/{species}/final_multiqc/multiqc.log"   
    conda:
        "envs/qc_check.yaml" 
    shell:
        "multiqc {input.raw_fastqc} {input.trim_fastqc} {input.fastp_fastqc} {input.fastp_metrics} -o {output} > {log} 2>&1"

"""
Mapping Module for AUTOVAR Pipeline

This module performs read mapping:
1. BWA mapping
2. SAM to BAM conversion
3. Sort and index BAM files
4. Mapping statistics
"""

# Index reference genome for BWA
rule bwa_index:
    input:
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"]
    output:
        touch("data/reference/{species}/bwa_index.done")
    log:
        "logs/Mapping/{species}/bwa_index/bwa.log"
    conda:
        "envs/mapping.yaml"
    shell:
        """
        bwa index {input.ref} > {log} 2>&1 && \
        mkdir -p data/reference/{wildcards.species} && \
        touch data/reference/{wildcards.species}/bwa_index.done
        """

# Map reads to reference using BWA-MEM
rule bwa_mem:
    input:
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"],
        index_done = "data/reference/{species}/bwa_index.done",
        fastp_fq1 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R1.fq.gz",
        fastp_fq2 = "results/QC/{species}/{sample}/fastp/{sample}_clean_R2.fq.gz"
    output:
        bam = temp("results/Mapping/{species}/{sample}/bwa/bwa_mem_{sample}.bam")
    log:
        "logs/Mapping/{species}/{sample}/bwa_mem/bwa.log"
    params:
        threads=8
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa mem -t {params.threads} -M -a -T 30 -B 4 {input.ref} {input.fastp_fq1} {input.fastp_fq2} | "
        "samtools view -bS -@ {params.threads} -o {output.bam} > {log} 2>&1"

# Fix mate information in BAM file
rule samtools_fixmate:
    input:
        bam = "results/Mapping/{species}/{sample}/bwa/bwa_mem_{sample}.bam"
    output:
        fixmate_bam = "results/Mapping/{species}/{sample}/samtools/fixmate_{sample}.bam"
    log:
        "logs/Mapping/{species}/{sample}/samtools_fixmate/samtools.log"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools fixmate -m {input.bam} {output.fixmate_bam} > {log} 2>&1"

# Sort BAM file by coordinates
rule samtools_sort:
    input:
        fixmate_bam = "results/Mapping/{species}/{sample}/samtools/fixmate_{sample}.bam"
    output:
        sorted_bam = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam"
    log:
        "logs/Mapping/{species}/{sample}/samtools_sorted/samtools.log"
    params:
        threads = 8
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools sort -@ {params.threads} -o {output.sorted_bam} {input.fixmate_bam} > {log} 2>&1"

# Index sorted BAM file
rule samtools_index:
    input:
        sorted_bam = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam"
    output:
        bai = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam.bai"
    log:
        "logs/Mapping/{species}/{sample}/samtools_index/samtools.log"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools index {input.sorted_bam} > {log} 2>&1"

# Generate mapping statistics
rule samtools_stats:
    input:
        sorted_bam = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam",
        bai = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam.bai"
    output:
        stats = "results/Mapping/{species}/{sample}/samtools/{sample}.stats"
    log:
        "logs/Mapping/{species}/{sample}/samtools_stats/samtools.log"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools stats {input.sorted_bam} > {output.stats} 2> {log}"

# Calculate read depth
rule samtools_depth:
    input:
        sorted_bam = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam",
        bai = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam.bai"
    output:
        depth = "results/Mapping/{species}/{sample}/samtools/{sample}.depth"
    log:
        "logs/Mapping/{species}/{sample}/samtools_depth/samtools.log"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools depth -q 20 -Q 30 {input.sorted_bam} > {output.depth} 2> {log}"






"""
Duplicate Removal Module for AUTOVAR Pipeline

This module removes PCR duplicates:
1. Mark duplicates using Picard
2. Remove duplicates
3. Index final BAM files
"""

# Add read groups to BAM file
rule gatk_add_read_groups:
    input:
        bam = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam"
    output:
        bam = "results/mk_rm_dup/{species}/{sample}/gatk/sorted_rg_{sample}.bam"
    log:
        "logs/mk_rm_dup/{species}/{sample}/gatk_add_read_groups/gatk.log"
    conda:
        "envs/variantcalling.yaml"
    params:
        RGID = lambda wildcards: samples_config["samples"][wildcards.species][wildcards.sample]["ReadGroup"]["RGID"],
        RGLB = lambda wildcards: samples_config["samples"][wildcards.species][wildcards.sample]["ReadGroup"]["RGLB"],
        RGPL = lambda wildcards: samples_config["samples"][wildcards.species][wildcards.sample]["ReadGroup"]["RGPL"],
        RGPU = lambda wildcards: samples_config["samples"][wildcards.species][wildcards.sample]["ReadGroup"]["RGPU"],
        RGSM = lambda wildcards: samples_config["samples"][wildcards.species][wildcards.sample]["ReadGroup"]["RGSM"]
    shell:
        """
        gatk AddOrReplaceReadGroups \
        -I {input.bam} -O {output.bam} \
        -ID {params.RGID} \
        -LB {params.RGLB} \
        -PL {params.RGPL} \
        -PU {params.RGPU} \
        -SM {params.RGSM} \
        > {log} 2>&1
        """

# Index BAM file with read groups
rule add_read_groups_index:
    input:
        bam = rules.gatk_add_read_groups.output.bam
    output:
        bai = "results/mk_rm_dup/{species}/{sample}/gatk/sorted_rg_{sample}.bam.bai"
    log:
        "logs/mk_rm_dup/{species}/{sample}/add_read_groups_index/samtools.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        "samtools index {input.bam} > {log} 2>&1"

# Mark and remove PCR duplicates
rule gatk_MarkDuplicates:
    input:
        bam = rules.gatk_add_read_groups.output.bam,
        bai = rules.add_read_groups_index.output.bai 
    output:
        bam = "results/mk_rm_dup/{species}/{sample}/gatk/rmdup_{sample}.bam",
        metrics = "results/mk_rm_dup/{species}/{sample}/gatk/rmdup_metrics_{sample}.txt"
    log:
        "logs/mk_rm_dup/{species}/{sample}/gatk_MarkDuplicates/gatk.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk MarkDuplicates \
        --REMOVE_DUPLICATES true \
        -MAX_FILE_HANDLES 1024 \
        -I {input.bam} \
        -O {output.bam} \
        -M {output.metrics} \
        >> {log} 2>&1
        """

# Index deduplicated BAM file
rule index_rmbam:
    input:
        bam = "results/mk_rm_dup/{species}/{sample}/gatk/rmdup_{sample}.bam"
    output:
        bai = "results/mk_rm_dup/{species}/{sample}/gatk/rmdup_{sample}.bam.bai"
    log:
        "logs/mk_rm_dup/{species}/{sample}/index_rmbam/samtools.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        "samtools index {input.bam} > {log} 2>&1"

# Create reference index and dictionary
rule fasidx_dict_ref:
    input:
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"]
    output:
        faidx = "data/reference/{species}/{species}.fa.fai",
        dictfa = "data/reference/{species}/{species}.dict"
    log:
        "logs/mk_rm_dup/{species}/fasidx_dict_ref/samtools_gatk.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        mkdir -p data/reference/{wildcards.species} && \
        samtools faidx {input.ref} && \
        gatk CreateSequenceDictionary -R {input.ref} \
        >> {log} 2>&1
        """
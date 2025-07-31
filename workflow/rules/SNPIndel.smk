"""
Variant Calling Module for AUTOVAR Pipeline

This module performs variant calling:
1. GATK HaplotypeCaller
2. Variant filtering
3. Variant statistics
"""

# Call variants using GATK HaplotypeCaller
rule gatk_haplotypecaller:
    input:
        bam = "results/mk_rm_dup/{species}/{sample}/gatk/rmdup_{sample}.bam",
        bai = "results/mk_rm_dup/{species}/{sample}/gatk/rmdup_{sample}.bam.bai",
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"]
    output:
        gvcf = "results/SNPIndel/{species}/{sample}/gatk/haplotypecaller_{sample}.g.vcf"
    log:
        "logs/SNPIndel/{species}/varcall/gatk_haplotypecaller/{sample}.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk HaplotypeCaller \
        -ERC GVCF \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        > {log} 2>&1
        """

# Create intervals for parallel processing
rule create_intervals:
    input:
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"],
        fai = lambda wildcards: references_config["reference"][wildcards.species]["path"] + ".fai"
    output:
        intervals = "data/reference/{species}/{species}.interval_list"  
    log:
        "logs/SNPIndel/{species}/create_intervals/gatk.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk ScatterIntervalsByNs \
            -R {input.ref} \
            -O {output.intervals} \
            -OT BOTH \
            > {log} 2>&1
        """

# Import GVCFs into GenomicsDB
rule genomicsdb_import:
    input:
        gvcfs = lambda wildcards: expand(
            "results/SNPIndel/{species}/{sample}/gatk/haplotypecaller_{sample}.g.vcf",
            species=wildcards.species,
            sample=samples_config["samples"][wildcards.species].keys()
        ),
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"],
        intervals = "data/reference/{species}/{species}.interval_list"
    output:
        db = directory("results/SNPIndel/{species}/gatk/db_{species}")
    params:
        vcf_args = lambda wildcards, input: " ".join([f"-V {vcf}" for vcf in input.gvcfs]),
        batch_size = 50,
        reader_threads = 4
    log:
        "logs/SNPIndel/{species}/genomicsdb_import/gatk.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk GenomicsDBImport \
            -R {input.ref} \
            --genomicsdb-workspace-path {output.db} \
            -L {input.intervals} \
            {params.vcf_args} \
            --batch-size {params.batch_size} \
            --reader-threads {params.reader_threads} \
            --merge-input-intervals \
            > {log} 2>&1
        """

# Genotype GVCFs to produce final VCF
rule genotype_gvcfs:
    input:
        db = "results/SNPIndel/{species}/gatk/db_{species}",
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"]
    output:
        vcf = "results/SNPIndel/{species}/gatk/combined_{species}.vcf"
    params:
        include_non_variant_sites = True,
        sites_only = False
    log:
        "logs/SNPIndel/{species}/genotype_gvcfs/gatk.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            --include-non-variant-sites {params.include_non_variant_sites} \
            > {log} 2>&1
        """

# Select SNPs and INDELs from combined VCF
rule gatk_selectvariants:
    input:
        vcf = "results/SNPIndel/{species}/gatk/combined_{species}.vcf"
    output:
        snp_vcf = "results/SNPIndel/{species}/gatk/SNP_{species}.vcf", 
        indel_vcf = "results/SNPIndel/{species}/gatk/Indel_{species}.vcf"
    params:
        snp_filter = "--select-type SNP",
        indel_filter = "--select-type INDEL"
    log:
        "logs/SNPIndel/{species}/gatk_selectvariants/gatk.log" 
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk SelectVariants \
            -V {input.vcf} \
            -O {output.snp_vcf} \
            {params.snp_filter} && \
        gatk SelectVariants \
            -V {input.vcf} \
            -O {output.indel_vcf} \
            {params.indel_filter} \
            > {log} 2>&1
        """

# Apply variant filtering criteria
rule gatk_mark_variantfiltration:
    input:
        snp_vcf = "results/SNPIndel/{species}/gatk/SNP_{species}.vcf", 
        indel_vcf = "results/SNPIndel/{species}/gatk/Indel_{species}.vcf",
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"]
    output:
        snp_mark_vcf = "results/SNPIndel/{species}/gatk/marked_SNP_{species}.vcf",
        indel_mark_vcf = "results/SNPIndel/{species}/gatk/marked_Indel_{species}.vcf"
    params:
        SNP_PARAM = lambda wildcards: config["gatk_VariantFiltration"]["SNP"]["filter_expression"],
        INDEL_PARAM = lambda wildcards: config["gatk_VariantFiltration"]["INDEL"]["filter_expression"],
        SNP_NAME = lambda wildcards: config["gatk_VariantFiltration"]["SNP"]["filter_name"],
        INDEL_NAME = lambda wildcards: config["gatk_VariantFiltration"]["INDEL"]["filter_name"]
    log:
        "logs/SNPIndel/{species}/gatk_mark_variantfiltration/gatk.log" 
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk VariantFiltration \
            -V {input.snp_vcf} \
            -O {output.snp_mark_vcf} \
            --filter-expression "{params.SNP_PARAM}" \
            --filter-name "{params.SNP_NAME}" && \
        gatk VariantFiltration \
            -V {input.indel_vcf} \
            -O {output.indel_mark_vcf} \
            --filter-expression "{params.INDEL_PARAM}" \
            --filter-name "{params.INDEL_NAME}" \
            > {log} 2>&1
        """

# Select filtered variants
rule gatk_selectvariants_filtermark:
    input:
        snp_mark_vcf = "results/SNPIndel/{species}/gatk/marked_SNP_{species}.vcf",
        indel_mark_vcf = "results/SNPIndel/{species}/gatk/marked_Indel_{species}.vcf"
    output:
        snp_filtered_vcf = "results/SNPIndel/{species}/SNP/filtered_SNP_{species}.vcf", 
        indel_filtered_vcf = "results/SNPIndel/{species}/Indel/filtered_Indel_{species}.vcf"
    log:
        "logs/SNPIndel/{species}/gatk_selectvariants_filtermark/gatk.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk SelectVariants \
            -V {input.snp_mark_vcf} \
            -O {output.snp_filtered_vcf} \
            --exclude-filtered && \
        gatk SelectVariants \
            -V {input.indel_mark_vcf} \
            -O {output.indel_filtered_vcf} \
            --exclude-filtered \
            > {log} 2>&1
        """

# Merge SNP and INDEL VCFs
rule merge_vcfs:
    input:
        snp_filtered_vcf = "results/SNPIndel/{species}/SNP/filtered_SNP_{species}.vcf", 
        indel_filtered_vcf = "results/SNPIndel/{species}/Indel/filtered_Indel_{species}.vcf"
    output:
        combined = "results/SNPIndel/{species}/Merge/variants_combined_{species}.vcf"
    log:
        "logs/SNPIndel/{species}/merge_vcfs/gatk.log"
    conda:
        "envs/variantcalling.yaml"
    shell:
        """
        gatk MergeVcfs \
            -I {input.snp_filtered_vcf} \
            -I {input.indel_filtered_vcf} \
            -O {output.combined} \
            > {log} 2>&1
        """



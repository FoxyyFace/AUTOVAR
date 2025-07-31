rule bowtie2_index:
    input:
        ref = lambda wildcards: references_config["reference"][wildcards.species]["path"]
    output:
        touch("data/reference/{species}/bowtie2_index.done"),
    threads: 8
    params:
        prefix = "data/reference/{species}/{species}",
    conda:
        "conda_env/rule_rna_single/map"
    log:
        "logs/Mapping/{species}/bowtie2_index.log"
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.ref} \
            {params.prefix} \
            > {log} 2>&1
        """

rule bowtie2_map:
    input:
        fq = rules.cutadapt.output.trimmed_fq,
        link = "data/reference/{species}/bowtie2_index.done"
    output:
        sam = "results/Mapping/{species}/{sample}/bowtie2/bowtie2_{sample}.sam",
        unpaired = "results/Mapping/{species}/{sample}/bowtie2/unpaired_{sample}.fastq.gz"
    threads: 8
    params:
        prefix = "data/reference/{species}/{species}"
    conda:
        "conda_env/rule_rna_single/map"
    log:
        "logs/Mapping/{species}/{sample}/bowtie2_map.log"
    shell:
        """
        bowtie2 \
            --threads {threads} \
            -x {params.prefix} \
            -U {input.fq} \
            -S {output.sam} \
            --un-gz {output.unpaired} \
            > {log} 2>&1
        """

rule samtools_view:
    input:
        sam = rules.bowtie2_map.output.sam,
    output:
        bam = temp("results/Mapping/{species}/{sample}/samtools/bam_{sample}.bam")
    threads: 8
    params:
        prefix = "data/reference/{species}/{species}"
    conda:
        "conda_env/rule_rna_single/map"
    log:
        "logs/Mapping/{species}/{sample}/samtools_view.log"
    shell:
        """
        samtools view \
            -@ {threads} \
            -bS {input.sam} \
            -o {output.bam} \
            > {log} 2>&1
        """

rule samtools_sort:
    input:
        bam = rules.samtools_view.output.bam,
    output:
        sorted_bam = "results/Mapping/{species}/{sample}/samtools/sorted_{sample}.bam"
    threads: 8
    params:
        prefix = "data/reference/{species}/{species}"
    conda:
        "conda_env/rule_rna_single/map"
    log:
        "logs/Mapping/{species}/{sample}/samtools_sort.log"
    shell:
        """
        samtools sort \
            {input.bam} \
            -@ {threads} \
            -o {output.sorted_bam} \
            > {log} 2>&1
        """


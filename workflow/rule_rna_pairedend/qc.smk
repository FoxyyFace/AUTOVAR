rule rawdata_fastqc:
    input:
        fq1 = get_raw_fq1,
        fq2 = get_raw_fq2
    output:
        report_dir = directory("results/QC/{species}/{sample}/fastqc/raw"),
        done_flag = "results/QC/{species}/{sample}/.status/raw_fastqc.done"
    threads: 8
    params:
        memory = 10000
    conda:
        "conda_env/rule_rna_pairedend/qc"
    log:
        "logs/QC/{species}/{sample}/raw_fastqc.log"
    shell:
        """
        mkdir -p {output.report_dir} && \
        fastqc -t {threads} \
               --memory {params.memory} \
               {input.fq1} {input.fq2}\
               -o {output.report_dir} \
               > {log} 2>&1
        touch {output.done_flag}
        """

rule cutadapt:
    input:
        fq1 = get_raw_fq1,
        fq2 = get_raw_fq2,
        qc_done = rules.rawdata_fastqc.output.done_flag,
    output:
        trimmed_fq1 = "results/QC/{species}/{sample}/trimmed/{sample}_1_val_1.fq.gz",
        trimmed_fq2 = "results/QC/{species}/{sample}/trimmed/{sample}_2_val_2.fq.gz",
        trimmed_dir = directory("results/QC/{species}/{sample}/trimmed/"),
    threads: 8
    params:
        adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        min_length = 22,
    conda:
        "conda_env/rule_rna_pairedend/qc"
    log:
        "logs/QC/{species}/{sample}/cutadapt.log"
    shell:
        """
        trim_galore \
            --paired \
            -j {threads} \
            {input.fq1} {input.fq2} \
            -o {output.trimmed_dir} \
            > {log} 2>&1
        """
        # """
        # cutadapt -j {threads} \
        #         -a {params.adapter} \
        #         -m {params.min_length} \
        #         -o {output.trimmed_fq} \
        #         {input.fastq} \
        #         > {log} 2>&1
        # """

rule fastqc_trimmed:
    input:
        fq1 = rules.cutadapt.output.trimmed_fq1,
        fq2 = rules.cutadapt.output.trimmed_fq2,
    output:
        report_dir = directory("results/QC/{species}/{sample}/fastqc/trimmed"),
        done_flag = "results/QC/{species}/{sample}/.status/trimmed_fastqc.done"
    threads: 8
    params:
        memory = 10000
    conda:
        "conda_env/rule_rna_pairedend/qc"
    log:
        "logs/QC/{species}/{sample}/trimmed_fastqc.log"
    shell:
        """
        mkdir -p {output.report_dir} && \
        fastqc -t {threads} \
               --memory {params.memory} \
               {input.fq1} {input.fq2} \
               -o {output.report_dir} \
               > {log} 2>&1
        touch {output.done_flag}
        """

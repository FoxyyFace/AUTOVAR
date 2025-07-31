rule rawdata_fastqc:
    input:
        fq = get_raw_fq,
    output:
        report_dir = directory("results/QC/{species}/{sample}/fastqc/raw"),
        done_flag = "results/QC/{species}/{sample}/.status/raw_fastqc.done"
    threads: 8
    params:
        memory = 10000
    conda:
        "conda_env/rule_rna_single/qc"
    log:
        "logs/QC/{species}/{sample}/raw_fastqc.log"
    shell:
        """
        mkdir -p {output.report_dir} && \
        fastqc -t {threads} \
               --memory {params.memory} \
               {input.fq} \
               -o {output.report_dir} \
               > {log} 2>&1
        touch {output.done_flag}
        """

rule cutadapt:
    input:
        fastq = get_raw_fq,
        qc_done = rules.rawdata_fastqc.output.done_flag,
    output:
        trimmed_fq = "results/QC/{species}/{sample}/trimmed/{sample}_trimmed.fq.gz",
        trimmed_dir = directory("results/QC/{species}/{sample}/trimmed/")
    threads: 8
    params:
        adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        min_length = 22,
    conda:
        "conda_env/rule_rna_single/qc"
    log:
        "logs/QC/{species}/{sample}/cutadapt.log"
    shell:
        """
        trim_galore \
            -j {threads} \
            {input.fastq} \
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
        fq = rules.cutadapt.output.trimmed_fq,
    output:
        report_dir = directory("results/QC/{species}/{sample}/fastqc/trimmed"),
    threads: 8
    params:
        memory = 10000
    conda:
        "conda_env/rule_rna_single/qc"
    log:
        "logs/QC/{species}/{sample}/trimmed_fastqc.log"
    shell:
        """
        mkdir -p {output.report_dir} && \
        fastqc -t {threads} \
               --memory {params.memory} \
               {input.fq} \
               -o {output.report_dir} \
               > {log} 2>&1
        """
rule FeatureCount:
    input:
        bam = rules.samtools_sort.output.sorted_bam,
        gff = "data/reference/{species}/{species}.gff"
        # "data/reference/{species}/sicX.gtf",
    output:
        counts = "results/Analyse/{species}/{sample}/featureCounts/counts_{sample}.txt"
    threads: 8
    log:
        "logs/Analyse/{species}/{sample}/featureCounts.log"
    conda:
        "conda_env/rule_rna_single/analyse"
    shell:
        """
        featureCounts \
            -F "GFF" \
            -a {input.gff} \
            -o {output.counts} \
            -g "locus_tag" \
            -t "CDS,sRNA" \
            -s 1 \
            -O \
            -T {threads} \
            {input.bam}
        """

# rule kraken2:
#     input:
#         fq = rules.bowtie2_map.output.unpaired,
#         db = lambda path: references_config["kraken"]["standard_16"]
#     output:
#         report = "results/kraken2/{species}/{sample}/kraken2_report.txt",
#         kraken = "results/kraken2/{species}/{sample}/kraken2_output.txt"
#     threads: 8
#     log:
#         "logs/kraken2/{species}/{sample}/kraken2.log"
#     conda:
#         "envs/analyse.yaml"
#     shell:
#         """
#         kraken2 \
#             --db {input.db} \
#             --threads {threads} \
#             --report {output.report} \
#             --output {output.kraken} \
#             --use-name \
#             {input.fq} \
#             > {log} 2>&1
#         """

# SPECIES, SAMPLES = glob_wildcards("results/kraken2/{species}/{sample}/kraken2_report.txt")

# rule aggregate_kraken_reports:
#     input:
#         expand("results/kraken2/{species}/{sample}/kraken2_report.txt", 
#                species=SPECIES, sample=SAMPLES)
#     output:
#         summary_table = "results/kraken2/summary/top10_species.tsv",
#         histogram = "results/kraken2/summary/top10_species_histogram.png"
#     params:
#         output_dir = directory("results/kraken2/summary"),
#         sample_info = samples_config["samples"]
#     conda:
#         "envs/analyse.yaml"
#     script:
#         "scripts/summarize_kraken_reports.py"

# rule extract_gene_expression:
#     input:
#         counts = "results/Analyse/{species}/{sample}/featureCounts/counts_{sample}.txt"
#     output:
#         expr = "results/Analyse/{species}/expression/{sample}/PA14_46160_expression.tsv"
#     params:
#         gene_id = "PA14_46160"  # 目标基因ID
#     conda:
#         "envs/analyse.yaml"
#     script:
#         "scripts/extract_gene_expression.py"

# rule aggregate_gene_expression:
#     input:
#         expr_files = expand(
#             "results/Analyse/{species}/expression/{sample}/PA14_46160_expression.tsv",
#             species=SPECIES,
#             sample=SAMPLES
#         ),
#         config = "configure/samples_config.yaml"  # 添加配置文件输入
#     output:
#         table = "results/Analyse/summary/PA14_46160_expression_table.tsv",
#         plot = "results/Analyse/summary/PA14_46160_volcano_plot.png"
#     conda:
#         "envs/analyse.yaml"
#     script:
#         "scripts/aggregate_gene_expression.py"
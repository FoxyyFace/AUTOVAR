"""
Visualization Module for AUTOVAR Pipeline

This module generates visualizations:
1. Depth coverage plots
2. Variant effect plots
3. AMR variant plots
"""

# Generate depth coverage plots
rule plots_depth:
    input:
        depth = "results/Mapping/{species}/{sample}/samtools/{sample}.depth"
    output:
        plot = "results/plots/{species}/{sample}/depth/{sample}_depth.svg"
    log:
        "logs/plots/{species}/plot_depth/{sample}.log"
    conda:
        "envs/plots.yaml"
    script:
        "scripts/depth_coverage_plots.py"

# Generate SnpEff effect plots
rule snpeff_plots:
    input:
        "results/annotation/{species}/snpeff/annotated_{species}.vcf"
    output:
        report = "results/plots/{species}/snpeff_report_{species}.txt",
        plot = "results/plots/{species}/snpeff_effects_{species}.png"
    conda:
        "envs/plots.yaml"
    log:
        "logs/plots/{species}/snpeff_plots/{species}.log"
    script:
        "scripts/snpeff_plots.py"

# Get Abricate files for all samples
def get_abricate_files(wildcards):
    species = wildcards.species
    if species not in samples_config["samples"]:
        raise ValueError(f"No samples found for species: {species}")
    files = []
    for sample in samples_config["samples"][species]:
        files.append(f"results/annotation/{species}/{sample}/abricate/abricate.tab")
    return files

# Integrate AMR variant information
rule integrate_amr:
    input:
        variants = "results/annotation/{species}/snpeff/annotated_{species}.vcf",
        abricate_files = get_abricate_files 
    output:
        "results/plots/{species}/annotation/{species}_amr_variants.tsv"
    conda:
        "envs/plots.yaml" 
    log:
        "logs/plots/{species}/integrate_amr/{species}.log"
    script:
        "scripts/map_abricate.py"


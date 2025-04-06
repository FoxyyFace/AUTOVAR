# AUTOVAR - Automated Bacterial Variant Analysis Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.29.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3-3776AB.svg?logo=python&logoColor=white)](https://www.python.org)

## Overview
AUTOVAR is a Snakemake pipeline designed for bacterial genomic variant analysis. It supports multi-species and multi-sample analysis.Paired-end data of NGS of bacteria have been tested currently. The pipeline provides comprehensive analysis including:
- Quality Control (QC)
   - Cycle removal of adapters, up to a maximum of 15 cycles by default
   - Automatically set fastp cutting parameters
- Read Mapping
- Genome Assembly
- Variant Calling (SNP&Indel)
- Variant Annotation (AMR detection using abricate)
- Simple Analyse
   - Depth distribution diagram of the bam file
   - A simple report diagram of the SnpEff annotation file
   - Detect if the variant occurs in the ARM region of the genome

Version: 0.8.0

## Quick Start Guide

### 1. Environment Setup(reference from Snakemake)
```bash
# Change to the project directory

cd AUTOVAR

# Create conda environment using mamba (recommended for faster installation)
# (conda install -n base -c conda-forge mamba)if you don't install mamba

mamba env create --name AUTOVAR --file environment.yaml

# Activate the environment

conda activate AUTOVAR
```

### 2. Configuration
Before running the pipeline, you need to configure the following files:
- Edit `configure/config.yaml` for general settings
- Edit `configure/samples_config.yaml` for sample information
- Edit `configure/references_config.yaml` for reference genome paths
DON'T FORGET:
## .fa
- Place the reference genome data in `data/reference/{species}/species.fa`

## .fq.gz
## If your fastqfile name is not end with fq.gz,you need to change the I/O of rule "trim_adapter" and "trim_adapter_fastqc" in "rules/qc.smk"
- Place the raw sequencing data in `data/samples/{samples}/samples_1.fq.gz` and samples_2.fq.gz

### 3. Running the Pipeline(The simple Command of Snakemake)
```bash
# Basic run with default settings
snakemake --use-conda

# Run with specific number of cores (replace X with desired number)
snakemake --use-conda --cores X

# Run with dry-run to check workflow
snakemake -np

```

### 4. Common Options
- `--cores X` or `-c X`: Use X CPU cores
- `-np`: Dry run (show what would be done without executing)
- `--resources mem_mb=X`: Limit memory usage to X MB per job

## Configuration Guide

### Configuration Files
1. `configure/config.yaml`: Main configuration file
2. `configure/samples_config.yaml`: Sample configuration file
3. `configure/references_config.yaml`: Reference genome configuration file

### Sample Configuration Requirements
- Species names must be listed in `species_list` and must match the SnpEff database name. If the names do not match, errors will occur during annotation. Alternatively, if you only need the VCF of variants, you can skip annotation and plots by commenting out the lines `include: "rules/annotation.smk"` and `include: "rules/plots.smk"` in the Snakefile using `##`.
- Sample IDs must be correctly configured under their respective species.
- Each sample requires paired-end sequencing file paths.

### Reference Genome Configuration
- Reference genome information is specified in `references_config.yaml`
- Supports multiple species reference genome configurations

## Troubleshooting Guide

### 1. Species/Sample Mismatch Errors
- Verify species names in `species_list` match directory structures.
- Ensure sample IDs are correctly configured in `samples_config.yaml`.

### 2. Missing Input Files
- Check if wildcard patterns match actual file paths.
- Validate regex patterns in `wildcard_constraints`.
- Any errors about the fastq and fasta file:Whether the file type style is ".fq.gz" and ".fa"

### 3. Configuration Issues
- Ensure all required fields exist in configuration files.
- Use online YAML validators to check syntax.

### 4. Chromosome Naming in SnpEff
When using SnpEff for annotation, inconsistencies in chromosome naming between the raw data or reference genome and the SnpEff database can lead to missing feature information in the annotation results. For example, in the case of *Corynebacterium glutamicum ATCC 13032*, the chromosome is labeled as `Chromosome` in the SnpEff database, whereas the raw data may use a format like `BAXXXX`. This mismatch can result in incomplete annotations.

#### Suggested Solutions:
1. **Manual Adjustment**: After running the pipeline, manually edit the chromosome names in the reference genome or annotation results to match the SnpEff database. You can find the expected chromosome names in the file `results/annotation/{species}/snpeff/annotated_{species}.txt` under the section:
   ```
   # Number of chromosomes      : 1
   # Chromosomes                : Format 'chromo_name size codon_table'
   #       'Chromosome' 3282708 Standard
   ```

   After modifying the chromosome information, re-run the snakemake command, and it will automatically start with the annotations.

2. **Custom Database Construction**: Build a custom SnpEff database using your reference genome to ensure compatibility with the chromosome names. This approach is particularly useful for genomes with multiple contigs, where manual adjustments may fail (e.g., *Janibacter melonis*).
3. **Annotation Summary Script**: When encountering annotation display issues, use the integrated `snpeff_plots.py` script to:
   - Generate quick statistical summaries of annotated VCF files
   - Create essential visualizations including:
     - Variant impact distribution
     - Functional region distribution
     - Variant type composition
     - Genome-wide variant density
     - Coverage depth distribution
   - Output location: `results/plots/{species}/snpeff_effects_{species}.png` (visualizations) and `results/plots/{species}/snpeff_report_{species}.txt` (statistical summary)
   - Automatically executes as part of the `snpeff_plots` rule in the pipeline
#### Known Challenges:
For reference genomes composed of multiple contigs, recognition issues may arise even after manual adjustments. In such cases, constructing a custom SnpEff database is recommended, but it requires careful preparation of the reference genome and associated files.

## Output Files
Main output files generated by the pipeline include:
- QC Results: `results/QC/{species}/multiqc/final/`
- Mapping Statistics: `results/Mapping/{species}/{sample}/samtools/{sample}.stats`
- Assembly Evaluation: `results/Assembly/{species}/{sample}/QUAST/.done`
- Visualization Results: `results/{species}/{sample}/depth/{sample}_depth.svg`
- Variant Annotation: `results/plots/{species}/annotation/{species}_amr_variants.tsv`

## Module Description
The pipeline includes the following main modules:
- `rules/qc.smk`: Quality Control
- `rules/mapping.smk`: Read Mapping
- `rules/assembly.smk`: Genome Assembly
- `rules/mk_rm_duplication.smk`: Duplicate Removal
- `rules/SNPIndel.smk`: Variant Detection
- `rules/annotation.smk`: Variant Annotation
- `rules/plots.smk`: Results Simple Analyse

## Project Structure
```
AUTOVAR/
├── configure/                      # Configuration files
│   ├── config.yaml                 # Main configuration
│   ├── samples_config.yaml         # Sample information
│   └── references_config.yaml      # Reference genome paths
│
├── workflow/                       # Pipeline workflow
│   ├── Snakefile                   # Main pipeline file
│   └── rules/                      # Rule modules
│       ├── qc.smk                  # Quality control
│       ├── mapping.smk             # Read mapping
│       ├── assembly.smk            # Genome assembly
│       ├── mk_rm_duplication.smk   # Duplicate removal
│       ├── SNPIndel.smk            # Variant detection
│       ├── annotation.smk          # Variant annotation
│       ├── plots.smk               # Visualization
│       ├── scripts/                # Helper scripts
│       └── envs/                   # Conda environments
│
├── data/                           # Input data directory
│   ├── reference/                  # Reference genomes
│   │   └── {species}/              # Species-specific references
│   │       └── ref.fa              # Reference genome file
│   └── samples/                    # Sample data
│       └── {sample}/               # Sample-specific data
│           ├── r1.fq.gz            # Forward reads
│           └── r2.fq.gz            # Reverse reads
│
├── logs/                           # Log files
├── results/                        # Output results
│   ├── QC/                         # Quality control results
│   ├── Mapping/                    # Mapping results
│   ├── Assembly/                   # Assembly results
│   └── plots/                      # Visualization results
│
├── environment.yaml                # Conda environment specification
├── README.md                       # This documentation
└── LICENSE.md                      # MIT LICENSE
``` 

## Decision to Skip Base Quality Score Recalibration (BQSR)

In this version of the pipeline, Base Quality Score Recalibration (BQSR) is not included because there is no suitable bacterial database available.

### Future Plans
The current `SNPIndel.smk` functionality will be used to build a database for BQSR, controlled by configuration settings. This process will require sequencing data with enough depth to ensure reliable recalibration.

Optimized according to snakemake's best practices.

## References

For detailed references of all tools used in this pipeline, see the [REFERENCES.md](REFERENCES.md) file.

## License
Copyright (c) 2025 FoxyyFace [MIT License](LICENSE.md).
import pandas as pd
from cyvcf2 import VCF
from intervaltree import IntervalTree
import sys

def load_amr_genes(abricate_files):
    """Load AMR gene information and return two structures:
    1. Interval tree for fast position queries
    2. DataFrame of all AMR loci (for reporting undetected cases)
    """
    amr_trees = {}
    all_amr_loci = []
    
    for f in abricate_files:
        try:
            # Skip empty files
            with open(f, 'r') as fh:
                if not any(line.strip() for line in fh if not line.startswith('#')):
                    print(f"Skipping empty file: {f}")
                    continue

            # Read header line
            with open(f, 'r') as fh:
                headers = []
                for line in fh:
                    if line.startswith('#'):
                        headers = line.lstrip('#').strip().split('\t')
                        break
                
                # Validate required columns
                required = ['SEQUENCE', 'START', 'END', 'GENE']
                if not all(col in headers for col in required):
                    print(f"Missing required columns in {f}")
                    continue
                
                # Load data with proper column names
                df = pd.read_csv(f, sep='\t', comment='#', names=headers)
                
                # Standardize column names and coordinates
                df = df.rename(columns={
                    'SEQUENCE': 'chrom',
                    'START': 'start',
                    'END': 'end',
                    'GENE': 'gene'
                })
                df['start'] = df[['start', 'end']].min(axis=1)
                df['end'] = df[['start', 'end']].max(axis=1)
                
                # Populate data structures
                for _, row in df.iterrows():
                    chrom = row['chrom']
                    start = row['start']
                    end = row['end']
                    gene = row['gene']
                    
                    # Add to interval tree
                    if chrom not in amr_trees:
                        amr_trees[chrom] = IntervalTree()
                    amr_trees[chrom].addi(start, end, data=gene)
                    
                    # Record all AMR loci
                    all_amr_loci.append({
                        'Chromosome': chrom,
                        'Start': start,
                        'End': end,
                        'Gene': gene,
                        'Detection': 'Not detected'
                    })
                    
        except Exception as e:
            print(f"Error processing {f}: {str(e)}")
            continue
            
    return amr_trees, pd.DataFrame(all_amr_loci).drop_duplicates()

def main():
    """Main workflow execution"""
    # Input/output paths from Snakemake
    vcf_path = snakemake.input.variants
    abricate_files = snakemake.input.abricate_files
    output_path = snakemake.output[0]
    
    # Load AMR gene information
    amr_trees, amr_loci_df = load_amr_genes(abricate_files)
    
    # Process VCF variants
    results = []
    if amr_trees:
        for variant in VCF(vcf_path):
            chrom = variant.CHROM
            pos = variant.POS
            
            if chrom in amr_trees:
                overlaps = amr_trees[chrom].overlap(pos)
                if overlaps:
                    results.append({
                        'CHROM': chrom,
                        'POS': pos,
                        'REF': variant.REF,
                        'ALT': ','.join(variant.ALT),
                        'FILTER': variant.FILTER or 'PASS',
                        'AMR_GENES': ','.join({iv.data for iv in overlaps})
                    })
    
    # Generate output
    if results:
        # Case 1: Found variants in AMR regions
        pd.DataFrame(results).to_csv(output_path, sep='\t', index=False)
    else:
        # Case 2: No variants found - report AMR loci
        if not amr_loci_df.empty:
            amr_loci_df['Note'] = 'No variants detected in this AMR region'
            amr_loci_df.to_csv(output_path, sep='\t', index=False)
        else:
            # Case 3: No AMR genes detected at all
            with open(output_path, 'w') as f:
                f.write("No AMR genes detected in any samples\n")

if __name__ == "__main__":
    main()
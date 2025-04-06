import pandas as pd
import matplotlib.pyplot as plt
import re

def extract_field(info, field_name):
    # Extract the value of a specific field from the INFO column.
    match = re.search(rf'{field_name}=([^;]+)', info)
    return match.group(1) if match else None

def get_variant_type(ref, alt):
    # Determine the type of variant (Transition, Transversion, or Indel).
    if len(ref) != 1 or len(alt) != 1:
        return 'Indel'
    transitions = [{'A', 'G'}, {'C', 'T'}]
    return 'Transition' if {ref.upper(), alt.upper()} in transitions else 'Transversion'

def plot_all_analysis(df, output_path):
    # Generate combined visualizations for variant analysis.
    plt.figure(figsize=(20, 25))
    
    # Variant Impact Distribution
    plt.subplot(3, 2, 1)
    if 'ANN' in df.columns and not df['ANN'].isnull().all():
        impact_levels = []
        for ann in df['ANN'].dropna():
            for entry in ann.split(','):
                parts = entry.split('|')
                if len(parts) > 2:
                    impact_levels.append(parts[2].strip())
        if impact_levels:
            impact_counts = pd.Series(impact_levels).value_counts().sort_index()
            impact_counts.plot(kind='bar', color='#1f77b4')
            plt.title('Variant Impact Distribution', fontsize=12, pad=20)
            plt.xlabel('Impact Level', fontsize=10)
            plt.ylabel('Count', fontsize=10)
            plt.xticks(rotation=45, ha='right')
        else:
            plt.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
            plt.title('Variant Impact Distribution', fontsize=12, pad=20)

    # Combined Variant Density and Coverage Depth
    plt.subplot(3, 2, 2)
    if 'POS' in df.columns and 'INFO' in df.columns:
        bin_size = 20000  # 0.02M bins
        max_pos = df['POS'].max()
        bins = range(0, max_pos + bin_size, bin_size)
        
        # Calculate variant density
        pos_counts = df.groupby(pd.cut(df['POS'], bins=bins)).size()
        pos_counts.index = [f'{int(b.left)/1e6:.2f}M-{int(b.right)/1e6:.2f}M' for b in pos_counts.index]
        
        # Calculate average coverage depth
        df['DP'] = df['INFO'].apply(lambda x: extract_field(x, 'DP')).astype(float)
        valid_dp = pd.to_numeric(df['DP'], errors='coerce')
        dp_means = valid_dp.groupby(pd.cut(df['POS'], bins=bins)).mean()
        dp_means.index = [f'{int(b.left)/1e6:.2f}M-{int(b.right)/1e6:.2f}M' for b in dp_means.index]

        ax1 = pos_counts.plot(kind='line', linewidth=2, color='#2ca02c', label='Variant Density')
        ax1.set_title('Variant Density and Coverage Depth (Bin Size: 0.02M)', fontsize=12, pad=20)
        ax1.set_xlabel('Genomic Position (M)', fontsize=10)
        ax1.set_ylabel('Variants per Bin', fontsize=10)
        ax1.grid(alpha=0.3)

        ax2 = ax1.twinx()
        dp_means.plot(kind='line', linewidth=2, color='red', label='Average DP', ax=ax2)
        ax2.set_ylabel('Average Depth (DP)', fontsize=10)

        fig = ax1.get_figure()
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='upper left')

    # Functional Regions Distribution
    plt.subplot(3, 2, 3)
    if 'ANN' in df.columns and not df['ANN'].isnull().all():
        feature_types = []
        for ann in df['ANN'].dropna():
            for entry in ann.split(','):
                parts = entry.split('|')
                if len(parts) > 5:
                    feature_types.append(parts[5].strip())
        if feature_types:
            feature_counts = pd.Series(feature_types).value_counts()
            feature_counts.sort_values().plot(kind='barh', color='#ff7f0e')
            plt.title('Functional Regions Distribution', fontsize=12, pad=20)
            plt.xlabel('Count', fontsize=10)
            plt.ylabel('Feature Type', fontsize=10)
        else:
            plt.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
            plt.title('Functional Regions Distribution', fontsize=12, pad=20)

    # Variant Type Composition
    plt.subplot(3, 2, 4)
    if 'REF' in df.columns and 'ALT' in df.columns:
        df['VariantType'] = df.apply(lambda x: get_variant_type(x['REF'], x['ALT']), axis=1)
        type_counts = df['VariantType'].value_counts()
        colors = ['#9467bd', '#8c564b', '#e377c2']
        if not type_counts.empty:
            type_counts.plot(kind='pie', autopct='%1.1f%%', colors=colors,
                            wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
            plt.title('Variant Type Composition', fontsize=12, pad=20)
            plt.ylabel('')
        else:
            plt.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
            plt.title('Variant Type Composition', fontsize=12, pad=20)

    # Coverage Depth Distribution
    plt.subplot(3, 2, 5)
    if 'INFO' in df.columns:
        df['DP'] = df['INFO'].apply(lambda x: extract_field(x, 'DP')).astype(float)
        valid_dp = pd.to_numeric(df['DP'], errors='coerce')
        dp_mean = valid_dp.mean(skipna=True)
        dp_median = valid_dp.median(skipna=True)
        dp_std = valid_dp.std(skipna=True)
        dp_min = valid_dp.min(skipna=True)
        dp_max = valid_dp.max(skipna=True)
        print(f"DP Mean: {dp_mean}, Median: {dp_median}, Std: {dp_std}, Min: {dp_min}, Max: {dp_max}")
        if not valid_dp.dropna().empty:
            valid_dp.dropna().plot(kind='hist', bins=30, color='#17becf', edgecolor='white')
            plt.axvline(dp_mean, color='red', linestyle='dashed', linewidth=2, label=f'Mean: {dp_mean:.1f}')
            plt.axvline(dp_median, color='green', linestyle='dashed', linewidth=2, label=f'Median: {dp_median:.1f}')
            plt.title('Coverage Depth Distribution', fontsize=12, pad=20)
            plt.xlabel('Depth (DP)', fontsize=10)
            plt.ylabel('Frequency', fontsize=10)
            plt.grid(alpha=0.3)
            plt.legend()
        else:
            plt.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
            plt.title('Coverage Depth Distribution', fontsize=12, pad=20)

    plt.tight_layout(pad=4.0)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

def generate_report(df, report_path):
    """Generate a statistical report for the variants."""
    stats = {
        "Total Variants": len(df),
        "Genes Affected": df['ANN'].dropna().str.split('|').str[3].nunique() if 'ANN' in df.columns else 0,
        "Average Depth": f"{df['DP'].mean(skipna=True):.1f}X" if 'DP' in df.columns else "N/A",
        "Median Depth": f"{df['DP'].median(skipna=True):.1f}X" if 'DP' in df.columns else "N/A",
        "Std Deviation Depth": f"{df['DP'].std(skipna=True):.1f}X" if 'DP' in df.columns else "N/A",
        "Min Depth": f"{df['DP'].min(skipna=True):.1f}X" if 'DP' in df.columns else "N/A",
        "Max Depth": f"{df['DP'].max(skipna=True):.1f}X" if 'DP' in df.columns else "N/A",
        "Most Common Effect": df['ANN'].dropna().str.split('|').str[1].value_counts().index[0] 
            if 'ANN' in df.columns and not df['ANN'].dropna().empty else "N/A"
    }

    with open(report_path, 'w') as f:
        f.write("Variant Analysis Report\n")
        f.write("======================\n\n")
        for k, v in stats.items():
            f.write(f"- {k}: {v}\n")

def main(input_vcf, report_path, plot_output_path):
    """Main function to process the VCF file, generate a report, and create visualizations."""
    df = pd.read_csv(input_vcf, sep='\t', comment='#', header=None)
    
    # Extract sample IDs
    sample_ids = df.iloc[0, 9:].tolist()
    
    # Set column names
    fixed_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    all_columns = fixed_columns + sample_ids
    df.columns = all_columns
    
    # Extract ANN field
    def extract_ann(info):
        match = re.search(r'ANN=([^;]+)', info)
        return match.group(1) if match else None
    
    df['ANN'] = df['INFO'].apply(extract_ann)
    
    # Extract DP field
    df['DP'] = df['INFO'].apply(lambda x: extract_field(x, 'DP'))
    
    # Convert DP field to numeric and handle invalid values
    df['DP'] = pd.to_numeric(df['DP'], errors='coerce')
    
    # Debugging: Print DP statistics
    print(f"Number of NaN DP values: {df['DP'].isna().sum()}")
    print(f"Number of non-NaN DP values: {df['DP'].notna().sum()}")
    
    # Debugging: Check ANN field content
    print("Sample ANN entries:")
    print(df['ANN'].head().tolist())
    
    # Generate report
    generate_report(df, report_path)
    
    # Generate visualizations
    plot_all_analysis(df, plot_output_path)

if __name__ == "__main__":
    import sys
    input_vcf = snakemake.input[0]
    report_path = snakemake.output.report
    plot_output_path = snakemake.output.plot
    main(input_vcf, report_path, plot_output_path)




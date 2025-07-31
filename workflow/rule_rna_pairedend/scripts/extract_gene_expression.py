import pandas as pd
import sys
from pathlib import Path

def extract_gene(input_file, output_file, target_gene):
    # 读取计数文件
    df = pd.read_csv(input_file, sep='\t', comment='#')
    
    # 计算总表达量
    total = df.iloc[:, -1].sum()
    
    # 提取目标基因并计算百分比
    target_row = df[df['Geneid'] == target_gene]
    count = target_row.iloc[0, -1] if not target_row.empty else 0
    percentage = (count / total) * 100 if total > 0 else 0.0
    
    # 保存结果
    pd.DataFrame({
        'Gene': [target_gene],
        'Sample': Path(input_file).parent.parent.name,  # 假设路径结构为 .../{sample}/...
        'Percentage': [percentage]
    }).to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    input_file = snakemake.input.counts
    output_file = snakemake.output.expr
    target_gene = snakemake.params.gene_id
    extract_gene(input_file, output_file, target_gene)
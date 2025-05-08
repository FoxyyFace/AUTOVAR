import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
from pathlib import Path
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def calculate_log2fc(merged_df, control_group):
    """计算所有样本相对于对照组的log2(fold change)"""
    # 提取对照组（Human）平均表达量
    control_mean = merged_df[merged_df['source'] == control_group]['Percentage'].mean()
    
    # 计算每个样本的log2FC（对照组自身为0）
    merged_df['log2FC'] = np.log2(merged_df['Percentage'] / control_mean)
    
    # 强制设置对照组log2FC为0（避免浮点误差）
    merged_df.loc[merged_df['source'] == control_group, 'log2FC'] = 0.0
    return merged_df

def calculate_pvalues(merged_df, control_group):
    """计算校正后的p值（假设有生物学重复）"""
    control_data = merged_df[merged_df['source'] == control_group]['Percentage'].values
    results = []
    
    # 按来源分组处理
    for source, group in merged_df[merged_df['source'] != control_group].groupby('source'):
        pvals = []
        for _, row in group.iterrows():
            _, p = ttest_ind(control_data, [row['Percentage']])
            pvals.append(p)
        
        # 校正p值（Benjamini-Hochberg）
        _, adj_pvals, _, _ = multipletests(pvals, method='fdr_bh')
        
        # 合并结果
        for sample, adj_p in zip(group['Sample'], adj_pvals):
            results.append({'Sample': sample, 'adj_pval': adj_p})
    
    # 合并回主表（对照组p值设为1）
    pval_df = pd.DataFrame(results)
    merged_df = pd.merge(merged_df, pval_df, on='Sample', how='left')
    merged_df['adj_pval'] = merged_df['adj_pval'].fillna(1.0)  # 对照组无p值
    return merged_df

def plot_volcano(df, control_group):
    # 定义颜色映射
    color_map = {
        'Human': 'gray',
        'Mouse': 'blue',
        'in_vitro': 'green'
    }
    
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        x='log2FC', 
        y=-np.log10(df['adj_pval']),
        hue='source',
        palette=color_map,
        data=df,
        s=100,
        edgecolor='black',
        linewidth=0.5
    )
    plt.axvline(0, color='black', linestyle='--', linewidth=1)
    plt.title("PA14_46160 Expression: Human vs Other Sources")
    plt.xlabel("log2(Fold Change vs Human)")
    plt.ylabel("-log10(Adjusted P-value)")
    plt.legend(title='Source')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(snakemake.output.plot)
    plt.close()

if __name__ == "__main__":
    # 读取数据
    expr_files = snakemake.input.expr_files
    config = load_config(snakemake.input.config)
    control_group = "Human"
    
    # 合并所有样本数据
    df_list = []
    for file in expr_files:
        data = pd.read_csv(file, sep='\t')
        sample = Path(file).parent.name
        species = Path(file).parent.parent.parent.name
        data['source'] = config['samples'][species][sample]['source']
        df_list.append(data)
    
    merged_df = pd.concat(df_list)
    
    # 计算log2FC和p值
    merged_df = calculate_log2fc(merged_df, control_group)
    merged_df = calculate_pvalues(merged_df, control_group)
    
    # 保存并绘图
    merged_df.to_csv(snakemake.output.table, sep='\t', index=False)
    plot_volcano(merged_df, control_group)
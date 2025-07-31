import os
from collections import defaultdict
import matplotlib.pyplot as plt

def parse_sample_source(filepath, sample_info):
    """从三级配置中解析样本来源"""
    # 从文件路径提取物种和样本ID
    path_parts = filepath.split('/')
    species = path_parts[-3]  # results/kraken2/{species}/{sample}/...
    sample_id = path_parts[-2]
    
    try:
        return sample_info[species][sample_id]['source']
    except KeyError:
        return 'unknown'

def process_reports(input_files, sample_info):
    """处理报告文件，返回增强数据结构"""
    species_data = defaultdict(lambda: {
        'total': 0,
        'samples': defaultdict(int),
        'sources': defaultdict(int)
    })
    
    for file_path in input_files:
        source = parse_sample_source(file_path, sample_info)
        
        with open(file_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6 and parts[3] == 'S':
                    taxon_name = parts[5].strip()
                    count = int(parts[2])
                    
                    # 更新数据
                    species_data[taxon_name]['total'] += count
                    species_data[taxon_name]['samples'][file_path.split('/')[-2]] += count
                    species_data[taxon_name]['sources'][source] += count
    return species_data

def generate_chart(species_names, total_reads, source_data):
    """生成优化版双坐标轴图表"""
    plt.figure(figsize=(18, 12))  # 增加宽度
    
    # ===== 名称处理 =====
    # 方法1：智能换行（每15字符加换行符）
    def wrap_labels(names, max_len=15):
        return ['\n'.join([n[i:i+max_len] for i in range(0, len(n), max_len)]) for n in names]
    
    # 方法2：截断长名称
    processed_names = [n[:20]+'...' if len(n)>20 else n for n in species_names]
    
    # 方法3：添加序号+缩写
    numbered_names = [f"{i+1}. {n[:15]}..." for i, n in enumerate(species_names)]
    
    # ===== 主图表设置 =====
    ax1 = plt.gca()
    bars = ax1.bar(numbered_names,  # 使用带序号的名称
                  [r/1e6 for r in total_reads],  # 转换为百万单位
                  color='#1f77b4',
                  alpha=0.7,
                  width=0.6)  # 调整柱宽
    
    ax1.set_ylabel('Total Reads (Millions)', fontsize=12, color='#1f77b4')
    ax1.tick_params(axis='y', labelcolor='#1f77b4')
    
    # ===== 横坐标优化 =====
    plt.xticks(
        rotation=75,  # 加大旋转角度
        ha='center',  # 水平对齐方式
        fontsize=9,  # 缩小字体
        wrap=True  # 启用自动换行
    )
    
    # ===== 次坐标轴设置 =====
    ax2 = ax1.twinx()
    all_sources = sorted({s for d in source_data for s in d.keys()})
    markers = ['o', 's', 'D', '^', 'v']
    
    # 绘制来源趋势线
    for idx, source in enumerate(all_sources):
        values = [d.get(source, 0)/1e6 for d in source_data]
        ax2.plot(numbered_names, values,
                marker=markers[idx%5],
                linestyle=':',
                linewidth=2,
                markersize=10,
                markeredgecolor='w',  # 标记描边
                label=f'{source}')
    
    ax2.set_ylabel('Source Contribution (M)', fontsize=12, color='#ff7f0e')
    ax2.tick_params(axis='y', labelcolor='#ff7f0e')
    
    # ===== 图表装饰 =====
    plt.title('Species Distribution with Source Contributions\n(Sorted by Total Reads)', 
             fontsize=14, pad=25)
    
    # 添加顶部序号解释
    plt.figtext(0.5, 0.95, 
               "Species Index:\n" + '\n'.join([f"{i+1}. {n}" for i, n in enumerate(species_names)]),
               ha="center",
               fontsize=8,
               bbox=dict(facecolor='#f8f8f8', alpha=0.9))
    
    # 合并图例
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, 
              bbox_to_anchor=(1.15, 1), 
              loc='upper left',
              frameon=False)
    
    # ===== 底部说明 =====
    explanation = (
        "Visual Guide:\n"
        "1. Blue bars: Total sequencing reads (left axis)\n"
        "2. Dotted lines: Source-specific contributions (right axis)\n"
        "3. Full species names see top index\n"
        "4. Mouse-over interactive version available in Jupyter"
    )
    
    plt.figtext(0.5, -0.25, 
               explanation,
               ha="center", 
               fontsize=10,
               bbox=dict(facecolor='#f0f0f0', alpha=0.8))
    
    # 调整布局
    plt.subplots_adjust(
        bottom=0.4,  # 增加底部空间
        top=0.8,     # 增加顶部空间
        left=0.1,    # 增加左侧空间
        right=0.85   # 增加右侧空间
    )
    
    return plt

if __name__ == "__main__":
    # 初始化
    os.makedirs(os.path.dirname(snakemake.output.summary_table), exist_ok=True)
    
    # 处理数据
    data = process_reports(snakemake.input, snakemake.params.sample_info)
    
    # 生成表格
    top_species = sorted(data.items(), key=lambda x: x[1]['total'], reverse=True)[:10]
    with open(snakemake.output.summary_table, 'w') as f:
        f.write("Species\tTotalReads\tTopSamples\tSampleSources\n")
        for species, info in top_species:
            samples = sorted(info['samples'].items(), key=lambda x: x[1], reverse=True)[:3]
            sources = [f"{k}:{v}" for k, v in info['sources'].items()]
            f.write(f"{species}\t{info['total']}\t{'|'.join(f'{s[0]}({s[1]})' for s in samples)}\t{','.join(sources)}\n")
    
    # 生成图表
    species_names = [s[0] for s in top_species]
    total_reads = [s[1]['total'] for s in top_species]
    source_data = [s[1]['sources'] for s in top_species]
    
    plt = generate_chart(species_names, total_reads, source_data)
    plt.savefig(snakemake.output.histogram, dpi=300, bbox_inches='tight')
    plt.close()
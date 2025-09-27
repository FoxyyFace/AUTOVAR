# AUTOVAR - 自动化基因分析流程  
版本: 0.8.0  

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.6.2-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  
[![Python](https://img.shields.io/badge/Python-3.12-3776AB.svg?logo=python&logoColor=white)](https://www.python.org)  

## 概述  
AUTOVAR是一个使用Snakemake搭建的自动化基因分析流程，支持RNA-seq和DNA-seq两种分析模式：  

### RNA-seq模式  
- 质量控制(QC)  
- 读段比对  
- 表达量分析(featureCounts)  

### DNA-seq模式  
- 质量控制(QC)  
- 读段比对  
- 基因组组装  
- 变异检测(SNP&Indel)  
- 变异注释(snpEff)/抗生素抗性检测(abricate)  
- 结果可视化  

## 安装指南  

### 1. 环境设置  
```bash  
# 克隆仓库  
git clone https://github.com/yourusername/AUTOVAR.git  
cd AUTOVAR  

# 创建conda环境  
mamba env create --name AUTOVAR --file environment.yaml  

# 激活环境  
conda activate AUTOVAR  

# 安装SLURM插件（集群环境下）  
pip install \  
  snakemake-executor-plugin-slurm \  
  snakemake-executor-plugin-slurm-jobstep==0.2.1  
```  

### 2. 配置文件准备  
在运行前需要配置以下文件：  

**configure/config.yaml** - 主配置文件  
```yaml  
global:  
  mode: "RNA"  # 或 "DNA"  
  # 其他全局参数...  
```  

**configure/samples_config.yaml** - 样本配置文件  
```yaml  
species_list:  
  - species1  
  - species2  

samples:  
  species1:  
    sample1:  
      fq_paths:  
        - "data/samples/species1/sample1/sample1_1.fq.gz"  
        - "data/samples/species1/sample1/sample1_2.fq.gz"  
    sample2:  
      fq_paths:  
        - "data/samples/species1/sample2/sample2_1.fq.gz"  
        - "data/samples/species1/sample2/sample2_2.fq.gz"  
  species2:  
    sample3:  
      fq_paths:  
        - "data/samples/species2/sample3/sample3_1.fq.gz"  
        - "data/samples/species2/sample3/sample3_2.fq.gz"  
```  

**configure/references_config.yaml** - 参考基因组配置  
```yaml  
reference:  
  species1:  
    fasta: "data/reference/species1/ref.fa"  
    gff: "data/reference/species1/ref.gff"  
  species2:  
    fasta: "data/reference/species2/ref.fa"  
    gff: "data/reference/species2/ref.gff"  
```  

### 3. 辅助脚本使用  
流程包含两个辅助脚本，用于环境设置和样本配置：  

**自动设置Conda环境**  
```bash  
# 安装所有环境  
python auto_set_conda_env.py -y  

# 选择性安装环境 (推荐)
python auto_set_conda_env.py  
```  
此脚本会扫描workflow目录下的环境配置文件，并允许用户选择安装所有或部分环境。  

DNA-seq模式

​无需提前配置环境​
Snakemake会自动在运行时下载所需环境
环境会保存在运行目录的.snakemake子目录中

**自动更新样本配置**  
```bash  
# 自动扫描并更新样本配置  
python auto_update_samples_config.py  
```  
此脚本会自动扫描data/samples目录下的样本文件，并更新samples_config.yaml文件。  

### 4. 数据目录结构  
```  
AUTOVAR/  
├── data/  
│   ├── reference/          # 参考基因组  
│   │   ├── species1/  
│   │   │   ├── ref.fa      # 参考基因组FASTA文件  
│   │   │   └── ref.gff     # 参考基因组GFF注释文件  
│   │   └── species2/  
│   │       ├── ref.fa  
│   │       └── ref.gff  
│   └── samples/            # 样本数据  
│       ├── species1/       # 物种1  
│       │   ├── sample1/    # 样本1  
│       │   │   ├── sample1_1.fq.gz  
│       │   │   └── sample1_2.fq.gz  
│       │   └── sample2/    # 样本2  
│       │       ├── sample2_1.fq.gz  
│       │       └── sample2_2.fq.gz  
│       └── species2/       # 物种2  
│           ├── sample3/    # 样本3  
│           │   ├── sample3_1.fq.gz  
│           │   └── sample3_2.fq.gz  
│           └── sample4/    # 样本4  
│               ├── sample4_1.fq.gz  
│               └── sample4_2.fq.gz  
```  

**重要文件命名规则**  
- 所有文件和目录名称只能使用**英文字母**、**数字**和**下划线(_)**  
- 禁止使用中文、空格、特殊字符或其他非ASCII字符
- 注意目录名称sampleX需要和文件名sampleX_1.fq.gz中的部分一致
- 示例：  
  - 有效：`sample_1.fq.gz`, `ref_genome.fa`, `annotation.gff`  
  - 无效：`样本1.fq.gz`, `ref genome.fa`, `参考基因组.fa`  

## 使用指南  

### 1. RNA-seq模式运行  
```bash  
# 在config.yaml中设置mode: "RNA"  
snakemake --use-conda --cores 16  
```  

### 2. DNA-seq模式运行  
```bash  
# 在config.yaml中设置mode: "DNA"  
snakemake --use-conda --cores 16  
```  

### 3. 集群运行(SLURM)  
```bash  
# 修改profile_slurm/config.yaml中的参数  
vim profile_slurm/config.yaml  

# 启动后台任务 (家目录直接运行,无需srun)
nohup snakemake --profile ./profile_slurm --use-conda --rerun-incomplete > snakemake.log 2>&1 &  

# 监控进度  
tail -f snakemake.log  

# 查看任务提交情况
squeue -u $USER  
```  

### 4. 常用命令选项  
- `--cores N`: 指定使用的CPU核心数  
- `-n`或`--dry-run`: 干运行（不实际执行）  
- `--rerun-incomplete`: 重新运行未完成的任务  
- `--unlock`: 解锁工作目录  
- `--cleanup-metadata`: 清理元数据  

## 输出结果  

### RNA-seq模式输出  
```  
results/  
├── QC/  
├── Mapping/  
└── Analyse/  
```  

### DNA-seq模式输出  
```  
results/  
├── QC/  
├── Mapping/  
├── Assembly/  
├── plots/  
└── annotation/  
```  

## 常见问题  

### 1. 物种/样本配置错误  
确保：  
- `species_list`中的物种名称与目录结构匹配  
- 参考基因组路径在`references_config.yaml`中正确设置  

### 2. 文件路径问题  
- FASTQ文件路径：`data/samples/{species}/{sample}/{sample}_1.fq.gz`  
- 参考基因组FASTA路径：`data/reference/{species}/ref.fa`  
- 参考基因组GFF路径：`data/reference/{species}/ref.gff`  

### 3. 文件名格式问题  
- **所有文件和目录名称必须使用英文字母、数字和下划线**  
- 避免使用中文、空格或其他特殊字符  
- 如果遇到错误，请检查文件名是否符合要求  

### 4. GFF文件缺失问题    
- 确保每个物种参考基因组目录下都有对应的GFF文件  
- GFF文件命名必须与配置文件中指定的名称一致  

### 5. SnpEff注释问题  
当出现染色体命名不一致时：  
1. 检查`results/annotation/{species}/snpeff/annotated_{species}.txt`中的染色体名称  
2. 手动调整参考基因组中的染色体名称  
3. 重新运行流程  

## 许可证  
MIT License - 详见[LICENSE](LICENSE)文件

# RNA-seq Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/Snakemake-7.0+-green)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](https://www.python.org)
[![R](https://img.shields.io/badge/R-4.3-blue)](https://www.r-project.org)

基于 **Snakemake** 构建的自动化 RNA-seq 分析流程，支持从原始 FASTQ 数据到差异表达基因的端到端分析。流程采用模块化设计，具有良好的可扩展性和可重复性。

---

## 流程图

```
FASTQ  -->  FastQC  -->  STAR  -->  featureCounts  -->  DESeq2  -->  差异基因
                |                                                |
                +-->  MultiQC (汇总报告)                          +--> 火山图/热图
```

## 环境要求

| 软件/工具 | 版本 | 用途 |
|-----------|------|------|
| Python | ≥ 3.9 | 流程框架 |
| Snakemake | ≥ 7.0 | 工作流管理 |
| FastQC | ≥ 0.12 | 测序数据质控 |
| STAR | ≥ 2.7 | 序列比对 |
| featureCounts (subread) | ≥ 2.0 | 基因定量 |
| MultiQC | ≥ 1.14 | 质控汇总 |
| R + DESeq2 | ≥ 4.3 | 差异表达分析 |
| samtools | ≥ 1.18 | BAM 处理 |

## 快速开始

### 1. 安装环境

```bash
# 使用 Conda 一键创建环境
conda env create -f envs/environment.yaml
conda activate rna_seq_pipeline
```

### 2. 配置

编辑 `config/config.yaml`，设置参考基因组路径和分析参数。

### 3. 准备数据

将 FASTQ 文件放入 `data/` 目录：

```
data/
├── control_1_R1.fastq.gz
├── control_2_R1.fastq.gz
├── treat_1_R1.fastq.gz
└── treat_2_R1.fastq.gz
```

### 4. 运行

```bash
# 预览流程
snakemake -n

# 运行（8 核）
snakemake -c 8
```

### 5. 查看结果

- `results/multiqc_report.html` — 质控汇总报告
- `results/differential/deseq2_results.csv` — 差异基因列表
- `results/differential/volcano_plot.png` — 火山图
- `results/differential/heatmap.png` — 表达热图

## 目录结构

```
├── Snakefile              # 主流程
├── config/
│   └── config.yaml        # 配置文件
├── envs/
│   └── environment.yaml   # Conda 环境
├── scripts/
│   ├── deseq2_analysis.R  # 差异分析 R 脚本
│   └── merge_counts.py    # 合并 counts
├── docs/
│   ├── usage.md           # 详细使用指南
│   └── developer.md       # 开发者文档
├── data/                  # FASTQ 数据
├── results/               # 分析结果
└── logs/                  # 运行日志
```

## 许可证

本项目采用 MIT 许可证 — 详见 [LICENSE](LICENSE)。

## 作者

- **Hu Yangcheng** — [HuYangcheng-hub](https://github.com/HuYangcheng-hub)

# RNA-seq 转录组分析流程

基于 Snakemake 构建的端到端 bulk RNA-seq 分析流程，覆盖从 FASTQ 质控到 DESeq2 差异表达分析的全流程。

## 流程步骤

```
FASTQ → 输入校验 → FastQC → STAR 比对 → featureCounts 定量 → DESeq2 差异分析
                  → MultiQC 汇总        → 表达矩阵合并        → 火山图/热图/PCA
```

## 快速开始

```bash
# 1. 安装环境
conda env create -f envs/environment.yaml
conda activate rna_seq_pipeline

# 2. 配置样本
# 编辑 config/samples.tsv，添加你的样本信息

# 3. 配置参考基因组路径
# 编辑 config/config.yaml，设置参考基因组和 GTF 路径

# 4. 预览流程
snakemake --use-conda -n

# 5. 运行
snakemake --use-conda -c 8
```

## 输出文件

- `results/multiqc_report.html` — 质控汇总报告
- `results/counts/all_samples_counts.txt` — 基因表达矩阵
- `results/differential/deseq2_results.csv` — 差异表达结果
- `results/differential/volcano_plot.png` — 火山图
- `results/differential/heatmap.png` — 表达热图
- `results/differential/pca_plot.png` — PCA 图

## 技术栈

Snakemake · FastQC · STAR · featureCounts · DESeq2 · R · Python

## 许可证

MIT © Hu Yangcheng

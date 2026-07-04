# RNA-seq 转录组分析流程

[![Test RNA-seq Pipeline](https://github.com/HuYangcheng-hub/rna-seq-pipeline/actions/workflows/test.yml/badge.svg)](https://github.com/HuYangcheng-hub/rna-seq-pipeline/actions/workflows/test.yml)

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

# 2. 下载测试数据（chr22 子集，约 200 MB）
bash scripts/download_test_data.sh

# 3. 准备参考基因组（仅需执行一次）
#    参见 reference/README.md

# 4. 配置样本
#    编辑 config/samples.tsv，添加你的样本信息

# 5. 预览流程
snakemake --use-conda -n

# 6. 运行（8 个核心）
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

Snakemake · FastQC · STAR · featureCounts · DESeq2 · R · Python · Docker · Conda

## 容器运行

```bash
# Docker
docker build -t rna_seq_pipeline .
docker run -it --rm -v $(pwd)/data:/work/data -v $(pwd)/reference:/work/reference -v $(pwd)/results:/work/results rna_seq_pipeline

# Singularity/Apptainer
snakemake --use-singularity -c 8
```

## 许可证

MIT © Hu Yangcheng

# RNA-seq Analysis Pipeline

基于 Snakemake 构建的 RNA-seq 自动分析流程。

## 流程步骤

1. **FastQC** — 原始数据质控
2. **STAR** — 比对到参考基因组
3. **samtools** — BAM 索引
4. **featureCounts** — 基因定量
5. **MultiQC** — 汇总报告

## 使用

```bash
# 预览流程
snakemake -n

# 运行（4 核）
snakemake -c4

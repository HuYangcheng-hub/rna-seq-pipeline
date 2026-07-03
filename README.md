# RNA-seq 转录组分析流程

这是一个基于 Snakemake 的 bulk RNA-seq 分析流程，覆盖从测序数据上游处理到差异表达分析的主要步骤。流程支持双端或单端 FASTQ 输入，包含样本与配置校验、FastQC 质控、可选 fastp 过滤、STAR 比对、featureCounts 定量、表达矩阵合并以及 DESeq2 差异分析。

## 流程概览

```text
samples.tsv
    |
输入与配置校验
    |
FastQC --------\
                > MultiQC
fastp（可选）
    |
STAR -> samtools index -> featureCounts -> 合并 counts -> DESeq2
```

## 主要结果

- `results/multiqc_report.html`：质控汇总报告
- `results/star/{sample}_Aligned.sortedByCoord.out.bam`：STAR 比对后的排序 BAM
- `results/counts/{sample}_counts.txt`：单样本 featureCounts 结果
- `results/counts/all_samples_counts.txt`：所有样本的基因表达矩阵
- `results/differential/deseq2_results.csv`：DESeq2 全量差异分析结果
- `results/differential/deseq2_sig_genes.csv`：显著差异基因列表
- `results/differential/normalized_counts.csv`：归一化表达矩阵
- `results/differential/volcano_plot.png`：火山图
- `results/differential/heatmap.png`：表达热图
- `results/differential/pca_plot.png`：PCA 图

## 快速开始

```bash
conda env create -f envs/environment.yaml
conda activate rna_seq_pipeline
snakemake --use-conda -n
snakemake --use-conda -c 8
```

运行真实数据前，请先修改：

- `config/config.yaml`：参考基因组、分析参数、是否双端、是否启用 fastp 等配置
- `config/samples.tsv`：样本名、分组信息和 FASTQ 路径

Snakemake 规则默认使用 `envs/` 下拆分好的小环境；`envs/environment.yaml` 提供给手动调试或一次性安装完整环境使用。

## 样本表格式

`config/samples.tsv` 是制表符分隔文件：

```text
sample	condition	fastq_1	fastq_2
control_1	control	data/control_1_R1.fastq.gz	data/control_1_R2.fastq.gz
control_2	control	data/control_2_R1.fastq.gz	data/control_2_R2.fastq.gz
treat_1	treat	data/treat_1_R1.fastq.gz	data/treat_1_R2.fastq.gz
treat_2	treat	data/treat_2_R1.fastq.gz	data/treat_2_R2.fastq.gz
```

如果是单端数据，请在 `config/config.yaml` 中设置：

```yaml
paired_end: false
```

同时将 `fastq_2` 留空。

## 参考基因组配置

请在 `config/config.yaml` 中配置 STAR index、GTF 和 genome FASTA：

```yaml
reference:
  star_index: "reference/GRCh38/star_index"
  gtf: "reference/GRCh38/gencode.v44.annotation.gtf"
  genome_fasta: "reference/GRCh38/GRCh38.p14.genome.fa"
```

大型 FASTQ、BAM、结果文件和参考基因组文件默认不会提交到 git。

## 注意事项

- 差异分析至少需要两个分组，建议每组至少 2 个生物学重复，最好 3 个以上。
- 如果希望重复数不足时直接终止流程，可在 `config/config.yaml` 中设置 `strict_replicates: true`。
- 如果 STAR 输出读入数为 0，请优先检查 FASTQ 路径、压缩格式、文件内容和 STAR index。
- 如果 featureCounts 分配 reads 为 0，请检查 BAM 与 GTF 的染色体命名是否一致，以及链特异性参数是否正确。
- 详细使用说明见 `docs/usage.md`，开发维护说明见 `docs/developer.md`。

## 许可证

本项目使用 MIT 许可证，详见 `LICENSE`。

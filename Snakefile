
"""
======================================================================
RNA-seq 分析流程 (RNA-seq Analysis Pipeline)
======================================================================

本流程实现从 FASTQ 原始数据到差异表达基因 (DEG) 的全自动分析。

流程步骤:
  1. FastQC       — 原始测序数据质量评估
  2. MultiQC      — 多样本质控报告汇总
  3. STAR         — 比对到参考基因组
  4. samtools     — BAM 文件排序与索引
  5. featureCounts — 基因水平定量
  6. merges_counts — 合并表达矩阵
  7. DESeq2       — 差异表达分析 & 可视化

使用方法:
  snakemake -c N         # 使用 N 个核心运行
  snakemake -n           # 预览流程（不实际运行）
  snakemake --dag | dot -Tpng > dag.png  # 生成流程图

作者: Hu Yangcheng
许可证: MIT
"""

configfile: "config/config.yaml"

# =======================
# 全局参数
# =======================
SAMPLES = config["samples"]
PAIRED_END = config.get("paired_end", True)
THREADS = config.get("resources", {})

# 辅助函数：获取 R1/R2 文件路径
def get_fastqs(sample, wildcards):
    r1 = f"data/{sample}_R1.fastq.gz"
    return r1

# =======================
# 最终目标
# =======================
rule all:
    """最终产出：质控报告 + 差异表达结果"""
    input:
        "results/multiqc_report.html",
        "results/differential/deseq2_results.csv"

# =======================
# Step 1: 原始数据质控
# =======================
rule fastqc:
    """FastQC: 对原始 FASTQ 进行质量评估
       输出 HTML 可视化报告 + ZIP 压缩包"""
    input:
        get_fastqs
    output:
        html = "results/fastqc/{sample}_fastqc.html",
        zip = "results/fastqc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    threads: THREADS.get("fastqc_threads", 4)
    shell:
        "fastqc -o results/fastqc --threads {threads} {input} 2> {log}"

# =======================
# Step 2: 质控汇总
# =======================
rule multiqc:
    """MultiQC: 汇总所有 FastQC 报告为一个交互式 HTML"""
    input:
        expand("results/fastqc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        "results/multiqc_report.html"
    log:
        "logs/multiqc.log"
    shell:
        "multiqc results/fastqc -o results/ 2> {log}"

# =======================
# Step 3: STAR 比对
# =======================
rule star:
    """STAR: 将 reads 比对到参考基因组
       输出按坐标排序的 BAM 文件"""
    input:
        r1 = "data/{sample}_R1.fastq.gz"
    output:
        bam = "results/star/{sample}_Aligned.sortedByCoord.out.bam",
        log_final = "results/star/{sample}_Log.final.out"
    params:
        index = config["star_index"],
        prefix = "results/star/{sample}_"
    log:
        "logs/star/{sample}.log"
    threads: THREADS.get("star_threads", 8)
    run:
        cmd = (
            f"STAR --genomeDir {params.index} "
            f"--readFilesIn {input.r1} "
        )
        if PAIRED_END:
            cmd += f"data/{{wildcards.sample}}_R2.fastq.gz "
        cmd += (
            f"--readFilesCommand zcat "
            f"--outFileNamePrefix {params.prefix} "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--runThreadN {threads} 2> {log}"
        )
        shell(cmd)

# =======================
# Step 4: BAM 索引
# =======================
rule samtools_index:
    """samtools: 为 BAM 文件建立索引 (.bai)"""
    input:
        "results/star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}_Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/samtools/{sample}.log"
    shell:
        "samtools index {input} 2> {log}"

# =======================
# Step 5: 基因定量
# =======================
rule feature_counts:
    """featureCounts: 在基因水平进行 reads 计数"""
    input:
        bam = "results/star/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "results/star/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        "results/counts/{sample}_counts.txt"
    params:
        gtf = config["gtf"],
        strand = config.get("feature_counts_params", {}).get("strand_specific", 0),
        min_qual = config.get("feature_counts_params", {}).get("min_mapping_quality", 10)
    log:
        "logs/featureCounts/{sample}.log"
    threads: THREADS.get("featurecounts_threads", 4)
    shell:
        "featureCounts -T {threads} -t exon -g gene_id "
        "-s {params.strand} -Q {params.min_qual} "
        "-a {params.gtf} -o {output} {input.bam} 2> {log}"

# =======================
# Step 6: 合并表达矩阵
# =======================
rule merge_counts:
    """合并所有样本的 count 数据为一个表达矩阵"""
    input:
        expand("results/counts/{sample}_counts.txt", sample=SAMPLES)
    output:
        "results/counts/all_samples_counts.txt"
    log:
        "logs/merge_counts.log"
    run:
        import pandas as pd
        import glob

        count_files = sorted(glob.glob("results/counts/*_counts.txt"))
        counts = []
        for f in count_files:
            # 跳过合并后的文件自身
            if "all_samples" in f:
                continue
            # 读取 featureCounts 输出（跳过前 2 行注释）
            df = pd.read_csv(f, sep="\t", comment="#", low_memory=False)
            sample_name = f.split("/")[-1].replace("_counts.txt", "")

            # 提取基因名和 counts
            if "Geneid" in df.columns and len(df.columns) >= 7:
                count_col = df.columns[6]
                counts.append(df[["Geneid", count_col]].rename(
                    columns={"Geneid": "gene_id", count_col: sample_name}
                ))

        if counts:
            merged = counts[0]
            for df in counts[1:]:
                merged = merged.merge(df, on="gene_id", how="outer")
            merged.fillna(0, inplace=True)
            merged.to_csv("results/counts/all_samples_counts.txt", sep="\\t", index=False)

# =======================
# Step 7: 差异表达分析
# =======================
rule deseq2_analysis:
    """DESeq2 R 脚本: 差异表达分析 + 火山图 + 热图"""
    input:
        "results/counts/all_samples_counts.txt"
    output:
        csv = "results/differential/deseq2_results.csv",
        sig = "results/differential/deseq2_sig_genes.csv",
        volcano = "results/differential/volcano_plot.png",
        heatmap = "results/differential/heatmap.png"
    log:
        "logs/deseq2/deseq2.log"
    script:
        "scripts/deseq2_analysis.R"

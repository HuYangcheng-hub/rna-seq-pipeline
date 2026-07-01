"""
RNA-seq 分析流程
=================
输入: FASTQ 文件 → 质控 → 比对 → 定量 → 质控汇总
"""

configfile: "config/config.yaml"

# 获取样本名列表
SAMPLES = config["samples"]

rule all:
    """最终输出：MultiQC 汇总报告"""
    input:
        expand("results/{sample}/{sample}_sorted.bam", sample=SAMPLES),
        "results/multiqc_report.html"

# ---------- step 1: 质控 ----------
rule fastqc:
    """对原始 FASTQ 进行质量检查"""
    input:
        "data/{sample}_R1.fastq.gz"
    output:
        html = "results/fastqc/{sample}_fastqc.html",
        zip = "results/fastqc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    threads: 4
    shell:
        "fastqc -o results/fastqc --threads {threads} {input} 2> {log}"

# ---------- step 2: 比对 ----------
rule star_alignment:
    """STAR 比对：将 reads 比对到参考基因组"""
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = lambda wildcards: f"data/{wildcards.sample}_R2.fastq.gz" \
            if config.get("paired_end", True) else []
    output:
        bam = "results/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        index = config["star_index"],
        prefix = "results/{sample}/{sample}_"
    log:
        "logs/star/{sample}.log"
    threads: 8
    shell:
        "STAR --genomeDir {params.index} "
        "--readFilesIn {input.r1} {input.r2} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM SortedByCoordinate "
        "--runThreadN {threads} "
        "2> {log}"

# ---------- step 3: 索引 BAM ----------
rule index_bam:
    """为 BAM 文件建立索引"""
    input:
        "results/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "results/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/samtools/{sample}.log"
    shell:
        "samtools index {input} 2> {log}"

# ---------- step 4: 基因定量 ----------
rule feature_counts:
    """featureCounts 进行基因水平定量"""
    input:
        bam = "results/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "results/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        "results/{sample}/{sample}_counts.txt"
    params:
        gtf = config["gtf"]
    log:
        "logs/featureCounts/{sample}.log"
    threads: 4
    shell:
        "featureCounts -T {threads} -t exon -g gene_id "
        "-a {params.gtf} -o {output} {input.bam} "
        "2> {log}"

# ---------- step 5: MultiQC 汇总 ----------
rule multiqc:
    """汇总所有质控报告"""
    input:
        expand("results/fastqc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        "results/multiqc_report.html"
    log:
        "logs/multiqc.log"
    shell:
        "multiqc results/fastqc -o results/ 2> {log}"

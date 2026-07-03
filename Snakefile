"""
RNA-seq upstream and differential expression workflow.

The workflow covers:
  1. input/config validation
  2. FastQC and MultiQC
  3. optional fastp trimming
  4. STAR alignment
  5. BAM indexing
  6. featureCounts quantification
  7. count matrix merge
  8. DESeq2 differential analysis
"""

import csv
from pathlib import Path

configfile: "config/config.yaml"


SAMPLE_SHEET = config.get("sample_sheet", "config/samples.tsv")
PAIRED_END = bool(config.get("paired_end", True))
TRIM_ENABLED = bool(config.get("trimming", {}).get("enabled", False))
READS = ["R1", "R2"] if PAIRED_END else ["R1"]


def read_samples(sample_sheet):
    samples = {}
    with open(sample_sheet, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "condition", "fastq_1"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{sample_sheet} is missing required columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            sample = row["sample"].strip()
            if not sample:
                continue
            samples[sample] = {
                "condition": row["condition"].strip(),
                "R1": row["fastq_1"].strip(),
                "R2": row.get("fastq_2", "").strip(),
            }
    if not samples:
        raise ValueError(f"No samples found in {sample_sheet}")
    return samples


SAMPLE_INFO = read_samples(SAMPLE_SHEET)
SAMPLES = list(SAMPLE_INFO)
CONDITIONS = sorted({SAMPLE_INFO[s]["condition"] for s in SAMPLES})


def raw_fastq(wildcards):
    return SAMPLE_INFO[wildcards.sample][wildcards.read]


def alignment_reads(wildcards):
    if TRIM_ENABLED:
        if PAIRED_END:
            return [
                f"results/trimmed/{wildcards.sample}_R1.fastq.gz",
                f"results/trimmed/{wildcards.sample}_R2.fastq.gz",
            ]
        return [f"results/trimmed/{wildcards.sample}_R1.fastq.gz"]
    reads = [SAMPLE_INFO[wildcards.sample]["R1"]]
    if PAIRED_END:
        reads.append(SAMPLE_INFO[wildcards.sample]["R2"])
    return reads


def star_read_command(wildcards, input):
    reads = list(input.reads)
    return "--readFilesCommand zcat" if any(str(read).endswith(".gz") for read in reads) else ""


def star_extra_params():
    params = config.get("star_params", {})
    return " ".join(f"--{key} {value}" for key, value in params.items())


def featurecounts_paired_params():
    return "-p --countReadPairs" if PAIRED_END else ""


FASTQC_OUTPUTS = expand(
    "results/fastqc/{sample}_{read}_fastqc.zip",
    sample=SAMPLES,
    read=READS,
)

TRIM_OUTPUTS = (
    expand("results/trimmed/{sample}_{read}.fastq.gz", sample=SAMPLES, read=READS)
    if TRIM_ENABLED
    else []
)


rule all:
    input:
        "results/validation/samples.validated.tsv",
        "results/multiqc_report.html",
        expand("results/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("results/star/{sample}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand("results/counts/{sample}_counts.txt", sample=SAMPLES),
        "results/counts/all_samples_counts.txt",
        "results/differential/deseq2_results.csv",
        "results/differential/deseq2_sig_genes.csv",
        "results/differential/normalized_counts.csv",
        "results/differential/volcano_plot.png",
        "results/differential/heatmap.png",
        "results/differential/pca_plot.png",


rule validate_inputs:
    input:
        sample_sheet=SAMPLE_SHEET,
        fastqs=[SAMPLE_INFO[s][read] for s in SAMPLES for read in READS],
    output:
        "results/validation/samples.validated.tsv",
    log:
        "logs/validation/validate_inputs.log",
    conda:
        "envs/validate.yaml",
    shell:
        "python scripts/validate_inputs.py "
        "--config config/config.yaml "
        "--sample-sheet {input.sample_sheet} "
        "--output {output} > {log} 2>&1"


rule fastqc:
    input:
        fastq=raw_fastq,
        validated="results/validation/samples.validated.tsv",
    output:
        html="results/fastqc/{sample}_{read}_fastqc.html",
        zip="results/fastqc/{sample}_{read}_fastqc.zip",
    log:
        "logs/fastqc/{sample}_{read}.log",
    threads: config.get("resources", {}).get("fastqc_threads", 2)
    conda:
        "envs/qc.yaml",
    shell:
        "fastqc -o results/fastqc --threads {threads} {input.fastq} > {log} 2>&1"


rule fastp_trim:
    input:
        r1=lambda wc: SAMPLE_INFO[wc.sample]["R1"],
        r2=lambda wc: SAMPLE_INFO[wc.sample]["R2"] if PAIRED_END else "",
        validated="results/validation/samples.validated.tsv",
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz" if PAIRED_END else [],
        html="results/trimmed/{sample}.fastp.html",
        json="results/trimmed/{sample}.fastp.json",
    log:
        "logs/fastp/{sample}.log",
    threads: config.get("resources", {}).get("trim_threads", 4)
    conda:
        "envs/qc.yaml",
    script:
        "scripts/run_fastp.py"


rule multiqc:
    input:
        fastqc=FASTQC_OUTPUTS,
        trimmed=TRIM_OUTPUTS,
    output:
        "results/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "envs/qc.yaml",
    shell:
        "multiqc results -o results --filename multiqc_report.html > {log} 2>&1"


rule star_align:
    input:
        reads=alignment_reads,
        validated="results/validation/samples.validated.tsv",
    output:
        bam="results/star/{sample}_Aligned.sortedByCoord.out.bam",
        log_final="results/star/{sample}_Log.final.out",
    params:
        index=lambda wc: config["reference"]["star_index"],
        prefix=lambda wc, output: str(output.bam).replace("Aligned.sortedByCoord.out.bam", ""),
        read_command=star_read_command,
        extra=star_extra_params(),
    log:
        "logs/star/{sample}.log",
    threads: config.get("resources", {}).get("star_threads", 8)
    conda:
        "envs/align.yaml",
    shell:
        "STAR --genomeDir {params.index} "
        "--readFilesIn {input.reads} "
        "{params.read_command} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM SortedByCoordinate "
        "--runThreadN {threads} "
        "{params.extra} > {log} 2>&1"


rule samtools_index:
    input:
        "results/star/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools/{sample}.log",
    conda:
        "envs/align.yaml",
    shell:
        "samtools index {input} > {log} 2>&1"


rule feature_counts:
    input:
        bam="results/star/{sample}_Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
        validated="results/validation/samples.validated.tsv",
    output:
        "results/counts/{sample}_counts.txt",
    params:
        gtf=lambda wc: config["reference"]["gtf"],
        strand=lambda wc: config.get("feature_counts_params", {}).get("strand_specific", 0),
        min_qual=lambda wc: config.get("feature_counts_params", {}).get("min_mapping_quality", 10),
        paired=featurecounts_paired_params(),
        extra=lambda wc: config.get("feature_counts_params", {}).get("extra", ""),
    log:
        "logs/featureCounts/{sample}.log",
    threads: config.get("resources", {}).get("featurecounts_threads", 4)
    conda:
        "envs/counts.yaml",
    shell:
        "featureCounts -T {threads} -t exon -g gene_id "
        "-s {params.strand} -Q {params.min_qual} "
        "{params.paired} {params.extra} "
        "-a {params.gtf} -o {output} {input.bam} > {log} 2>&1"


rule merge_counts:
    input:
        counts=expand("results/counts/{sample}_counts.txt", sample=SAMPLES),
        sample_sheet=SAMPLE_SHEET,
    output:
        matrix="results/counts/all_samples_counts.txt",
    log:
        "logs/merge_counts.log",
    conda:
        "envs/counts.yaml",
    script:
        "scripts/merge_counts.py"


rule deseq2_analysis:
    input:
        counts="results/counts/all_samples_counts.txt",
        sample_sheet=SAMPLE_SHEET,
    output:
        csv="results/differential/deseq2_results.csv",
        sig="results/differential/deseq2_sig_genes.csv",
        normalized="results/differential/normalized_counts.csv",
        volcano="results/differential/volcano_plot.png",
        heatmap="results/differential/heatmap.png",
        pca="results/differential/pca_plot.png",
    log:
        "logs/deseq2/deseq2.log",
    conda:
        "envs/deseq2.yaml",
    script:
        "scripts/deseq2_analysis.R"

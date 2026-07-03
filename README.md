# RNA-seq Analysis Pipeline

Snakemake workflow for bulk RNA-seq upstream analysis through differential
expression. It supports paired-end or single-end FASTQ input, quality control,
optional trimming, STAR alignment, featureCounts quantification, count matrix
generation, and DESeq2 analysis.

## Workflow

```text
samples.tsv
    |
validate inputs
    |
FastQC --------\
                > MultiQC
fastp (optional)
    |
STAR -> samtools index -> featureCounts -> merge counts -> DESeq2
```

## Main Outputs

- `results/multiqc_report.html`
- `results/star/{sample}_Aligned.sortedByCoord.out.bam`
- `results/counts/{sample}_counts.txt`
- `results/counts/all_samples_counts.txt`
- `results/differential/deseq2_results.csv`
- `results/differential/deseq2_sig_genes.csv`
- `results/differential/normalized_counts.csv`
- `results/differential/volcano_plot.png`
- `results/differential/heatmap.png`
- `results/differential/pca_plot.png`

## Quick Start

```bash
conda env create -f envs/environment.yaml
conda activate rna_seq_pipeline
snakemake --use-conda -n
snakemake --use-conda -c 8
```

Edit `config/config.yaml` and `config/samples.tsv` before running real data.
Snakemake rules use smaller per-step environments under `envs/`; the combined
`envs/environment.yaml` is provided for interactive use.

## Sample Sheet

`config/samples.tsv` is tab-delimited:

```text
sample	condition	fastq_1	fastq_2
control_1	control	data/control_1_R1.fastq.gz	data/control_1_R2.fastq.gz
control_2	control	data/control_2_R1.fastq.gz	data/control_2_R2.fastq.gz
treat_1	treat	data/treat_1_R1.fastq.gz	data/treat_1_R2.fastq.gz
treat_2	treat	data/treat_2_R1.fastq.gz	data/treat_2_R2.fastq.gz
```

For single-end data, set `paired_end: false` in `config/config.yaml` and leave
`fastq_2` empty.

## Notes

- Differential analysis needs at least two conditions. Biological replicates
  are strongly recommended; set `strict_replicates: true` if you want the
  validator to fail early when a group has too few replicates.
- Large FASTQ, BAM, result, and reference files are ignored by git.
- See `docs/usage.md` for detailed setup and `docs/developer.md` for
  maintenance notes.

## License

MIT. See `LICENSE`.

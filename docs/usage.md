# Usage Guide

## 1. Prepare the Environment

```bash
conda env create -f envs/environment.yaml
conda activate rna_seq_pipeline
```

Run with Snakemake-managed environments:

```bash
snakemake --use-conda -n
snakemake --use-conda -c 8
```

The workflow uses per-rule environments in `envs/validate.yaml`, `envs/qc.yaml`,
`envs/align.yaml`, `envs/counts.yaml`, and `envs/deseq2.yaml`. The combined
`envs/environment.yaml` is useful for manual debugging.

## 2. Prepare Reference Files

Create a STAR index before alignment:

```bash
STAR \
  --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir reference/GRCh38/star_index \
  --genomeFastaFiles reference/GRCh38/GRCh38.p14.genome.fa \
  --sjdbGTFfile reference/GRCh38/gencode.v44.annotation.gtf \
  --sjdbOverhang 100
```

Then update `config/config.yaml`:

```yaml
reference:
  star_index: "reference/GRCh38/star_index"
  gtf: "reference/GRCh38/gencode.v44.annotation.gtf"
  genome_fasta: "reference/GRCh38/GRCh38.p14.genome.fa"
```

## 3. Prepare the Sample Sheet

Use a tab-delimited `config/samples.tsv`:

```text
sample	condition	fastq_1	fastq_2
control_1	control	data/control_1_R1.fastq.gz	data/control_1_R2.fastq.gz
control_2	control	data/control_2_R1.fastq.gz	data/control_2_R2.fastq.gz
treat_1	treat	data/treat_1_R1.fastq.gz	data/treat_1_R2.fastq.gz
treat_2	treat	data/treat_2_R1.fastq.gz	data/treat_2_R2.fastq.gz
```

Sample names may contain letters, numbers, `_`, `-`, and `.`.

## 4. Configure Analysis Parameters

Important options in `config/config.yaml`:

- `paired_end`: set to `true` or `false`.
- `trimming.enabled`: set to `true` to run fastp before STAR.
- `feature_counts_params.strand_specific`: use `0`, `1`, or `2`.
- `deseq2_params.control_group`: baseline condition.
- `deseq2_params.treat_group`: treatment condition.
- `deseq2_params.strict_replicates`: fail validation if replicate count is low.

## 5. Run the Workflow

Preview:

```bash
snakemake --use-conda -n -p
```

Run:

```bash
snakemake --use-conda -c 8 --latency-wait 60
```

Generate a report:

```bash
snakemake --report results/snakemake_report.html
```

## 6. Troubleshooting

- If validation fails, inspect `logs/validation/validate_inputs.log`.
- If STAR outputs zero input reads, verify FASTQ paths, compression suffixes,
  and read formatting.
- If featureCounts assigns zero reads, verify chromosome naming consistency
  between BAM and GTF, strandedness, and paired-end settings.
- If DESeq2 writes placeholder plots, inspect the message in the plot and
  confirm the count matrix has assigned reads and enough samples.

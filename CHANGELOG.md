# Changelog

## 2026-07-03

- Converted the pipeline to a sample-sheet-driven workflow.
- Added input/config validation before analysis steps run.
- Added paired-end aware FastQC, STAR, samtools, and featureCounts rules.
- Added optional fastp trimming.
- Reworked count merging to use explicit Snakemake inputs.
- Hardened DESeq2 analysis against empty/invalid count matrices and added normalized counts and PCA output.
- Added per-rule conda environment declarations for reproducible execution.
- Split conda environments by workflow step to make validation and debugging lighter.
- Updated documentation for usage, development, and reference management.

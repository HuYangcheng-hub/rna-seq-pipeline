# Developer Notes

## Design Principles

- Keep sample metadata in `config/samples.tsv`, not in the Snakefile.
- Fail early on invalid paths, references, and sample metadata.
- Use explicit Snakemake inputs instead of globbing runtime directories.
- Keep large data, results, and references out of git.
- Prefer small external scripts over long `run:` blocks.

## Important Files

- `Snakefile`: workflow graph and rule wiring.
- `config/config.yaml`: project-level parameters.
- `config/samples.tsv`: sample metadata and FASTQ paths.
- `scripts/validate_inputs.py`: early validation.
- `scripts/merge_counts.py`: featureCounts matrix merge.
- `scripts/deseq2_analysis.R`: DESeq2 statistics and plots.
- `envs/*.yaml`: per-step reproducible software environments.

## Validation Checklist

Run these before committing:

```bash
snakemake --use-conda -n -p
snakemake --lint text
python scripts/validate_inputs.py \
  --config config/config.yaml \
  --sample-sheet config/samples.tsv \
  --output /tmp/rna_seq_pipeline_samples.validated.tsv
```

For a clean end-to-end test, remove old outputs first:

```bash
rm -rf results logs .snakemake
snakemake --use-conda -c 8 --latency-wait 60
```

Only do this when you intentionally want to delete local outputs.

## Release Notes

Record user-facing workflow changes in `CHANGELOG.md`.

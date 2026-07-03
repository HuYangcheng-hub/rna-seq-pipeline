#!/usr/bin/env python3
"""Merge featureCounts output files into a gene-by-sample count matrix."""

from pathlib import Path

import pandas as pd


def sample_from_count_path(path):
    name = Path(path).name
    if name.endswith("_counts.txt"):
        return name[: -len("_counts.txt")]
    return Path(path).stem


def read_featurecounts(path, sample):
    df = pd.read_csv(path, sep="\t", comment="#")
    if "Geneid" not in df.columns:
        raise ValueError(f"{path} does not look like featureCounts output: missing Geneid")
    if len(df.columns) < 7:
        raise ValueError(f"{path} does not contain a count column")
    count_col = df.columns[-1]
    out = df.loc[:, ["Geneid", count_col]].copy()
    out.columns = ["gene_id", sample]
    out[sample] = pd.to_numeric(out[sample], errors="raise").astype(int)
    return out


def main():
    if "snakemake" in globals():
        count_files = list(snakemake.input.counts)
        output = snakemake.output.matrix
    else:
        count_files = sorted(str(path) for path in Path("results/counts").glob("*_counts.txt"))
        count_files = [path for path in count_files if "all_samples" not in path]
        output = "results/counts/all_samples_counts.txt"

    if not count_files:
        raise ValueError("No featureCounts files were provided")

    matrices = []
    for path in count_files:
        sample = sample_from_count_path(path)
        matrices.append(read_featurecounts(path, sample))

    merged = matrices[0]
    for matrix in matrices[1:]:
        merged = merged.merge(matrix, on="gene_id", how="outer")

    merged = merged.fillna(0)
    sample_cols = [col for col in merged.columns if col != "gene_id"]
    merged[sample_cols] = merged[sample_cols].astype(int)

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, sep="\t", index=False)

    print(f"Merged {len(count_files)} count files into {output_path}")
    print(f"Genes: {merged.shape[0]}")
    print(f"Samples: {', '.join(sample_cols)}")
    print(f"Total assigned counts: {int(merged[sample_cols].sum().sum())}")


if __name__ == "__main__":
    main()

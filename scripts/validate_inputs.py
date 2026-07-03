#!/usr/bin/env python3
"""Validate RNA-seq workflow configuration and sample sheet."""

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

import yaml


SAMPLE_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
FASTQ_SUFFIXES = (
    ".fastq",
    ".fq",
    ".fastq.gz",
    ".fq.gz",
)


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def warn(message):
    print(f"WARNING: {message}", file=sys.stderr)


def load_samples(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "condition", "fastq_1"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        rows = []
        for line_no, row in enumerate(reader, start=2):
            sample = row.get("sample", "").strip()
            condition = row.get("condition", "").strip()
            fastq_1 = row.get("fastq_1", "").strip()
            fastq_2 = row.get("fastq_2", "").strip()
            if not any([sample, condition, fastq_1, fastq_2]):
                continue
            if not sample:
                fail(f"Missing sample name at line {line_no}")
            if not SAMPLE_RE.match(sample):
                fail(f"Sample '{sample}' contains unsupported characters")
            if not condition:
                fail(f"Missing condition for sample '{sample}'")
            if not fastq_1:
                fail(f"Missing fastq_1 for sample '{sample}'")
            rows.append(
                {
                    "sample": sample,
                    "condition": condition,
                    "fastq_1": fastq_1,
                    "fastq_2": fastq_2,
                }
            )
    if not rows:
        fail(f"No samples found in {path}")
    return rows


def check_fastq(path, sample, label):
    fastq = Path(path)
    if not fastq.exists():
        fail(f"{label} for sample '{sample}' does not exist: {path}")
    if not fastq.is_file():
        fail(f"{label} for sample '{sample}' is not a file: {path}")
    if not str(fastq).endswith(FASTQ_SUFFIXES):
        warn(f"{label} for sample '{sample}' does not use a standard FASTQ suffix: {path}")
    if fastq.stat().st_size == 0:
        fail(f"{label} for sample '{sample}' is empty: {path}")


def check_reference(config):
    reference = config.get("reference", {})
    for key in ("star_index", "gtf", "genome_fasta"):
        value = reference.get(key)
        if not value:
            fail(f"Missing reference.{key} in config")
        path = Path(value)
        if not path.exists():
            fail(f"reference.{key} does not exist: {value}")
    index = Path(reference["star_index"])
    required_index_files = ["Genome", "SA", "SAindex", "genomeParameters.txt"]
    missing = [name for name in required_index_files if not (index / name).exists()]
    if missing:
        fail(
            "STAR index is incomplete; missing "
            + ", ".join(missing)
            + f" under {index}"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-sheet", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.config) as handle:
        config = yaml.safe_load(handle) or {}

    paired_end = bool(config.get("paired_end", True))
    rows = load_samples(args.sample_sheet)

    sample_counts = Counter(row["sample"] for row in rows)
    duplicates = sorted(sample for sample, count in sample_counts.items() if count > 1)
    if duplicates:
        fail(f"Duplicate sample names found: {', '.join(duplicates)}")

    for row in rows:
        check_fastq(row["fastq_1"], row["sample"], "fastq_1")
        if paired_end:
            if not row["fastq_2"]:
                fail(f"Missing fastq_2 for paired-end sample '{row['sample']}'")
            check_fastq(row["fastq_2"], row["sample"], "fastq_2")

    check_reference(config)

    condition_counts = Counter(row["condition"] for row in rows)
    if len(condition_counts) < 2:
        fail("At least two conditions are required for differential expression")

    de_params = config.get("deseq2_params", {})
    min_reps = int(de_params.get("min_replicates_per_group", 2))
    strict_reps = bool(de_params.get("strict_replicates", False))
    low_rep_groups = {
        condition: count
        for condition, count in condition_counts.items()
        if count < min_reps
    }
    if low_rep_groups:
        message = (
            "Some conditions have fewer than "
            f"{min_reps} replicates: "
            + ", ".join(f"{k}={v}" for k, v in sorted(low_rep_groups.items()))
        )
        if strict_reps:
            fail(message)
        warn(message)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["sample", "condition", "fastq_1", "fastq_2"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Validated {len(rows)} samples across {len(condition_counts)} conditions.")


if __name__ == "__main__":
    main()

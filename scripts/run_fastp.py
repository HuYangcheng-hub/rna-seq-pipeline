#!/usr/bin/env python3
"""Snakemake wrapper for paired-end or single-end fastp trimming."""

from snakemake.shell import shell


paired = bool(snakemake.config.get("paired_end", True))
extra = snakemake.config.get("trimming", {}).get("extra", "")

if paired:
    shell(
        "fastp "
        "-i {snakemake.input.r1} -I {snakemake.input.r2} "
        "-o {snakemake.output.r1} -O {snakemake.output.r2} "
        "--thread {snakemake.threads} "
        "--html {snakemake.output.html} --json {snakemake.output.json} "
        "{extra} > {snakemake.log} 2>&1"
    )
else:
    shell(
        "fastp "
        "-i {snakemake.input.r1} "
        "-o {snakemake.output.r1} "
        "--thread {snakemake.threads} "
        "--html {snakemake.output.html} --json {snakemake.output.json} "
        "{extra} > {snakemake.log} 2>&1"
    )

#!/usr/bin/env python3
"""Generate synthetic FASTQ reads for CI testing."""
import gzip
import os
import random
import subprocess
import sys

random.seed(42)
n_reads = 10000
samples = ["control_1", "control_2", "treat_1", "treat_2"]
os.makedirs("data", exist_ok=True)
ref = "reference/chr22/Homo_sapiens.GRCh38.dna.chromosome.22.fa"

# Try wgsim first (better reads), fall back to random
use_wgsim = os.path.exists(ref) and os.path.exists("/usr/bin/wgsim")
if use_wgsim:
    print("Using wgsim for read generation...")
    for sample in samples:
        subprocess.run([
            "wgsim", "-N", str(n_reads), "-1", "150", "-2", "150",
            ref,
            f"data/{sample}_R1.fastq", f"data/{sample}_R2.fastq"
        ], check=True)
        subprocess.run(["gzip", "-f", f"data/{sample}_R1.fastq", f"data/{sample}_R2.fastq"])
else:
    print("wgsim not found, generating random reads...")
    for sample in samples:
        with gzip.open(f"data/{sample}_R1.fastq.gz", "wt") as f1, \
             gzip.open(f"data/{sample}_R2.fastq.gz", "wt") as f2:
            for i in range(n_reads):
                seq1 = "".join(random.choices("ACGT", k=150))
                seq2 = "".join(random.choices("ACGT", k=150))
                qual = "I" * 150
                f1.write(f"@{sample}/{1}\n{seq1}\n+\n{qual}\n")
                f2.write(f"@{sample}/{1}\n{seq2}\n+\n{qual}\n")

print(f"Generated {len(samples)} x2 FASTQ files ({n_reads} reads each)")

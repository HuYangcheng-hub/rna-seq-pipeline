#!/usr/bin/env bash
# Download small chr22-enriched FASTQ files for end-to-end pipeline testing.
#
# This script downloads a subset of paired-end RNA-seq reads from a public
# dataset, sufficient to run the full pipeline.  Total download ~200 MB.
#
# Requirements: fasterq-dump (conda install -c bioconda sra-tools)
#               or wget + pigz
#
# Usage:
#   bash scripts/download_test_data.sh          # download to data/
#   bash scripts/download_test_data.sh /tmp/out  # custom output directory

set -euo pipefail

OUTDIR="${1:-data}"
mkdir -p "$OUTDIR"

# ── Method 1: fasterq-dump (preferred) ──────────────────────────────────
# Downloads a small number of reads from a public human RNA-seq run.
# SRR1039508: human cell line RNA-seq, paired-end, ~20M reads total.

if command -v fasterq-dump &>/dev/null; then
  echo "[download_test_data] Using fasterq-dump to fetch reads..."
  cd "$OUTDIR"

  # Fetch 300k spots (600k reads) – enough for a meaningful chr22 test
  fasterq-dump SRR1039508 \
    --max-spot 300000 \
    --split-files \
    --threads 4 \
    --progress \
    --outdir .

  # Compress and rename to match sample sheet expectations
  echo "[download_test_data] Compressing..."
  pigz -p 4 -f SRR1039508_1.fastq SRR1039508_2.fastq

  # Create symlinks for sample_1 (control rep 1)
  ln -sf SRR1039508_1.fastq.gz sample_1_R1.fastq.gz
  ln -sf SRR1039508_2.fastq.gz sample_1_R2.fastq.gz

  cd - >/dev/null
  echo "[download_test_data] Done: $OUTDIR/sample_1_R{1,2}.fastq.gz"

# ── Method 2: ENA direct download (fallback) ────────────────────────────
elif command -v wget &>/dev/null; then
  echo "[download_test_data] fasterq-dump not found; using wget to fetch from ENA..."
  BASE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508"

  # Download only the first 4M lines (~1M reads per file)
  # and compress on-the-fly to keep things small
  for READ in 1 2; do
    echo "[download_test_data] Fetching R${READ}..."
    wget -q -O - "${BASE}/SRR1039508_${READ}.fastq.gz" \
      | zcat \
      | head -n 4000000 \
      | pigz -p 4 -c \
      > "${OUTDIR}/sample_1_R${READ}.fastq.gz"
  done
  echo "[download_test_data] Done: $OUTDIR/sample_1_R{1,2}.fastq.gz"

else
  echo "[download_test_data] ERROR: neither fasterq-dump nor wget found."
  echo "  Install one of them first:"
  echo "    conda install -c bioconda sra-tools"
  echo "    or: brew install wget"
  exit 1
fi

# ── Verify ──────────────────────────────────────────────────────────────
echo ""
echo "[download_test_data] File sizes:"
ls -lh "$OUTDIR"/sample_1_R{1,2}.fastq.gz

echo ""
echo "Test data ready. To run a single-sample test, update config/samples.tsv"
echo "to reference the downloaded files and simplify to one condition pair."
echo ""
echo "See docs/usage.md for next steps."

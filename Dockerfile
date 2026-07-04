# ── RNA-seq Pipeline Docker Image ─────────────────────────────────────────
# Build:
#   docker build -t rna_seq_pipeline .
#
# Run interactively (mount data, reference, and results):
#   docker run -it --rm \
#     -v $(pwd)/data:/work/data \
#     -v $(pwd)/reference:/work/reference \
#     -v $(pwd)/results:/work/results \
#     rna_seq_pipeline
#
# Or use Singularity/Apptainer:
#   singularity build rna_seq_pipeline.sif docker-daemon://rna_seq_pipeline:latest
#   singularity exec --bind data,reference,results rna_seq_pipeline.sif snakemake -c 8
#
# Or let Snakemake pull containers per-rule via --use-singularity.

FROM mambaorg/micromamba:2.0-jammy

LABEL org.opencontainers.image.title="RNA-seq Pipeline"
LABEL org.opencontainers.image.description="End-to-end bulk RNA-seq workflow with Snakemake"
LABEL org.opencontainers.image.licenses="MIT"

WORKDIR /work

# Install the full environment
COPY envs/environment.yaml /tmp/environment.yaml
RUN micromamba install -y -n base -f /tmp/environment.yaml \
    && micromamba clean -afy

# Copy workflow files (excludes data/, reference/, results/ via .dockerignore)
COPY . /work

# Snakemake uses per-rule conda envs; the base env has core tooling.
# For container-only runs (no conda), set up a shared env:
#   snakemake --use-conda --conda-prefix /opt/conda-envs -c 8

ENTRYPOINT ["snakemake"]
CMD ["--use-conda", "-c", "8"]

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

get_input <- function(name, fallback) {
  if (exists("snakemake")) {
    value <- snakemake@input[[name]]
    if (!is.null(value)) return(as.character(value))
  }
  fallback
}

get_output <- function(name, fallback) {
  if (exists("snakemake")) {
    value <- snakemake@output[[name]]
    if (!is.null(value)) return(as.character(value))
  }
  fallback
}

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

plot_message <- function(path, title, message) {
  ensure_parent(path)
  p <- ggplot(data.frame(x = 0, y = 0, label = message), aes(x, y)) +
    geom_text(aes(label = label), size = 4) +
    theme_void() +
    labs(title = title) +
    xlim(-1, 1) +
    ylim(-1, 1)
  ggsave(path, p, width = 7, height = 5, dpi = 180)
}

write_empty_outputs <- function(reason, out_csv, out_sig, out_norm, out_volcano, out_heatmap, out_pca) {
  message(reason)
  empty_res <- data.frame(
    gene_id = character(),
    baseMean = numeric(),
    log2FoldChange = numeric(),
    lfcSE = numeric(),
    stat = numeric(),
    pvalue = numeric(),
    padj = numeric(),
    note = character()
  )
  ensure_parent(out_csv)
  write.csv(empty_res, out_csv, row.names = FALSE)
  write.csv(empty_res, out_sig, row.names = FALSE)
  write.csv(data.frame(gene_id = character()), out_norm, row.names = FALSE)
  plot_message(out_volcano, "Volcano Plot", reason)
  plot_message(out_heatmap, "Expression Heatmap", reason)
  plot_message(out_pca, "PCA Plot", reason)
}

config <- if (exists("snakemake")) snakemake@config else yaml::read_yaml("config/config.yaml")

count_file <- get_input("counts", "results/counts/all_samples_counts.txt")
sample_sheet <- get_input("sample_sheet", config$sample_sheet %||% "config/samples.tsv")

out_csv <- get_output("csv", "results/differential/deseq2_results.csv")
out_sig <- get_output("sig", "results/differential/deseq2_sig_genes.csv")
out_norm <- get_output("normalized", "results/differential/normalized_counts.csv")
out_volcano <- get_output("volcano", "results/differential/volcano_plot.png")
out_heatmap <- get_output("heatmap", "results/differential/heatmap.png")
out_pca <- get_output("pca", "results/differential/pca_plot.png")

for (path in c(out_csv, out_sig, out_norm, out_volcano, out_heatmap, out_pca)) {
  ensure_parent(path)
}

params <- config$deseq2_params %||% list()
padj_cutoff <- as.numeric(params$padj_cutoff %||% 0.05)
lfc_cutoff <- as.numeric(params$lfc_cutoff %||% 1.0)
control_group <- as.character(params$control_group %||% "control")
treat_group <- as.character(params$treat_group %||% "treat")
min_count <- as.numeric(params$min_count %||% 10)
min_samples <- as.numeric(params$min_samples %||% 2)

count_data <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
sample_data <- read.table(sample_sheet, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("sample", "condition")
missing_cols <- setdiff(required_cols, colnames(sample_data))
if (length(missing_cols) > 0) {
  stop("Sample sheet missing required columns: ", paste(missing_cols, collapse = ", "))
}

sample_data <- sample_data[sample_data$sample %in% colnames(count_data), , drop = FALSE]
sample_data <- sample_data[match(colnames(count_data), sample_data$sample), , drop = FALSE]
if (any(is.na(sample_data$sample))) {
  missing <- colnames(count_data)[is.na(sample_data$sample)]
  stop("Count matrix samples missing from sample sheet: ", paste(missing, collapse = ", "))
}

keep_groups <- sample_data$condition %in% c(control_group, treat_group)
sample_data <- sample_data[keep_groups, , drop = FALSE]
count_data <- count_data[, sample_data$sample, drop = FALSE]

if (ncol(count_data) < 2 || length(unique(sample_data$condition)) < 2) {
  write_empty_outputs(
    "DESeq2 skipped: at least two groups and two samples are required.",
    out_csv, out_sig, out_norm, out_volcano, out_heatmap, out_pca
  )
  quit(save = "no", status = 0)
}

count_data <- round(as.matrix(count_data))
storage.mode(count_data) <- "integer"

if (sum(count_data) == 0) {
  write_empty_outputs(
    "DESeq2 skipped: count matrix contains zero assigned reads.",
    out_csv, out_sig, out_norm, out_volcano, out_heatmap, out_pca
  )
  quit(save = "no", status = 0)
}

sample_data$condition <- factor(sample_data$condition, levels = c(control_group, treat_group))
rownames(sample_data) <- sample_data$sample

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_data,
  design = ~ condition
)

keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]

if (nrow(dds) == 0) {
  write_empty_outputs(
    "DESeq2 skipped: no genes passed the expression filter.",
    out_csv, out_sig, out_norm, out_volcano, out_heatmap, out_pca
  )
  quit(save = "no", status = 0)
}

dds <- tryCatch(
  DESeq(dds),
  error = function(e) {
    write_empty_outputs(
      paste("DESeq2 skipped:", conditionMessage(e)),
      out_csv, out_sig, out_norm, out_volcano, out_heatmap, out_pca
    )
    quit(save = "no", status = 0)
  }
)

res <- results(dds, contrast = c("condition", treat_group, control_group))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.csv(res_df, out_csv, row.names = FALSE)

normalized <- as.data.frame(counts(dds, normalized = TRUE))
normalized$gene_id <- rownames(normalized)
normalized <- normalized[, c("gene_id", setdiff(colnames(normalized), "gene_id"))]
write.csv(normalized, out_norm, row.names = FALSE)

sig_genes <- subset(
  res_df,
  !is.na(padj) & padj < padj_cutoff & !is.na(log2FoldChange) & abs(log2FoldChange) > lfc_cutoff
)
write.csv(sig_genes, out_sig, row.names = FALSE)

res_df$sig <- "NS"
res_df$sig[!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange > lfc_cutoff] <- "Up"
res_df$sig[!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff] <- "Down"
res_df$padj_plot <- pmax(res_df$padj, .Machine$double.xmin, na.rm = FALSE)

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj_plot), color = sig)) +
  geom_point(alpha = 0.65, size = 1.6, na.rm = TRUE) +
  scale_color_manual(values = c("NS" = "grey70", "Up" = "#D73027", "Down" = "#4575B4")) +
  theme_minimal(base_size = 12) +
  labs(
    title = paste(treat_group, "vs", control_group),
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Status"
  ) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", alpha = 0.5)
ggsave(out_volcano, volcano, width = 8, height = 6, dpi = 300)

vsd <- tryCatch(vst(dds, blind = FALSE), error = function(e) NULL)
if (is.null(vsd) || nrow(vsd) < 2 || ncol(vsd) < 2) {
  plot_message(out_heatmap, "Expression Heatmap", "Not enough data for heatmap.")
  plot_message(out_pca, "PCA Plot", "Not enough data for PCA.")
} else {
  ranked <- res_df[!is.na(res_df$padj), , drop = FALSE]
  if (nrow(sig_genes) >= 2) {
    top_genes <- head(sig_genes$gene_id, 30)
  } else if (nrow(ranked) >= 2) {
    top_genes <- head(ranked$gene_id, 30)
  } else {
    top_genes <- head(rownames(vsd), min(30, nrow(vsd)))
  }
  top_genes <- intersect(top_genes, rownames(vsd))

  if (length(top_genes) < 2) {
    plot_message(out_heatmap, "Expression Heatmap", "Not enough genes for heatmap.")
  } else {
    mat <- assay(vsd)[top_genes, , drop = FALSE]
    mat <- mat - rowMeans(mat)
    annotation_col <- data.frame(condition = sample_data$condition)
    rownames(annotation_col) <- sample_data$sample
    pheatmap(
      mat,
      annotation_col = annotation_col,
      main = "Top variable / differential genes",
      filename = out_heatmap,
      width = 8,
      height = 10
    )
  }

  pca <- prcomp(t(assay(vsd)))
  percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  pca_df <- data.frame(
    sample = rownames(pca$x),
    condition = sample_data[rownames(pca$x), "condition"],
    PC1 = pca$x[, 1],
    PC2 = if (ncol(pca$x) >= 2) pca$x[, 2] else 0
  )
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, show.legend = FALSE) +
    theme_minimal(base_size = 12) +
    labs(
      title = "PCA of variance-stabilized counts",
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[min(2, length(percent_var))], "% variance")
    )
  ggsave(out_pca, p, width = 7, height = 5, dpi = 300)
}

message("DESeq2 analysis completed.")
message("Total genes tested: ", nrow(res_df))
message("Significant genes: ", nrow(sig_genes))

# ── Reproducibility ──────────────────────────────────────────────────────
session_file <- file.path(dirname(out_csv), "session_info.txt")
writeLines(capture.output(sessionInfo()), session_file)
message("Session info written to ", session_file)

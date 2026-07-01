#!/usr/bin/env Rscript
#
# DESeq2 差异表达分析
# 输入: all_samples_counts.txt
# 输出: 差异基因列表、火山图、热图
#

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)

# ---- 读取配置文件 ----
library(yaml)
config <- yaml::yaml.load_file("config/config.yaml")
padj_cutoff <- config$deseq2_params$padj_cutoff
lfc_cutoff <- config$deseq2_params$lfc_cutoff
control_group <- config$deseq2_params$control_group
treat_group <- config$deseq2_params$treat_group

# ---- 读取表达矩阵 ----
count_data <- read.table("results/counts/all_samples_counts.txt",
                         header = TRUE, row.names = 1, sep = "\t")

# ---- 构建样本信息 ----
sample_names <- colnames(count_data)
condition <- ifelse(grepl(control_group, sample_names), "control", "treat")
condition <- factor(condition, levels = c("control", "treat"))

col_data <- data.frame(row.names = sample_names, condition = condition)

# ---- DESeq2 分析 ----
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ condition
)

# 过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treat", "control"))
res_df <- as.data.frame(res)

# ---- 保存结果 ----
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, "results/differential/deseq2_results.csv")

# 筛选显著差异基因
sig_genes <- subset(res_df, padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)
write.csv(sig_genes, "results/differential/deseq2_sig_genes.csv")

# ---- 火山图 ----
res_df$sig <- "NS"
res_df$sig[res_df$padj < padj_cutoff & res_df$log2FoldChange > lfc_cutoff] <- "Up"
res_df$sig[res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff] <- "Down"

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("NS" = "grey", "Up" = "red", "Down" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot (DEGs)",
       x = "log2(Fold Change)", y = "-log10(adjusted p-value)") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", alpha = 0.5)

ggsave("results/differential/volcano_plot.png", p, width = 8, height = 6, dpi = 300)

# ---- 热图 ----
if (nrow(sig_genes) > 2) {
  top_genes <- head(rownames(sig_genes), 30)
  mat <- assay(rlog(dds))[top_genes, ]
  mat <- mat - rowMeans(mat)

  annotation_col <- data.frame(condition = condition)
  rownames(annotation_col) <- colnames(mat)

  pheatmap(mat,
           annotation_col = annotation_col,
           main = "Top 30 DEGs Heatmap",
           filename = "results/differential/heatmap.png",
           width = 8, height = 10, dpi = 300)
}

cat("DESeq2 analysis completed!\\n")
cat(paste("  Total significant DEGs:", nrow(sig_genes), "\\n"))

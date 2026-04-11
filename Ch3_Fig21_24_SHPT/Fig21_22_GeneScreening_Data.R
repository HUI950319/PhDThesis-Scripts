# ============================================================
# CH3_Fig21_22_Data.R
# SHPT 候选基因三维度筛选: 轨迹相关 + 差异表达 + 比例相关
# 比例相关: Bulk 用反卷积 GEP 提取 Cluster 2 亚群表达 z-score
#          scRNA 用样本级均值 z-score
# 输出: ./out/gene_screening/res_merged_SH.csv
#       ./out/gene_screening/venn_gene_lists_SH.Rdata
# ============================================================

library(scMMR)
library(Seurat)
library(qs)
library(tidyverse)

use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

out_dir    <- "./out/gene_screening"
bulk_dir   <- "~/project/bulk/counts"
deconv_dir <- "./out/deconv"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cluster2_col <- "Cluster 2"

# --- 加载数据 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime_V3.qs")

# ============================================================
# 1. RunTraceGene: 轨迹相关基因 (SHPT → Lineage2)
# ============================================================

seu_SH <- subset(seu, group == "SHPT")
res_SH <- scMMR::RunTraceGene(seu_SH, lineages = c("Lineage2"))
trace_genes_SH <- unique(res_SH$gene)

# ============================================================
# 2. RunDE: SHPT vs Normal (全基因, 无过滤)
# ============================================================

seu_SH_Normal <- subset(seu, group %in% c("SHPT", "Normal"))
markers_SH <- scMMR::RunDE(
  seu_SH_Normal,
  group.by     = "group",
  group1       = "SHPT",
  group2       = "Normal",
  only.pos     = FALSE,
  fc.threshold = 1,
  min.pct      = 0
)

# ============================================================
# 3. 比例相关: Cluster 2 亚群表达 z-score ~ Cluster 2 proportion
#    Bulk: 反卷积 adaptive GEP 提取亚群特异性表达
#    scRNA: 样本级均值
# ============================================================

# --- 3a. 加载反卷积模型 & Bulk counts ---
seu_ref    <- qs::qread("./out/seu.qs")
model_path <- file.path(deconv_dir, "deconv_model.pt")

read_featurecounts <- function(file_path) {
  raw <- read.delim(file_path, header = TRUE, row.names = 1,
                    check.names = FALSE, comment.char = "#")
  if ("gene_name" %in% colnames(raw)) {
    gn <- raw$gene_name; raw$gene_name <- NULL
    keep <- !duplicated(gn)
    raw <- raw[keep, , drop = FALSE]; rownames(raw) <- gn[keep]
  }
  counts <- as.matrix(raw[, sapply(raw, is.numeric), drop = FALSE])
  colnames(counts) <- gsub(".*/(\\w+)_Aligned.*", "\\1", colnames(counts))
  counts
}

counts_shpt <- read_featurecounts(file.path(bulk_dir, "self_SHPT", "self_SHPT_counts.txt"))

counts_norm <- read_featurecounts(file.path(bulk_dir, "self_Normal", "self_Normal_counts.txt"))
counts_bulk <- cbind(counts_shpt, counts_norm)

# --- 3b. Bulk: adaptive 反卷积 → Cluster 2 亚群特异性表达 z-score ---
deconv_res <- scMMR::DNN_deconv_predict(
  counts_bulk, model_path, device = "auto",
  adaptive = TRUE, adaptive_mode = "high-resolution"
)
props_bulk  <- deconv_res$proportions       # data.frame: samples × cell_types
sig_bulk    <- deconv_res$sigmatrix          # 3D array: samples × cell_types × genes

genes_use <- dimnames(sig_bulk)[[3]]   # 全基因: sigmatrix 中所有基因
cl2_expr  <- sig_bulk[, cluster2_col, genes_use, drop = FALSE]
cl2_expr  <- matrix(cl2_expr, nrow = dim(sig_bulk)[1],
                    dimnames = list(dimnames(sig_bulk)[[1]], genes_use))
gene_z_bulk <- as.data.frame(scale(cl2_expr), check.names = FALSE)

# --- 3c. scRNA: Cluster 2 样本级均值 z-score ---
seu_cl2 <- subset(seu, annotation_final == cluster2_col)
sc_mat  <- AverageExpression(
  seu_cl2, features = genes_use, group.by = "sample_id",
  assays = "RNA", slot = "data"
)[["RNA"]]
gene_z_sc <- as.data.frame(t(scale(t(sc_mat))), check.names = FALSE)

# --- 3d. 组装 Spearman 相关: gene z-score ~ Cluster 2 proportion ---
# Bulk 比例
prop_cl2_bulk <- props_bulk[[cluster2_col]]
names(prop_cl2_bulk) <- rownames(props_bulk)

# scRNA 比例 (prop.table)
ct_tab <- table(seu$sample_id, seu$annotation_final)
prop_sc <- as.data.frame.matrix(prop.table(ct_tab, margin = 1))
prop_cl2_sc <- prop_sc[[cluster2_col]]
names(prop_cl2_sc) <- rownames(prop_sc)

# Spearman 相关 (Bulk + scRNA 合并)
cor_list <- lapply(genes_use, function(g) {
  # Bulk
  bulk_samples <- intersect(names(prop_cl2_bulk), rownames(gene_z_bulk))
  cor_b <- cor.test(gene_z_bulk[bulk_samples, g],
                    prop_cl2_bulk[bulk_samples], method = "spearman")
  # scRNA
  sc_samples <- intersect(names(prop_cl2_sc), colnames(gene_z_sc))
  cor_s <- cor.test(as.numeric(gene_z_sc[g, sc_samples]),
                    prop_cl2_sc[sc_samples], method = "spearman")
  data.frame(
    gene       = g,
    prop_cor   = mean(c(cor_b$estimate, cor_s$estimate)),
    prop_pval  = max(cor_b$p.value, cor_s$p.value),
    bulk_rho   = cor_b$estimate,
    sc_rho     = cor_s$estimate,
    stringsAsFactors = FALSE
  )
})
prop_cor_df_SH <- dplyr::bind_rows(cor_list)

# ============================================================
# 4. 三维度合并 + Composite Score
# ============================================================

df_trace <- res_SH %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    trace_cor  = mean(cor, na.rm = TRUE),
    trace_pval = min(pval, na.rm = TRUE),
    .groups    = "drop"
  )

df_de <- markers_SH %>%
  dplyr::select(gene, de_logFC = avg_log2FC, de_pval = p_val_adj)

df_merged <- df_trace %>%
  dplyr::full_join(df_de, by = "gene") %>%
  dplyr::full_join(prop_cor_df_SH, by = "gene") %>%
  dplyr::mutate(
    score = abs(trace_cor) + abs(de_logFC) + abs(prop_cor)
  ) %>%
  dplyr::arrange(dplyr::desc(score))

# ============================================================
# 5. Venn 基因集: 三维度各自显著基因
# ============================================================

genes_trace <- df_merged %>%
  dplyr::filter(trace_pval < 0.05) %>% dplyr::pull(gene)
genes_de    <- df_merged %>%
  dplyr::filter(de_pval < 0.05, abs(de_logFC) > 0.5) %>% dplyr::pull(gene)
genes_prop  <- df_merged %>%
  dplyr::filter(prop_pval < 0.05) %>% dplyr::pull(gene)

venn_gene_lists_SH <- list(
  "Trajectory-associated genes\n(padj < 0.05)"        = genes_trace,
  "Differentially expressed genes\n(log2FC > 0.5, padj < 0.05)" = genes_de,
  "Proportion-correlated genes\n(Spearman, padj < 0.05)"        = genes_prop
)

# ============================================================
# 6. 保存输出
# ============================================================

write.csv(df_merged, file.path(out_dir, "res_merged_SH.csv"),
          row.names = FALSE)
save(venn_gene_lists_SH, file = file.path(out_dir, "venn_gene_lists_SH.Rdata"))

# 同时保存 gene z-score 及 proportion (供 Fig22 Gene Validation Plot 使用)
save(gene_z_bulk, gene_z_sc, props_bulk, prop_cl2_sc,
     file = file.path(out_dir, "gene_zscore_prop_SH.Rdata"))


# ============================================================
# CH4_Fig25_26_Data.R
# PHPT 候选基因三维度筛选: 轨迹相关 + 差异表达 + 比例相关
# 比例相关: Bulk 用反卷积 GEP 提取 Cluster 3 亚群表达 z-score
#          scRNA 用样本级均值 z-score
# 含外部验证数据集: PRJNA516535 (Bulk), GSE190773 (scRNA)
# 输出: ./out/gene_screening/res_merged_PH.csv
#       ./out/gene_screening/venn_gene_lists_PH.Rdata
#       ./out/gene_screening/gene_zscore_prop_PH.Rdata
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

cluster3_col <- "Cluster 3"

# --- 加载数据 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime_V3.qs")

# ============================================================
# 1. RunTraceGene: 轨迹相关基因 (PHPT → Lineage3)
# ============================================================

seu_PH <- subset(seu, group == "PHPT")
res_PH <- scMMR::RunTraceGene(seu_PH, lineages = c("Lineage3"))
trace_genes_PH <- unique(res_PH$gene)

# ============================================================
# 2. RunDE: PHPT vs Normal (全基因, 无过滤)
# ============================================================

seu_PH_Normal <- subset(seu, group %in% c("PHPT", "Normal"))
markers_PH <- scMMR::RunDE(
  seu_PH_Normal,
  group.by     = "group",
  group1       = "PHPT",
  group2       = "Normal",
  only.pos     = FALSE,
  fc.threshold = 1,
  min.pct      = 0
)

# ============================================================
# 3. 比例相关: Cluster 3 亚群表达 z-score ~ Cluster 3 proportion
#    Bulk: 反卷积 adaptive GEP 提取亚群特异性表达
#    scRNA: 样本级均值
#    含外部数据集: PRJNA516535 (Bulk), GSE190773 (scRNA)
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

# Our Bulk: PHPT + Normal
counts_phpt <- read_featurecounts(file.path(bulk_dir, "self_PHPT", "self_PHPT_counts.txt"))
counts_norm <- read_featurecounts(file.path(bulk_dir, "self_Normal", "self_Normal_counts.txt"))
counts_our  <- cbind(counts_phpt, counts_norm)

# External Bulk: PRJNA516535
counts_ext <- read_featurecounts(file.path(bulk_dir, "PRJNA516535", "PRJNA516535_counts.txt"))

# --- 3b. Bulk: adaptive 反卷积 → Cluster 3 亚群特异性表达 z-score ---
# Our Bulk
deconv_our <- scMMR::DNN_deconv_predict(
  counts_our, model_path, device = "auto",
  adaptive = TRUE, adaptive_mode = "high-resolution"
)
props_our  <- deconv_our$proportions
sig_our    <- deconv_our$sigmatrix

# External Bulk (PRJNA516535)
deconv_ext <- scMMR::DNN_deconv_predict(
  counts_ext, model_path, device = "auto",
  adaptive = TRUE, adaptive_mode = "high-resolution"
)
props_ext  <- deconv_ext$proportions
sig_ext    <- deconv_ext$sigmatrix

genes_use <- dimnames(sig_our)[[3]]   # 全基因: sigmatrix 中所有基因

# Our Bulk → Cluster 3 expression z-score
cl3_our  <- sig_our[, cluster3_col, genes_use, drop = FALSE]
cl3_our  <- matrix(cl3_our, nrow = dim(sig_our)[1],
                   dimnames = list(dimnames(sig_our)[[1]], genes_use))
gene_z_bulk_our <- as.data.frame(scale(cl3_our), check.names = FALSE)

# External Bulk → Cluster 3 expression z-score
genes_ext <- intersect(genes_use, dimnames(sig_ext)[[3]])
cl3_ext   <- sig_ext[, cluster3_col, genes_ext, drop = FALSE]
cl3_ext   <- matrix(cl3_ext, nrow = dim(sig_ext)[1],
                    dimnames = list(dimnames(sig_ext)[[1]], genes_ext))
gene_z_bulk_ext <- as.data.frame(scale(cl3_ext), check.names = FALSE)

# --- 3c. scRNA: Cluster 3 样本级均值 z-score ---
# Our scRNA
seu_cl2 <- subset(seu, annotation_final == cluster3_col)
sc_mat_our <- AverageExpression(
  subset(seu_cl2, group %in% c("PHPT", "Normal")),
  features = genes_use, group.by = "sample_id",
  assays = "RNA", slot = "data"
)[["RNA"]]
gene_z_sc_our <- as.data.frame(t(scale(t(sc_mat_our))), check.names = FALSE)

# External scRNA (GSE190773)
seu_ext <- qs::qread("./out/external/GSE190773_seu.qs")
seu_ext_cl2 <- subset(seu_ext, annotation_final == cluster3_col)
sc_mat_ext <- AverageExpression(
  seu_ext_cl2, features = intersect(genes_use, rownames(seu_ext)),
  group.by = "sample_id", assays = "RNA", slot = "data"
)[["RNA"]]
gene_z_sc_ext <- as.data.frame(t(scale(t(sc_mat_ext))), check.names = FALSE)

# --- 3d. 组装 Spearman 相关: gene z-score ~ Cluster 3 proportion ---
# Bulk 比例 (our)
prop_cl2_our <- props_our[[cluster3_col]]
names(prop_cl2_our) <- rownames(props_our)

# Bulk 比例 (external)
prop_cl2_ext <- props_ext[[cluster3_col]]
names(prop_cl2_ext) <- rownames(props_ext)

# scRNA 比例 (our)
ct_tab <- table(seu$sample_id, seu$annotation_final)
prop_sc <- as.data.frame.matrix(prop.table(ct_tab, margin = 1))
prop_cl2_sc_our <- prop_sc[[cluster3_col]]
names(prop_cl2_sc_our) <- rownames(prop_sc)

# scRNA 比例 (external)
ct_tab_ext <- table(seu_ext$sample_id, seu_ext$annotation_final)
prop_sc_ext <- as.data.frame.matrix(prop.table(ct_tab_ext, margin = 1))
prop_cl2_sc_ext <- prop_sc_ext[[cluster3_col]]
names(prop_cl2_sc_ext) <- rownames(prop_sc_ext)

# Spearman 相关 (4 来源合并: our Bulk + our scRNA + ext Bulk + ext scRNA)
cor_list <- lapply(genes_use, function(g) {
  rhos <- c(); pvals <- c()
  # Our Bulk
  bs <- intersect(names(prop_cl2_our), rownames(gene_z_bulk_our))
  if (length(bs) >= 3 && g %in% colnames(gene_z_bulk_our)) {
    ct <- cor.test(gene_z_bulk_our[bs, g], prop_cl2_our[bs], method = "spearman")
    rhos <- c(rhos, ct$estimate); pvals <- c(pvals, ct$p.value)
  }
  # Our scRNA
  ss <- intersect(names(prop_cl2_sc_our), colnames(gene_z_sc_our))
  if (length(ss) >= 3 && g %in% rownames(gene_z_sc_our)) {
    ct <- cor.test(as.numeric(gene_z_sc_our[g, ss]), prop_cl2_sc_our[ss], method = "spearman")
    rhos <- c(rhos, ct$estimate); pvals <- c(pvals, ct$p.value)
  }
  # External Bulk (PRJNA516535)
  bs2 <- intersect(names(prop_cl2_ext), rownames(gene_z_bulk_ext))
  if (length(bs2) >= 3 && g %in% colnames(gene_z_bulk_ext)) {
    ct <- cor.test(gene_z_bulk_ext[bs2, g], prop_cl2_ext[bs2], method = "spearman")
    rhos <- c(rhos, ct$estimate); pvals <- c(pvals, ct$p.value)
  }
  # External scRNA (GSE190773)
  ss2 <- intersect(names(prop_cl2_sc_ext), colnames(gene_z_sc_ext))
  if (length(ss2) >= 3 && g %in% rownames(gene_z_sc_ext)) {
    ct <- cor.test(as.numeric(gene_z_sc_ext[g, ss2]), prop_cl2_sc_ext[ss2], method = "spearman")
    rhos <- c(rhos, ct$estimate); pvals <- c(pvals, ct$p.value)
  }
  if (length(rhos) == 0) return(NULL)
  data.frame(gene = g, prop_cor = mean(rhos), prop_pval = max(pvals),
             stringsAsFactors = FALSE)
})
prop_cor_df_PH <- dplyr::bind_rows(cor_list)

# ============================================================
# 4. 三维度合并 + Composite Score
# ============================================================

df_trace <- res_PH %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    trace_cor  = mean(cor, na.rm = TRUE),
    trace_pval = min(pval, na.rm = TRUE),
    .groups    = "drop"
  )

df_de <- markers_PH %>%
  dplyr::select(gene, de_logFC = avg_log2FC, de_pval = p_val_adj)

df_merged <- df_trace %>%
  dplyr::full_join(df_de, by = "gene") %>%
  dplyr::full_join(prop_cor_df_PH, by = "gene") %>%
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

venn_gene_lists_PH <- list(
  "Trajectory-associated genes\n(padj < 0.05)"        = genes_trace,
  "Differentially expressed genes\n(log2FC > 0.5, padj < 0.05)" = genes_de,
  "Proportion-correlated genes\n(Spearman, padj < 0.05)"        = genes_prop
)

# ============================================================
# 6. 保存
# ============================================================

write.csv(df_merged, file.path(out_dir, "res_merged_PH.csv"),
          row.names = FALSE)
save(venn_gene_lists_PH, file = file.path(out_dir, "venn_gene_lists_PH.Rdata"))

# 保存 gene z-score 及 proportion (供 Fig26 Gene Validation Plot 使用)
# 4 来源: our Bulk / our scRNA / ext Bulk (PRJNA516535) / ext scRNA (GSE190773)
save(gene_z_bulk_our, gene_z_sc_our, props_our, prop_cl2_sc_our,
     gene_z_bulk_ext, gene_z_sc_ext, props_ext, prop_cl2_sc_ext,
     file = file.path(out_dir, "gene_zscore_prop_PH.Rdata"))


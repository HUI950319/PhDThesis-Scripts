# ==============================================================================
# Figure 20 Data: BayesPrism 解卷积验证 — 与 scMMR DNN 比例的一致性
# ==============================================================================
#
# 描述：
#   使用 BayesPrism 对 Bulk RNA-seq 独立解卷积，验证 scMMR DNN 结果。
#
# 输入：
#   - ./out/seu.qs (scRNA-seq 参考)
#   - Bulk counts: self_PHPT (n=12), self_SHPT (n=12), PRJNA516535 (n=10)
#   - ./out/Figures/Figure19/plot_data.Rdata (scMMR 解卷积的 df1-df6)
#
# 输出：
#   - ./out/deconv/bp_*.rds (BayesPrism 完整结果)
#   - ./out/Figures/Figure20/plot_data2.Rdata (df_PH_bp, df_SH_bp)
#
# 依赖：BayesPrism, Seurat, qs, tidyverse
# ==============================================================================

library(BayesPrism)
library(Seurat)
library(qs)
library(tidyverse)

bulk_dir <- "/home/oyh/project/bulk/counts"
out_dir  <- "./out/deconv"
fig_dir  <- "./out/Figures/Figure20"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. 构建 scRNA-seq 参考矩阵
# ==============================================================================

seu_ref <- qs::qread("./out/seu.qs")

# BayesPrism 要求：细胞 x 基因
sc_counts         <- t(as.matrix(Seurat::GetAssayData(seu_ref, assay = "RNA", layer = "counts")))
cell_type_labels  <- seu_ref$celltype
cell_state_labels <- seu_ref$celltype

# ==============================================================================
# 2. 读取 Bulk counts
# ==============================================================================

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
  t(counts)  # 转置为 样本 x 基因
}

bulk_phpt  <- read_featurecounts(file.path(bulk_dir, "self_PHPT",   "self_PHPT_counts.txt"))
bulk_shpt  <- read_featurecounts(file.path(bulk_dir, "self_SHPT",   "self_SHPT_counts.txt"))
bulk_prjna <- read_featurecounts(file.path(bulk_dir, "PRJNA516535", "PRJNA516535_counts.txt"))

# ==============================================================================
# 3. BayesPrism 解卷积
# ==============================================================================

run_bayesprism <- function(bulk_sxg, sc_sxg, cell_types, cell_states,
                           save_path = NULL) {
  bp_obj <- new.prism(
    reference = sc_sxg, mixture = bulk_sxg,
    input.type = "count.matrix",
    cell.type.labels = cell_types, cell.state.labels = cell_states,
    key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1
  )
  bp_res <- run.prism(prism = bp_obj, n.cores = 10)
  theta  <- get.fraction(bp_res, which = "final", state.or.type = "type")

  if (!is.null(save_path)) saveRDS(bp_res, save_path)

  props <- as.data.frame(theta)
  props$sample_id <- rownames(props)
  props[, c("sample_id", setdiff(colnames(props), "sample_id"))]
}

bp_props_phpt <- run_bayesprism(
  bulk_phpt, sc_counts, cell_type_labels, cell_state_labels,
  save_path = file.path(out_dir, "bp_self_PHPT.rds"))

bp_props_shpt <- run_bayesprism(
  bulk_shpt, sc_counts, cell_type_labels, cell_state_labels,
  save_path = file.path(out_dir, "bp_self_SHPT.rds"))

bp_props_prjna <- run_bayesprism(
  bulk_prjna, sc_counts, cell_type_labels, cell_state_labels,
  save_path = file.path(out_dir, "bp_PRJNA516535.rds"))

# ==============================================================================
# 4. 合并 scMMR 比例 & BayesPrism 比例
# ==============================================================================

load("./out/Figures/Figure19/plot_data.Rdata")

cluster3_col <- "Cluster 3"
cluster2_col <- "Cluster 2"

df1$bp_pct <- bp_props_prjna[[cluster3_col]]
df3$bp_pct <- bp_props_phpt[[cluster3_col]]
df5$bp_pct <- bp_props_shpt[[cluster2_col]]

# ==============================================================================
# 5. 组装 & 保存
# ==============================================================================

df_PH_bp <- dplyr::bind_rows(
  df1 %>% dplyr::mutate(group = "Bulk (PRJNA516535, n = 10)"),
  df3 %>% dplyr::mutate(group = "Bulk (our, n = 12)")
) %>%
  dplyr::mutate(group = factor(group, levels = c("Bulk (our, n = 12)",
                                          "Bulk (PRJNA516535, n = 10)")))

df_SH_bp <- df5 %>%
  dplyr::mutate(group = factor("Bulk (our, n = 12)"))

save(df_PH_bp, df_SH_bp,
     file = file.path(fig_dir, "plot_data2.Rdata"))

save(bp_props_phpt, bp_props_shpt, bp_props_prjna,
     file = file.path(out_dir, "all_bayesprism_proportions.Rdata"))
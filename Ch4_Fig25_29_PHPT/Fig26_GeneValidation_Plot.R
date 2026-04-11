# ============================================================
# CH4_Fig26_Plot.R
# PHPT 候选基因验证: 轨迹相关 + 比例相关散点 + 差异表达小提琴
# 输入: ./out/gene_screening/gene_zscore_prop_PH.Rdata
#       ./out/trajectory/seu_with_pseudotime_V3.qs
# 输出: ./figures/CH4_Fig26_PHPT_Gene_Validation.pdf
# ============================================================

library(scMMR)
library(UtilsR)
library(Seurat)
library(qs)
library(tidyverse)
library(patchwork)

out_dir      <- "./out/gene_screening"
fig_dir      <- "./figures"
cluster3_col <- "Cluster 3"
top_genes    <- c("CHGB", "TNFRSF12A", "NCAM2")

pal_lancet <- UtilsR::get_palette("lancet")

# --- 加载数据 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime_V3.qs")
seu_PH_Normal <- subset(seu, group %in% c("PHPT", "Normal"))
seu_PH_Normal$celltype <- factor(seu_PH_Normal$celltype,
                                 levels = paste("cluster", c(0, 1, 3)))
seu_PH_Normal$group <- factor(seu_PH_Normal$group,
                              levels = c("PHPT", "Normal"))

load(file.path(out_dir, "gene_zscore_prop_PH.Rdata"))
# → gene_z_bulk_our, gene_z_sc_our, props_our, prop_cl2_sc_our
# → gene_z_bulk_ext, gene_z_sc_ext, props_ext, prop_cl2_sc_ext

# ============================================================
# Row 1: PlotDynamicFeatures — 轨迹相关
# ============================================================

p2 <- scMMR::PlotDynamicFeatures(
  subset(seu_PH_Normal, group == "PHPT"),
  pseudotime     = "Lineage3",
  features       = top_genes,
  group.by       = "celltype",
  smooth_k       = 2,
  pt.size        = 0.6,
  layer          = "data",
  exp_method     = "raw",
  stat_method    = "spearman",
  line_palcolor  = pal_lancet[2],
  point_palcolor = pal_lancet,
  raster         = TRUE,
  ncol           = 3
)

# ============================================================
# Row 2: PlotScatter — 比例相关 (4 来源反卷积 z-score)
# ============================================================

# --- 构建 4 来源 scatter df ---
build_bulk_df <- function(gene_z, props, genes, label) {
  df <- data.frame(
    sample      = rownames(gene_z),
    cluster_pct = props[[cluster3_col]],
    gene_z[, genes, drop = FALSE],
    check.names = FALSE
  )
  colnames(df)[3:ncol(df)] <- paste0("gene_", genes)
  df$group <- label
  df
}

build_sc_df <- function(gene_z, prop_cl2, genes, label) {
  ss <- intersect(names(prop_cl2), colnames(gene_z))
  df <- data.frame(
    sample      = ss,
    cluster_pct = prop_cl2[ss],
    t(gene_z[genes, ss, drop = FALSE]),
    check.names = FALSE
  )
  colnames(df)[3:ncol(df)] <- paste0("gene_", genes)
  df$group <- label
  df
}

df_PH <- dplyr::bind_rows(
  build_sc_df(gene_z_sc_our, prop_cl2_sc_our, top_genes, "SC"),
  build_bulk_df(gene_z_bulk_our, props_our, top_genes, "Bulk"),
  build_sc_df(gene_z_sc_ext, prop_cl2_sc_ext, top_genes, "GSE190773"),
  build_bulk_df(gene_z_bulk_ext, props_ext, top_genes, "PRJNA516535")
)

n_sc   <- sum(df_PH$group == "SC")
n_bulk <- sum(df_PH$group == "Bulk")
n_ext1 <- sum(df_PH$group == "GSE190773")
n_ext2 <- sum(df_PH$group == "PRJNA516535")

df_PH <- df_PH %>%
  dplyr::mutate(group = factor(
    group,
    levels = c("SC", "Bulk", "GSE190773", "PRJNA516535"),
    labels = c(
      paste0("scRNA (our, n = ", n_sc, ")"),
      paste0("Bulk (our, n = ", n_bulk, ")"),
      paste0("scRNA (GSE190773, n = ", n_ext1, ")"),
      paste0("Bulk (PRJNA516535, n = ", n_ext2, ")")
    )
  ))

p3 <- scMMR::PlotScatter(
  df_PH, "cluster_pct",
  paste0("gene_", top_genes),
  group.by     = "group",
  point.size   = 1.5,
  point.alpha  = 1,
  smooth.size  = 1,
  smooth.color = "#925E9F",
  palette      = "lancet",
  global.cor   = TRUE
)

# ============================================================
# Row 3: FeaturePlot3 — 差异表达小提琴
# ============================================================

p4 <- scMMR::FeaturePlot3(
  seu_PH_Normal,
  features   = top_genes,
  group.by   = "group",
  split.by   = "celltype",
  palcolor   = c("#00468BFF", "#ED0000FF", "#0099B4FF"),
  type       = "violin",
  add_box    = TRUE,
  ylab       = "Expression level (log-normalized)",
  stack      = FALSE
)

# ============================================================
# 组合三行面板
# ============================================================

P_row1 <- p2 %>%
  fmt_text(ylab = c("Normalized expression", "", ""), title = "") %>%
  fmt_legend(collect = TRUE, legend.position = "right",
             legend_theme = leg1()) %>%
  fmt_strip(label = paste(top_genes, "(trajectory correlation)"),
            label_color = "black") %>%
  fmt_panel()

P_row2 <- p3 %>%
  fmt_legend(collect = TRUE, legend.position = "right",
             legend_theme = leg1()) %>%
  fmt_text(xlab = "Percentage of Cluster 3",
           ylab = c("Expression (zscore)", "", ""),
           title = "") %>%
  fmt_strip(label = paste(top_genes, "(proportion correlation)"),
            label_color = "black", label_fill = "grey75") %>%
  fmt_panel()

P_row3 <- p4 %>%
  fmt_legend(collect = TRUE, legend_theme = leg1()) %>%
  fmt_text(ylab = c("Normalized expression", "", "")) %>%
  fmt_strip(label = paste(top_genes, "(differential expression)"),
            label_color = "black", label_fill = "grey60")

P_final <- P_row1 / P_row2 / P_row3

ggsave(
  P_final,
  filename = file.path(fig_dir, "CH4_Fig26_PHPT_Gene_Validation.pdf"),
  width = 11, height = 7.5
)


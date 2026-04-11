# ============================================================
# CH3_Fig22_Plot.R
# SHPT 候选基因验证: 轨迹相关 + 比例相关散点 + 差异表达小提琴
# 输入: ./out/gene_screening/gene_zscore_prop_SH.Rdata
#       ./out/trajectory/seu_with_pseudotime_V3.qs
# 输出: ./figures/CH3_Fig22_SHPT_Gene_Validation.pdf
# ============================================================

library(scMMR)
library(UtilsR)
library(Seurat)
library(qs)
library(tidyverse)
library(patchwork)

out_dir      <- "./out/gene_screening"
fig_dir      <- "./figures"
cluster2_col <- "Cluster 2"
top_genes    <- c("DUSP6", "GLRX", "HSPB1")

pal_lancet <- UtilsR::get_palette("lancet")

# --- 加载数据 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime_V3.qs")
seu_SH_Normal <- subset(seu, group %in% c("SHPT", "Normal"))
load(file.path(out_dir, "gene_zscore_prop_SH.Rdata"))
# → gene_z_bulk, gene_z_sc, props_bulk, prop_cl2_sc

# ============================================================
# Row 1: PlotDynamicFeatures — 轨迹相关
# ============================================================

p2 <- scMMR::PlotDynamicFeatures(
  subset(seu_SH_Normal, group == "SHPT"),
  pseudotime     = "Lineage2",
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
# Row 2: PlotScatter — 比例相关 (反卷积 z-score ~ Cluster 2 %)
# ============================================================

# --- 构建 Bulk scatter df (反卷积 Cluster 2 表达 z-score) ---
df_bulk <- data.frame(
  sample      = rownames(gene_z_bulk),
  cluster_pct = props_bulk[[cluster2_col]],
  gene_z_bulk[, top_genes, drop = FALSE],
  check.names = FALSE
)
colnames(df_bulk)[3:ncol(df_bulk)] <- paste0("gene_", top_genes)

# --- 构建 scRNA scatter df (样本级均值 z-score) ---
sc_samples <- intersect(names(prop_cl2_sc), colnames(gene_z_sc))
df_sc <- data.frame(
  sample      = sc_samples,
  cluster_pct = prop_cl2_sc[sc_samples],
  t(gene_z_sc[top_genes, sc_samples, drop = FALSE]),
  check.names = FALSE
)
colnames(df_sc)[3:ncol(df_sc)] <- paste0("gene_", top_genes)

n_bulk <- nrow(df_bulk)
n_sc   <- nrow(df_sc)

df_SH <- dplyr::bind_rows(
  df_bulk %>% dplyr::mutate(group = "Bulk"),
  df_sc   %>% dplyr::mutate(group = "SC")
) %>%
  dplyr::mutate(group = factor(
    group,
    levels = c("SC", "Bulk"),
    labels = c(
      paste0("scRNA (our, n = ", n_sc, ")"),
      paste0("Bulk (our, n = ", n_bulk, ")")
    )
  ))

p3 <- scMMR::PlotScatter(
  df_SH, "cluster_pct",
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
  seu_SH_Normal,
  features   = top_genes,
  group.by   = "group",
  split.by   = "celltype",
  palcolor   = c("#00468BFF", "#ED0000FF", "#42B540FF"),
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
  fmt_text(xlab = "Percentage of Cluster 2",
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

# --- 拼接 & 保存 ---
P_final <- P_row1 / P_row2 / P_row3

ggsave(
  P_final,
  filename = file.path(fig_dir, "CH3_Fig22_SHPT_Gene_Validation.pdf"),
  width = 10, height = 7.5
)


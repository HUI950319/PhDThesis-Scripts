# ============================================================
# CH3_Fig24_Plot.R
# SHPT 通路分析可视化: 小提琴 + 轨迹相关 + UMAP blend
# 输入: ./out/trajectory/seu_with_pseudotime_V3.qs
#       (AUCell scores 已由 CH3_Fig23_Data.R 写入 metadata)
# 输出: ./figures/CH3_Fig24_SHPT_Pathway_Analysis.pdf
# ============================================================

library(scMMR)
library(UtilsR)
library(Seurat)
library(qs)
library(patchwork)

fig_dir <- "./figures"

pal_lancet <- UtilsR::get_palette("lancet")

# --- 加载数据 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime_V3.qs")
seu_SH_Normal <- subset(seu, group %in% c("SHPT", "Normal"))

fea_pi3k <- "AUCell_REACTOME_PI3K_AKT_SIGNALING_IN_CANCER"
fea_mapk <- "AUCell_KEGG_MAPK_SIGNALING_PATHWAY"

# ============================================================
# Row 1: FeaturePlot3 — 通路活性小提琴 (分组 × 细胞类型)
# ============================================================

p4 <- scMMR::FeaturePlot3(
  seu_SH_Normal,
  features = c(fea_pi3k, fea_mapk),
  group.by = "group",
  split.by = "celltype",
  palcolor = c("#00468BFF", "#ED0000FF", "#42B540FF"),
  type     = "violin",
  add_box  = TRUE,
  ylab     = "AUCell score",
  stack    = FALSE
)

# ============================================================
# Row 2: PlotDynamicFeatures — 通路活性轨迹相关
# ============================================================

p5 <- scMMR::PlotDynamicFeatures(
  subset(seu, group == "SHPT"),
  pseudotime     = "Lineage2",
  features       = c(fea_pi3k, fea_mapk),
  group.by       = "celltype",
  smooth_k       = 2,
  pt.size        = 0.4,
  layer          = "data",
  exp_method     = "raw",
  stat_method    = "spearman",
  line_palcolor  = pal_lancet[2],
  point_palcolor = c("#00468BFF", "#42B540FF"),
  raster         = TRUE
)

# ============================================================
# Row 3: FeaturePlot2 — UMAP blend (DUSP6 + 通路活性)
# ============================================================

p6 <- scMMR::FeaturePlot2(
  subset(seu, group == "SHPT"),
  features          = c("DUSP6", fea_pi3k, fea_mapk),
  dims              = c(1, 2),
  color_blend_mode  = "blend",
  raster            = TRUE,
  raster_method     = "rasterise",
  legend.title      = ""
)

# ============================================================
# 组合三行面板
# ============================================================

P_row1 <- p4 %>%
  fmt_legend(collect = TRUE, legend_theme = leg1()) %>%
  fmt_text(ylab = c("AUCell score", ""),
           xlab = "Pathway activity across cell types and groups",
           title = "") %>%
  fmt_strip(label = c("PI3K AKT Pathway activity",
                       "MAPK Pathway activity"),
            label_color = "black")

P_row2 <- p5 %>%
  fmt_legend(collect = TRUE, legend_theme = leg1()) %>%
  fmt_text(ylab = c("AUCell score", ""), title = "") %>%
  fmt_strip(label = c("PI3K AKT Pathway activity",
                       "MAPK Pathway activity"),
            label_color = "black")

P_final <- P_row1 / P_row2 / p6

ggsave(
  P_final,
  filename = file.path(fig_dir, "CH3_Fig24_SHPT_Pathway_Analysis.pdf"),
  width = 10, height = 8
)


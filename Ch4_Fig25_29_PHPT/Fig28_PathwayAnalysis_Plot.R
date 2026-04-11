# ============================================================
# CH4_Fig28_Plot.R
# PHPT 通路分析可视化: 小提琴 + 轨迹相关 + UMAP blend
# 聚焦通路: TNFA_SIGNALING_VIA_NFKB
# 输入: ./out/trajectory/seu_with_pseudotime_V3.qs
#       (AUCell scores 已由 CH4_Fig27_Data.R 写入 metadata)
# 输出: ./figures/CH4_Fig28_PHPT_Pathway_Analysis.pdf
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
seu_PH_Normal <- subset(seu, group %in% c("PHPT", "Normal"))
seu_PH_Normal$celltype <- factor(seu_PH_Normal$celltype,
                                 levels = paste("cluster", c(0, 1, 3)))
seu_PH_Normal$group <- factor(seu_PH_Normal$group,
                              levels = c("PHPT", "Normal"))

fea_tnfa <- "AUCell_HALLMARK_TNFA_SIGNALING_VIA_NFKB"

# ============================================================
# Panel 1 (左上): FeaturePlot3 — 通路活性小提琴 (分组 × 细胞类型)
# ============================================================

p4 <- scMMR::FeaturePlot3(
  seu_PH_Normal,
  features = fea_tnfa,
  group.by = "group",
  split.by = "celltype",
  palcolor = c("#00468BFF", "#ED0000FF", "#0099B4FF"),
  type     = "violin",
  add_box  = TRUE,
  ylab     = "AUCell score",
  stack    = FALSE
)

# ============================================================
# Panel 2 (右上): PlotDynamicFeatures — 通路活性轨迹相关
# ============================================================

p5 <- scMMR::PlotDynamicFeatures(
  subset(seu, group == "PHPT" & Lineage3 < 11.5),
  pseudotime     = "Lineage3",
  features       = fea_tnfa,
  group.by       = "celltype",
  smooth_k       = 2,
  pt.size        = 0.4,
  layer          = "data",
  exp_method     = "raw",
  stat_method    = "spearman",
  line_palcolor  = pal_lancet[2],
  point_palcolor = c("#00468BFF", "#0099B4FF"),
  raster         = TRUE
)

# ============================================================
# Panel 3 (下): FeaturePlot2 — UMAP blend (ZFP36 + TNFA 通路)
# ============================================================

p6 <- scMMR::FeaturePlot2(
  subset(seu, group == "PHPT"),
  features         = c("ZFP36", fea_tnfa),
  dims             = c(1, 2),
  color_blend_mode = "blend",
  raster           = TRUE,
  raster_method    = "rasterise",
  legend.title     = ""
)

# ============================================================
# 组合: (小提琴 | 轨迹) / UMAP blend
# ============================================================

P_top_left <- p4 %>%
  fmt_legend(collect = TRUE, legend_theme = leg1()) %>%
  fmt_text(ylab = c("AUCell score", ""), xlab = "Group", title = "") %>%
  fmt_strip(label = "TNFA_SIGNALING_VIA_NFKB", label_color = "black")

P_top_right <- p5 %>%
  fmt_legend(collect = TRUE, legend_theme = leg1()) %>%
  fmt_text(ylab = c("AUCell score", ""), title = "") %>%
  fmt_strip(label = "TNFA_SIGNALING_VIA_NFKB", label_color = "black")

P_final <- (P_top_left | P_top_right) / p6

ggsave(
  P_final,
  filename = file.path(fig_dir, "CH4_Fig28_PHPT_Pathway_Analysis.pdf"),
  width = 10, height = 6
)


# ============================================================
# Fig16_18_Trajectory.R
# Figure 16: Monocle3 Pseudotime Lineages UMAP
# Figure 17: CytoTRACE2 / Monocle3 / Velocity pseudotime
# Figure 18: Dynamic gene expression (CASR, VDR, PTH)
# ============================================================

library(Seurat)
library(scMMR)
library(UtilsR)
library(qs)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

# --- 加载数据 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime.qs")
trajectory <- seu@tools$Monocle3$trajectory

# ============================================================
# Figure 16A: Lineage UMAP (all groups)
# ============================================================

scMMR::DimPlot2(  seu, group.by = "celltype", reduction = "umap",
  cells.highlight = TRUE,
  lineages          = c("Lineage1", "Lineage2", "Lineage3"),
  palcolor          = UtilsR::pal_lancet,
  lineages_palcolor = UtilsR::pal_lancet[2:6],
  lineages_arrow    = grid::arrow(length = grid::unit(0.2, "inches")),
  lineages_span     = 1,
  lineages_line_bg_stroke = 0.8,
  label = TRUE, label_insitu = TRUE, label.size = 4,
  title     = "Monocle3 Pseudotime Lineages",
  theme_use = "theme_blank"
)
ggplot2::ggsave("./figures/Fig16A_Lineage_UMAP.pdf", width = 6, height = 6)

# ============================================================
# Figure 16C: Lineage UMAP split by group
# ============================================================

scMMR::DimPlot2(
  seu, group.by = "celltype", reduction = "umap",
  split.by = "group",
  cells.highlight = TRUE,
  lineages          = c("Lineage1", "Lineage2", "Lineage3"),
  palcolor          = UtilsR::pal_lancet,
  lineages_palcolor = UtilsR::pal_lancet[2:6],
  lineages_arrow    = grid::arrow(length = grid::unit(0.2, "inches")),
  lineages_span     = 1,  lineages_line_bg_stroke = 0.8,
  label = TRUE, label_insitu = TRUE, label.size = 4,
  title     = "Monocle3 Pseudotime Lineages",
  theme_use = "theme_blank"
) %>% UtilsR::fmt_legend(merge_legends = TRUE)
ggplot2::ggsave("./figures/Fig16C_Lineage_UMAP_split.pdf", width = 14, height = 6)

# ============================================================
# Figure 17: Pseudotime comparison (CytoTRACE2 / Monocle3 / Velocity)
# ============================================================

spectral_pal     <- RColorBrewer::brewer.pal(11, "Spectral")
spectral_pal_rev <- rev(spectral_pal)

p1 <- scMMR::FeaturePlot2(
  seu, features = "CytoTRACE2_Score",
  palcolor  = spectral_pal,
  theme_use = UtilsR::theme_my() & ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
)

p2 <- scMMR::FeaturePlot2(
  seu, features = "Monocle3_Pseudotime",
  palcolor  = spectral_pal_rev,
  theme_use = UtilsR::theme_my() & ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
) + trajectory

p3 <- scMMR::FeaturePlot2(
  seu, features = "velocity_pseudotime_rev",  palcolor        = spectral_pal_rev,
  lower_quantile  = 0.1,
  upper_quantile  = 0.6,
  theme_use = UtilsR::theme_my() & ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
)

P1 <- (patchwork::wrap_plots(p1, p2, p3, ncol = 3) & UtilsR::theme_my()) %>%
  UtilsR::fmt_strip(
    label      = c("CytoTRACE2 Score", "Monocle3 Pseudotime", "RNA Velocity Pseudotime"),
    label_color = "black",
    label_fill  = c("grey60", "grey75", "grey90")
  ) %>%
  UtilsR::fmt_legend(
    collect = TRUE,
    title   = c("CytoTRACE2 Score", "Monocle3 Pseudotime", "RNA Velocity Pseudotime")
  )

P2 <- scMMR::FeaturePlot3(
  seu,
  features = c("CytoTRACE2_Score", "Monocle3_Pseudotime", "velocity_pseudotime_rev"),
  group.by = "celltype",
  type     = "box",
  stack    = FALSE
) & UtilsR::theme_my() & ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

P2 <- P2 %>%
  UtilsR::fmt_strip(
    label       = c("CytoTRACE2 Score", "Monocle3 Pseudotime", "RNA Velocity Pseudotime"),    label_color = "black",
    label_fill  = c("grey60", "grey75", "grey90")
  ) %>%
  UtilsR::fmt_axisText(x = 45) %>%
  UtilsR::fmt_text(ylab = c("CytoTRACE2 Score", "Monocle3 Pseudotime", "RNA Velocity Pseudotime")) %>%
  UtilsR::fmt_com(
    com_method   = list(c("cluster 0", "cluster 2"), c("cluster 0", "cluster 3")),
    label.y.prop = 0.85,
    tip.length   = 0.01
  ) %>%
  UtilsR::fmt_legend(collect = TRUE)

P1 / P2
ggplot2::ggsave("./figures/Fig17_Pseudotime.pdf", width = 14, height = 7)

# ============================================================
# Figure 18: Dynamic gene expression along lineages
# ============================================================

p_shpt <- scMMR::PlotDynamicFeatures(
  subset(seu, group == "SHPT"),
  pseudotime     = "Lineage2",
  features       = c("CASR", "VDR", "PTH"),
  group.by       = "celltype",
  smooth_k       = 2,
  pt.size        = 0.4,
  line_palcolor  = UtilsR::pal_lancet[2],
  point_palcolor = UtilsR::pal_lancet,
  ncol           = 3)

p_phpt <- scMMR::PlotDynamicFeatures(
  subset(seu, group == "PHPT"),
  pseudotime     = "Lineage3",
  features       = c("CASR", "VDR", "PTH"),
  group.by       = "celltype",
  smooth_k       = 2,
  pt.size        = 0.4,
  line_palcolor  = UtilsR::pal_lancet[2],
  point_palcolor = c("#00468BFF", "#0099B4FF"),
  ncol           = 3
)

(p_shpt / p_phpt & UtilsR::theme_my(panel.grid.minor = ggplot2::element_blank())) %>%
  UtilsR::fmt_text(xlab = "Monocele3 Pseudotime", ylab = "Expression (LogNormalized)") %>%
  UtilsR::fmt_legend(legend_theme = UtilsR::leg1())
ggplot2::ggsave("./figures/Fig18_DynamicFeatures.pdf", width = 12, height = 7)
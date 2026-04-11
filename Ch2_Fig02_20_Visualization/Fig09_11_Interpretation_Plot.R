# ============================================================
# Fig09_11_Interpretation.R
# Figure 9:  Gene Importance    (A: CircleLollipop, B: Jaccard)
# Figure 10:  KEGG Importance    (A: CircleLollipop, B: Jaccard)
# Figure 11: Regulon Importance (A: CircleLollipop, B: Jaccard)
# ============================================================

library(qs)
library(ggplot2)
library(dplyr)
library(patchwork)
library(UtilsR)

# --- 加载中间数据 ---
out_dir <- "./out/09_model_interpretation"
dat <- qs::qread(file.path(out_dir, "fig9_intermediate.qs"))

gene_pc      <- dat$gene_pc
hallmark_top <- dat$hallmark_top
kegg_top     <- dat$kegg_top
regulon_top  <- dat$regulon_top
top_k_class  <- dat$params$top_k_class
# ============================================================
# Jaccard Heatmaps (Figure 9B / 10B / 11B)
# ============================================================

p1 <- UtilsR::PlotHeatmapJaccard(
  data         = gene_pc,
  group_col    = "cell_type",
  name_col     = "gene",
  group_levels = names(UtilsR::pal_paraSC),
  title        = sprintf("Gene Importance Jaccard (top %d IG genes)", top_k_class),
  show_strip   = TRUE,
  strip_colors = UtilsR::pal_paraSC,
  strip_size   = 0.25,
  width = 8, height = 8
)

p2 <- UtilsR::PlotHeatmapJaccard(
  data         = kegg_top,
  group_col    = "cell_type",
  name_col     = "name",
  group_levels = names(UtilsR::pal_paraSC),
  title        = sprintf("KEGG Pathway Jaccard (top %d)", top_k_class),
  show_strip   = TRUE,
  strip_colors = UtilsR::pal_paraSC,
  strip_size   = 0.25,
  width = 8, height = 8
)
p3 <- UtilsR::PlotHeatmapJaccard(
  data         = regulon_top,
  group_col    = "cell_type",
  name_col     = "name",
  group_levels = names(UtilsR::pal_paraSC),
  title        = sprintf("Regulon Jaccard (top %d SCENIC regulons)", top_k_class),
  show_strip   = TRUE,
  strip_colors = UtilsR::pal_paraSC,
  strip_size   = 0.25,
  width = 8, height = 8
)

patchwork::wrap_plots(p1, p2, p3, nrow = 1) %>% UtilsR::fmt_legend(merge_legends = TRUE)
ggplot2::ggsave("./figures/Fig9_Jaccard_Heatmaps.pdf", width = 21, height = 7.5)

# ============================================================
# CircleLollipop Polar Plots (Figure 9A / 10A / 11A)
# ============================================================

# --- Figure 9A: Gene importance ---
gene_pc %>%
  UtilsR::PlotCircleLollipop(
    group_col    = "cell_type",
    name_col     = "gene",
    value_col    = "importance",
    group_levels = names(UtilsR::pal_paraSC),    top_n        = 5,
    value_scale  = "group",
    label_side   = "outer",
    yaxis_label  = "Gene Importance (IG)",
    label_cex    = 0.65,
    sector_cex   = 0.65,
    track_height = 0.3,
    gap_after    = 15,
    filename     = "./figures/Fig9A_Gene_Polar.pdf",
    width = 8, height = 8
  )

# --- Figure 10A: KEGG importance ---
kegg_top %>%
  UtilsR::PlotCircleLollipop(
    group_col    = "cell_type",
    name_col     = "name_clean",
    value_col    = "importance",
    group_levels = names(UtilsR::pal_paraSC),
    top_n        = 5,
    value_scale  = "group",
    label_side   = "outer",
    yaxis_label  = "KEGG Importance (IG)",
    label_cex    = 0.65,
    sector_cex   = 0.65,
    track_height = 0.3,
    gap_after    = 15,
    filename     = "./figures/Fig10A_KEGG_Polar.pdf",
    width = 8, height = 8
  )
# --- Figure 11A: Regulon importance ---
regulon_top %>%
  UtilsR::PlotCircleLollipop(
    group_col    = "cell_type",
    name_col     = "name_clean",
    value_col    = "importance",
    group_levels = names(UtilsR::pal_paraSC),
    top_n        = 5,
    value_scale  = "group",
    label_side   = "outer",
    yaxis_label  = "Regulon Importance (IG)",
    label_cex    = 0.65,
    sector_cex   = 0.65,
    track_height = 0.3,
    gap_after    = 15,
    filename     = "./figures/Fig11A_Regulon_Polar.pdf",
    width = 8, height = 8
  )
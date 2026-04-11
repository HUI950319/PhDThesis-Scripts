# ============================================================
# CH3_Fig21_Plot.R
# Figure 21: SHPT 轨迹相关火山图 + 三集韦恩图 (inset)
# A: PlotCorrelation 火山图 (Spearman rho vs -log10 padj)
#    点大小 = |DE logFC|, 透明度 = prop_cor
# B: Venn inset (TraceGene / DE / PropCor 三集交集)
# ============================================================

library(tidyverse)
library(patchwork)
library(scMMR)
library(UtilsR)
library(LabR)

pal_lancet <- as.character(UtilsR::pal_lancet)

# --- 加载中间数据 ---
out_dir <- "./out/gene_screening"
res_merged_SH <- read.csv(file.path(out_dir, "res_merged_SH.csv"))
load(file.path(out_dir, "venn_gene_lists_SH.Rdata"))

# ============================================================
# Panel A: 轨迹相关火山图
# ============================================================

# 正相关 + 上调基因子集
plot_df <- res_merged_SH %>%
  dplyr::filter(prop_cor > 0 & de_logFC > 0)

p1 <- UtilsR::PlotCorrelation(
  plot_df,
  p.col      = "padjust",
  p.cutoff   = 0.00001,
  size.by    = "de_abs_logFC",
  alpha.by   = "prop_cor",
  title      = ""
) %>%
  LabR::fmt_legend(legend.position = "right")

# ============================================================
# Panel B (inset): 三集韦恩图
# ============================================================

p2 <- UtilsR::plt_upset(
  venn_gene_lists_SH,
  output     = "venn",
  label_size    = 2.5,
  set.name.size = 2.5,
  label         = "count",
  label_color   = "white",
  colors        = pal_lancet[3]
) & Seurat::NoLegend()

# ============================================================
# 组合: 火山图 + 韦恩图 inset
# ============================================================

p_fig21 <- (
  p1 %>%
    LabR::fmt_strip(label = "Trajectory-correlated Volcano - SHPT",
                    label_color = "black") %>%
    LabR::fmt_text(ylab = "Spearman correlation", title = "") &
    LabR::theme_my() %>%
    LabR::fmt_legend(legend_theme = LabR::leg1())
) +
  patchwork::inset_element(p2, left = 0.53, right = 1, bottom = 0.3, top = 0.7)

ggplot2::ggsave("./figures/CH3_Fig21_SHPT_Trajectory_Volcano_Venn.pdf",
                p_fig21, width = 9, height = 6)


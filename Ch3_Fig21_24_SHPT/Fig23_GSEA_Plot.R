# ============================================================
# CH3_Fig23_Plot.R
# SHPT 通路筛选可视化: GSEA Volcano + DUSP6~通路活性 散点内嵌
# 输入: ./out/pathway_screening/enrich_plot_SH.csv
#       ./out/pathway_screening/scatter_df_SH.csv
# 输出: ./figures/CH3_Fig23_SHPT_Pathway_Screening.pdf
# ============================================================

library(scMMR)
library(UtilsR)
library(qs)
library(tidyverse)
library(patchwork)

out_dir <- "./out/pathway_screening"
fig_dir <- "./figures"

pal_lancet <- UtilsR::get_palette("lancet")

# --- 加载数据 ---
res     <- qs::qread(file.path(out_dir, "gsea_res_SH.qs"))
df_long <- read.csv(file.path(out_dir, "scatter_df_SH.csv"),
                    check.names = FALSE)

# ============================================================
# Panel A: GSEA Volcano (PlotCorrelation)
# ============================================================

p1 <- UtilsR::PlotCorrelation(
  res$enrichment,
  name.col  = "Description",
  score.col = "NES",
  p.col     = "p.adjust",
  p.cutoff  = 0.05,
  size.by   = "pvalue",
  ylab      = "Normalized Enrichment Score"
) %>%
  fmt_legend(legend.position = "right")

# ============================================================
# Panel B (inset): DUSP6 ~ 通路活性 散点
# ============================================================

p2 <- scMMR::PlotScatter(
  df_long, "DUSP6", "value",
  group.by    = "Pathways",
  split.by    = NULL,
  show.rug    = TRUE,
  point.size  = 0.8,
  point.alpha = 0.6,
  palcolor    = pal_lancet[1:2],
  point.color = "#984ea3",
  rug.color   = "#7fc97f",
  raster      = TRUE
)

p2_fmt <- p2 %>%
  fmt_legend(legend.position = c(0.95, 0.1),
             legend_theme = leg1(), scale = 0.4) %>%
  fmt_text(xlab  = "Expression of DUSP6",
           ylab  = "Pathway activity (AUCell score)",
           title = "")

# ============================================================
# 组合: Volcano + 散点 inset
# ============================================================

P_final <- (p1 %>%
  fmt_strip(label = "GSEA Volcano - SHPT", label_color = "black") %>%
  fmt_text(ylab = "Normalized Enrichment Score", title = "") &
  theme_my() %>%
  fmt_legend(legend_theme = leg1())) +
  patchwork::inset_element(p2_fmt,
                           left = 0.45, right = 0.95,
                           bottom = 0.05, top = 0.75)

ggsave(
  P_final,
  filename = file.path(fig_dir, "CH3_Fig23_SHPT_Pathway_Screening.pdf"),
  width = 10, height = 6
)


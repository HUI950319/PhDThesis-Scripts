# ============================================================
# CH4_Fig25_Plot.R
# PHPT 候选基因三维度筛选: Volcano + Venn
# 输入: ./out/gene_screening/res_merged_PH.csv
#       ./out/gene_screening/venn_gene_lists_PH.Rdata
# 输出: ./figures/CH4_Fig25_PHPT_Trajectory_Volcano_Venn.pdf
# ============================================================

library(scMMR)
library(UtilsR)
library(tidyverse)
library(patchwork)

out_dir <- "./out/gene_screening"
fig_dir <- "./figures"

pal_lancet <- UtilsR::get_palette("lancet")

# --- 加载数据 ---
res_merged_PH <- read.csv(file.path(out_dir, "res_merged_PH.csv"),
                          check.names = FALSE)
load(file.path(out_dir, "venn_gene_lists_PH.Rdata"))

# ============================================================
# Panel A: Volcano (正向 prop_cor + 正向 de_logFC 子集)
# ============================================================

p1 <- UtilsR::PlotCorrelation(
  res_merged_PH %>% dplyr::filter(prop_cor > 0 & de_logFC > 0),
  p.col     = "prop_pval",
  p.cutoff  = 0.00001,
  size.by   = "de_logFC",
  alpha.by  = "prop_cor",
  col.pos   = pal_lancet[2],
  col.neg   = pal_lancet[1],
  title     = "Correlation Volcano - PHPT"
) %>%
  fmt_legend(legend.position = "right")

# ============================================================
# Panel B (inset): Venn diagram
# ============================================================

p2 <- UtilsR::plt_upset(
  venn_gene_lists_PH,
  output         = "venn",
  label_size     = 2.5,
  set.name.size  = 2.5,
  label          = "count",
  label_color    = "white",
  colors         = pal_lancet[3]
) & NoLegend()

# ============================================================
# 组合: Volcano + Venn inset
# ============================================================

P_final <- (p1 %>%
  fmt_strip(label = "Trajectory-correlated Volcano - PHPT",
            label_color = "black") %>%
  fmt_text(ylab = "Spearman correlation", title = "") &
  theme_my() %>%
  fmt_legend(legend_theme = leg1())) +
  patchwork::inset_element(p2, left = 0.53, right = 1,
                           bottom = 0.3, top = 0.7)

ggsave(
  P_final,
  filename = file.path(fig_dir, "CH4_Fig25_PHPT_Trajectory_Volcano_Venn.pdf"),
  width = 9, height = 6
)


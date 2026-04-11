# ==============================================================================
# Figure 20: scMMR vs BayesPrism 解卷积一致性验证
# ==============================================================================
#
# 描述：
#   散点图展示 scMMR DNN 与 BayesPrism 两种方法估计的细胞亚群比例的一致性。
#   上图: PHPT Cluster 3（自有 Bulk + PRJNA516535）
#   下图: SHPT Cluster 2（自有 Bulk）
#
# 输入：
#   - ./out/Figures/Figure20/plot_data2.Rdata (df_PH_bp, df_SH_bp)
#
# 输出：
#   - ./out/Figures/Figure20/bayesprism_validation.pdf (7 x 8 inch)
#
# 依赖：
#   - scMMR::PlotScatter()
#   - LabR: theme_my(), fmt_text(), fmt_legend(), fmt_panel(), fmt_strip(), leg1()
#   - UtilsR: pal_lancet
#   - patchwork
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scMMR)
library(LabR)

# ==============================================================================
# 1. 加载数据
# ==============================================================================

load("./out/Figures/Figure20/plot_data2.Rdata")

# ==============================================================================
# 2. PHPT: scMMR vs BayesPrism (Cluster 3)
# ==============================================================================

p_bp1 <- scMMR::PlotScatter(
  df_PH_bp, "cluster_pct", "bp_pct",
  group.by = "group",
  point.size = 3, point.alpha = 1,
  smooth.size = 1.5,
  palette = "lancet",
  global.cor = TRUE
)

# ==============================================================================
# 3. SHPT: scMMR vs BayesPrism (Cluster 2)
# ==============================================================================

p_bp2 <- scMMR::PlotScatter(
  df_SH_bp, "cluster_pct", "bp_pct",
  group.by = "group",
  point.size = 3, point.alpha = 1,
  smooth.size = 1.5,
  palette = "lancet",
  smooth.color = "#925E9F",
  global.cor = TRUE
)

# ==============================================================================
# 4. 组合排版
# ==============================================================================

P_bp1 <-
  (p_bp1 & theme_my()) %>%
  fmt_text(
    xlab = "Percentage of Cluster 3 (scMMR)",
    ylab = "Percentage of Cluster 3 (BayesPrism)",
    title = "",
    legend_title = "Datasets"
  ) %>%
  fmt_legend(
    legend_theme = leg1(),
    legend.position = c(0.99, 0.02),
    scale = 0.7
  ) %>%
  fmt_panel(grid = "major") %>%
  fmt_strip(
    label = "scMMR vs BayesPrism (PHPT Cluster 3)",
    label_color = "black",
    label_fill = "grey60"
  )

P_bp2 <-
  (p_bp2 & theme_my()) %>%
  fmt_text(
    xlab = "Percentage of Cluster 2 (scMMR)",
    ylab = "Percentage of Cluster 2 (BayesPrism)",
    title = "",
    legend_title = "Datasets"
  ) %>%
  fmt_legend(
    legend_theme = leg1(),
    legend.position = c(0.99, 0.02),
    scale = 0.8
  ) %>%
  fmt_panel(grid = "major") %>%
  fmt_strip(
    label = "scMMR vs BayesPrism (SHPT Cluster 2)",
    label_color = "black",
    label_fill = "grey90"
  )

# ==============================================================================
# 5. 保存
# ==============================================================================

ggplot2::ggsave(
  "./out/Figures/Figure20/bayesprism_validation.pdf",
  plot = P_bp1 / P_bp2,
  width = 7, height = 8
)
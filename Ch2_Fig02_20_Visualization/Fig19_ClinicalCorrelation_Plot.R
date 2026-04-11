# ==============================================================================
# Figure 19: 临床相关性分析 — 细胞亚群比例与临床指标的相关性
# ==============================================================================
#
# 描述：
#   绘制甲状旁腺细胞亚群比例与临床指标（肿瘤大小、PTH水平）的散点相关图。
#   Bulk RNA-seq 的细胞比例通过 scMMR 包中的解卷积函数计算得到。
#
# 数据来源：
#   - scRNA-seq (our, n=6)
#   - Bulk RNA-seq (our, n=12, 解卷积估计)
#   - scRNA-seq (GSE190773, n=5)
#   - Bulk RNA-seq (PRJNA516535, n=10, 解卷积估计)
#
# 输入：
#   - ./out/Figures/Figure19/plot_data.Rdata
#     包含 df1-df6 数据框（各数据集的亚群比例 + 临床信息）
#
# 输出：
#   - ./out/Figures/Figure19/percent_clinical.pdf (11 x 8 inch)
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

# --- 调色板 ---
pal_lancet <- UtilsR::pal_lancet

# ==============================================================================
# 1. 数据准备
# ==============================================================================

load("./out/Figures/Figure19/plot_data.Rdata")

# PHPT: 4个数据集 (自有scRNA + 自有Bulk + GSE190773 + PRJNA516535)
df_PH <- 
  dplyr::bind_rows(
    df1 %>% dplyr::mutate(group = "PRJNA516535"),
    df2 %>% dplyr::mutate(group = "GSE190773"),
    df3 %>% dplyr::mutate(group = "Bulk"),
    df4 %>% dplyr::mutate(group = "SC")
  ) %>% 
  dplyr::mutate(group = factor(
    group, 
    levels = c("SC", "Bulk", "GSE190773", "PRJNA516535"),
    labels = c("scRNA (our, n = 6)", "Bulk (our, n = 12)", 
               "scRNA (GSE190773, n = 5)", "Bulk (PRJNA516535, n = 10)")
  ))

# SHPT: 2个数据集 (自有scRNA + 自有Bulk)
df_SH <- 
  dplyr::bind_rows(
    df5 %>% dplyr::mutate(group = "Bulk"),
    df6 %>% dplyr::mutate(group = "SC")
  ) %>% 
  dplyr::mutate(group = factor(
    group, 
    levels = c("SC", "Bulk"),
    labels = c("scRNA (our, n = 6)", "Bulk (our, n = 12)")
  ))

# ==============================================================================
# 2. PHPT 散点图: Cluster 3 比例 vs 肿瘤大小 / PTH
# ==============================================================================

p1 <- scMMR::PlotScatter(
  df_PH, "cluster_pct", "tumor_size", 
  group.by = "group", 
  point.size = 3, point.alpha = 1,
  smooth.size = 1.5,
  palette = "lancet",
  global.cor = TRUE
)

p2 <- scMMR::PlotScatter(
  df_PH, "cluster_pct", "PTH",
  group.by = "group", 
  point.size = 3, point.alpha = 1,
  smooth.size = 1.5,
  palcolor = pal_lancet,
  global.cor = TRUE
)

# ==============================================================================
# 3. SHPT 散点图: Cluster 2 比例 vs 肿瘤大小 / PTH
# ==============================================================================

p3 <- scMMR::PlotScatter(
  df_SH, "cluster_pct", "tumor_size", 
  group.by = "group", 
  point.size = 3, point.alpha = 1,
  smooth.size = 1.5,
  palette = "lancet", 
  smooth.color = "#925E9F",
  global.cor = TRUE
) & ggplot2::scale_y_continuous(limits = c(1, 3), breaks = c(1, 2, 3))

p4 <- scMMR::PlotScatter(
  df_SH, "cluster_pct", "PTH",
  group.by = "group", 
  point.size = 3, point.alpha = 1,
  smooth.size = 1.5,
  smooth.color = "#925E9F",
  palette = "lancet",
  global.cor = TRUE
)

# ==============================================================================
# 4. 组合排版
# ==============================================================================

# 上半部分: PHPT — Cluster 3
P1 <- 
  ((p1 | p2) & theme_my()) %>% 
  fmt_text(
    xlab = "Percentage of Cluster 3", 
    ylab = c("Tumor size (cm)", "PTH levels(pg/mL)"), 
    title = "", 
    legend_title = "Datasets"
  ) %>% 
  fmt_legend(
    legend_theme = leg1(), 
    legend.position = c(0.99, 0.02), 
    scale = 0.5
  ) %>%
  fmt_panel(grid = "major") %>%
  fmt_strip(
    label = c("Para_cluster 3 vs tumor size in PHPT", 
              "Para_cluster 3 vs PTH in PHPT"), 
    label_color = "black",
    label_fill = "grey60"
  )

# 下半部分: SHPT — Cluster 2
P2 <- 
  ((p3 | p4) & theme_my()) %>% 
  fmt_text(
    xlab = "Percentage of Cluster 2", 
    ylab = c("Tumor size (cm)", "PTH levels (pg/mL)"), 
    title = "", 
    legend_title = "Datasets"
  ) %>% 
  fmt_legend(
    legend_theme = leg1(), 
    legend.position = c(0.99, 0.02), 
    scale = 0.5
  ) %>%
  fmt_panel(grid = "major") %>%
  fmt_strip(
    label = c("Para_cluster 2 vs tumor size in SHPT", 
              "Para_cluster 2 vs PTH in SHPT"), 
    label_color = "black",
    label_fill = "grey90"
  )

# ==============================================================================
# 5. 保存
# ==============================================================================

ggplot2::ggsave(
  "./out/Figures/Figure19/percent_clinical.pdf",
  plot = P1 / P2,
  width = 11, height = 8
)
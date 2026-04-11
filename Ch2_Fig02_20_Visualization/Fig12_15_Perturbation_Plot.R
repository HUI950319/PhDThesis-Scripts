# ============================================================
# Fig12_15_Perturbation.R
# Figure 12: Alluvial plot (cell type proportions by group)
# Figure 13: Proportion perturbation (scMMR vs miloR)
# Figure 14: Expression perturbation (scMMR vs Augur)
# Figure 15: Cross-comparison (PHPT vs SHPT)
# ============================================================

library(qs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(UtilsR)
library(scMMR)

# --- 加载扰动分析结果 ---
base_dir  <- "./out/perturbation"

# Proportion: scMMR RankPercent + miloR
pct_phpt <- read.csv(file.path(base_dir, "scMMR_percent",
                               "percent_PHPT_vs_Normal_summary.csv"))
pct_shpt <- read.csv(file.path(base_dir, "scMMR_percent",
                               "percent_SHPT_vs_Normal_summary.csv"))
pct_nhood_phpt <- read.csv(file.path(base_dir, "scMMR_percent",
                                     "percent_PHPT_vs_Normal_da_results.csv"))
pct_nhood_shpt <- read.csv(file.path(base_dir, "scMMR_percent",
                                     "percent_SHPT_vs_Normal_da_results.csv"))

milo_phpt <- read.csv(file.path(base_dir, "miloR",
                                "miloR_PHPT_vs_Normal_summary.csv"))
milo_shpt <- read.csv(file.path(base_dir, "miloR",
                                "miloR_SHPT_vs_Normal_summary.csv"))
milo_nhood_phpt <- read.csv(file.path(base_dir, "miloR",
                                      "miloR_PHPT_vs_Normal_da_results.csv"))
milo_nhood_shpt <- read.csv(file.path(base_dir, "miloR",
                                      "miloR_SHPT_vs_Normal_da_results.csv"))

# Expression: scMMR RankPerturbation + Augur
pert_phpt <- read.csv(file.path(base_dir, "scMMR_perturbation",
                                "perturbation_PHPT_vs_Normal_results.csv"))
pert_shpt <- read.csv(file.path(base_dir, "scMMR_perturbation",
                                "perturbation_SHPT_vs_Normal_results.csv"))
augur_phpt <- read.csv(file.path(base_dir, "augur",
                                 "augur_PHPT_vs_Normal_auc.csv"))
augur_shpt <- read.csv(file.path(base_dir, "augur",
                                 "augur_SHPT_vs_Normal_auc.csv"))

# --- 构建汇总数据 ---
seu_meta <- qs::qread(file.path(base_dir, "seu_meta.data.qs"))
ct_pct <- seu_meta %>%
  dplyr::count(celltype) %>%
  dplyr::mutate(pct = n / sum(n) * 100) %>%
  dplyr::select(celltype, pct)

ct_diff <- seu_meta %>%
  dplyr::count(group, celltype) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(pct_group = n / sum(n) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::select(group, celltype, pct_group) %>%
  tidyr::pivot_wider(names_from = group, values_from = pct_group, values_fill = 0) %>%
  dplyr::transmute(
    celltype,
    diff_PHPT = PHPT - Normal,
    diff_SHPT = SHPT - Normal
  ) %>%
  tidyr::pivot_longer(cols = c(diff_PHPT, diff_SHPT),
               names_to = "comparison", values_to = "diff_pct") %>%
  dplyr::mutate(
    comparison = dplyr::recode(comparison,
                        diff_PHPT = "PHPT_vs_Normal",
                        diff_SHPT = "SHPT_vs_Normal"),
    abs_diff_pct = abs(diff_pct)
  )

prop_data <- dplyr::bind_rows(
  dplyr::bind_rows(pct_phpt, pct_shpt) %>%
    dplyr::transmute(celltype = cell_type, comparison,              value = abs(mean_logFC), method = "scMMR"),
  dplyr::bind_rows(milo_phpt, milo_shpt) %>%
    dplyr::transmute(celltype, comparison,
              value = abs(mean_logFC), method = "miloR")
) %>%
  dplyr::left_join(ct_pct, by = "celltype") %>%
  dplyr::left_join(ct_diff, by = c("celltype", "comparison"))

nhood_long <- dplyr::bind_rows(
  pct_nhood_phpt %>%
    dplyr::transmute(celltype = cell_type_majority, logFC,
              comparison = "PHPT_vs_Normal", method = "scMMR RankPercent"),
  pct_nhood_shpt %>%
    dplyr::transmute(celltype = cell_type_majority, logFC,
              comparison = "SHPT_vs_Normal", method = "scMMR RankPercent"),
  milo_nhood_phpt %>%
    dplyr::transmute(celltype, logFC,
              comparison = "PHPT_vs_Normal", method = "miloR"),
  milo_nhood_shpt %>%
    dplyr::transmute(celltype, logFC,
              comparison = "SHPT_vs_Normal", method = "miloR")
) %>%
  dplyr::filter(celltype != "Mixed")

expr_data <- dplyr::bind_rows(
  dplyr::bind_rows(pert_phpt, pert_shpt) %>%
    dplyr::transmute(celltype = cell_type, comparison,
              value = score, method = "scMMR"),
  dplyr::bind_rows(augur_phpt, augur_shpt) %>%    dplyr::transmute(celltype = cell_type, comparison,
              value = auc, method = "Augur")
) %>%
  dplyr::left_join(ct_pct, by = "celltype") %>%
  dplyr::left_join(ct_diff, by = c("celltype", "comparison"))

# ============================================================
# Figure 12: Alluvial plot
# ============================================================

p1 <- scMMR::PlotAlluvia2(
  seu_meta,
  by         = "group",
  fill       = "celltype",
  palcolor   = UtilsR::pal_paraSC,
  flow.alpha = 0.25,
  bar.width  = 0.52,
  show.label = TRUE, label.size = 2.6,
  show.pct   = TRUE,
  curve.type = "linear"
)

(p1 & UtilsR::theme_alluvia()) %>%
  UtilsR::fmt_legend(legend_theme = UtilsR::leg1(), scale = 1.3)
ggplot2::ggsave("./figures/Fig12_Alluvial.pdf", width = 15, height = 9)

# ============================================================
# Figure 13: Proportion perturbation (scMMR vs miloR)
# ============================================================
p2 <- UtilsR::PlotButterfly2(
  data        = nhood_long,
  stat.by     = "celltype",
  type        = "violin_box",
  value.by    = "logFC",
  group.by    = "comparison",
  fill.by     = "method",
  levels      = "up",
  palette     = c("#B2182B", "#2166AC"),
  ref.line    = 0,
  color_by    = "fixed",
  gradient_colors = c("#2166AC", "white", "#B2182B"),
  base_size       = 14,
  violin.scale    = "width",
  legend.position = c(0.1, 0.2),
  legend_theme    = UtilsR::theme_legend1()
)

p3 <- UtilsR::PlotRankCor(
  data      = prop_data,
  stat.by   = "celltype",
  value.by  = "value",
  group.by  = "comparison",
  method.by = "method",
  size_by     = "pct",
  size_range  = c(2, 6),
  alpha_by    = "diff_pct",
  alpha_range = c(0.25, 1),
  palette     = UtilsR::pal_lancet[c(3:4)],  title           = "Proportion Perturbation: scMMR RankPercent vs miloR",
  label.size      = 4,
  base_size       = 14,
  legend.position = "right",
  use_rank        = TRUE
)

p2 %>% UtilsR::fmt_legend(legend.position = "right", collect = TRUE, legend_theme = UtilsR::leg1()) /
  p3 %>% UtilsR::fmt_legend(legend.position = "right", collect = TRUE, legend_theme = UtilsR::leg1())
ggplot2::ggsave("./figures/Fig13_Proportion.pdf", width = 13, height = 10)

# ============================================================
# Figure 14: Expression perturbation (scMMR vs Augur)
# ============================================================

p4 <- expr_data %>%
  dplyr::mutate(value = dplyr::if_else(
    celltype %in% c("Mast cells", "Cycling T cells", "Monocytes"),
    0.5 * value, value)) %>%
  UtilsR::PlotButterfly(
    type     = "bar_dodge",
    stat.by  = "celltype",
    fill.by  = "method",
    value.by = "value",
    group.by = "comparison",
    palette  = c("#B2182B", "#2166AC"),
    levels   = "up",
    xlab     = "|log2FC| (Proportion Perturbation)",
    base_size       = 12,    legend.position = c(0.1, 0.2),
    legend_theme    = UtilsR::theme_legend1(),
    minmax          = TRUE
  )

p5 <- expr_data %>%
  dplyr::mutate(value = dplyr::if_else(
    celltype %in% c("Mast cells", "Cycling T cells", "Monocytes"),
    0.5 * value, value)) %>%
  dplyr::filter(celltype %in% intersect(pert_phpt$cell_type, augur_phpt$cell_type)) %>%
  UtilsR::PlotRankCor(
    stat.by   = "celltype",
    value.by  = "value",
    group.by  = "comparison",
    method.by = "method",
    size_by     = "pct",
    size_range  = c(2, 6),
    alpha_by    = "diff_pct",
    alpha_range = c(0.2, 1),
    palette     = UtilsR::pal_lancet[c(3:4)],
    title           = "Proportion Perturbation: scMMR RankPercent vs miloR",
    label.size      = 3,
    legend.position = "right",
    use_rank        = TRUE
  )

p4 %>% UtilsR::fmt_legend(legend.position = "right", collect = TRUE, legend_theme = UtilsR::leg1()) /
  p5 %>% UtilsR::fmt_legend(legend.position = "right", collect = TRUE, legend_theme = UtilsR::leg1())
ggplot2::ggsave("./figures/Fig14_Expression.pdf", width = 13, height = 10)
# ============================================================
# Figure 15: Cross-comparison (PHPT vs SHPT)
# ============================================================

# --- 14A: Proportion perturbation PHPT vs SHPT ---
p7 <- nhood_long %>%
  dplyr::group_by(celltype, comparison, method) %>%
  dplyr::summarise(value = mean(logFC), .groups = "drop") %>%
  dplyr::filter(method == "scMMR RankPercent") %>%
  as.data.frame() %>%
  dplyr::left_join(ct_pct, by = "celltype") %>%
  dplyr::mutate(value = dplyr::if_else(celltype == "iTAFs", 0.5 * value, value)) %>%
  UtilsR::PlotRankCor(
    stat.by   = "celltype",
    value.by  = "value",
    method.by = "comparison",
    size_by     = "pct",
    size_range  = c(2, 6),
    alpha_by    = "pct",
    alpha_range = c(0.3, 1),
    palette     = "#925E9F",
    title       = "Proportion Perturbation: PHPT vs SHPT",
    xlab        = "PHPT vs Normal (log\u2082FC)",
    ylab        = "SHPT vs Normal (log\u2082FC)",
    label.size      = 4,
    base_size       = 14,
    legend.position = "right",
    use_rank   = FALSE,
    show.diag  = FALSE  ) %>%
  UtilsR::fmt_ref(x = 0, y = 0, color = "#B2182B")

# --- 14B: Expression perturbation PHPT vs SHPT ---
tmp <- dplyr::bind_rows(
  prop_data %>% dplyr::filter(method == "scMMR") %>% dplyr::mutate(type = "prop"),
  expr_data %>% dplyr::filter(method == "scMMR") %>% dplyr::mutate(type = "expr")
)

p8 <- UtilsR::PlotRankCor(
  data      = tmp %>% dplyr::filter(type == "expr"),
  stat.by   = "celltype",
  value.by  = "value",
  method.by = "comparison",
  size_by     = "pct",
  size_range  = c(2, 6),
  alpha_by    = "pct",
  alpha_range = c(0.25, 1),
  palette     = "#79AF97",
  title       = "Expression Perturbation: PHPT vs SHPT",
  xlab        = "PHPT vs Normal (normalized score)",
  ylab        = "SHPT vs Normal (normalized score)",
  label.size      = 4,
  base_size       = 14,
  legend.position = "right",
  use_rank  = FALSE,
  scale_x   = "minmax",
  scale_y   = "minmax"
)
(p7 %>% UtilsR::fmt_legend(legend.position = "right", collect = TRUE, legend_theme = UtilsR::leg1()) /
  p8 %>% UtilsR::fmt_legend(legend.position = "right", collect = TRUE, legend_theme = UtilsR::leg1())) %>%
  UtilsR::fmt_text(title = "") %>%
  UtilsR::fmt_strip(
    label      = c("Proportion Perturbation: PHPT vs SHPT",
                    "Expression Perturbation: PHPT vs SHPT"),
    label_color = "black",
    label_fill  = c("#925E9F66", "#79AF9766")
  ) %>%
  UtilsR::flatten_patchwork(ncol = 1)
ggplot2::ggsave("./figures/Fig15_CrossComparison.pdf", width = 9, height = 11)
# ============================================================
# Fig05_EvalEmbedding.R
# Embedding Quality Evaluation
# Panel A: Dumbbell Chart + inset paired boxplot (Silhouette)
# Panel B: Standardized distance scatter (PlotScatter2)
# Panel C: PCA vs DNN UMAP comparison
# ============================================================

library(scMMR)
library(UtilsR)
library(ggplot2)
library(patchwork)
library(qs)

# --- 路径 ---
eval_dir <- "./out/evaluate_embedding"

# --- 加载数据 ---
sil_df <- read.csv(file.path(eval_dir, "silhouette_per_type.csv"),
                   stringsAsFactors = FALSE)
sil_df$delta     <- sil_df$sil_emb - sil_df$sil_ref
sil_mean_emb     <- mean(sil_df$sil_emb)
sil_mean_ref     <- mean(sil_df$sil_ref)

all_res   <- qs::qread(file.path(eval_dir, "eval_results.qs"))
eval_res  <- all_res$eval_res
dist_data <- eval_res$consistency$dist_cor
d_emb_all <- dist_data$d_emb
d_ref_all <- dist_data$d_ref
rho_val   <- dist_data$correlation

umap_data   <- qs::qread(file.path(eval_dir, "umap_comparison_data.qs"))

# ============================================================
# Panel A: Dumbbell Chart + inset paired boxplot
# ============================================================
p1 <- UtilsR::PlotDumbbell(sil_df,
                   x_col = "sil_ref",
                   xend_col = "sil_emb",
                   y_col = "cell_type",
                   label_x = "PCA",
                   label_xend = "scMMR",
                   legend_theme = UtilsR::theme_legend1(),
                   show_mean_lines = TRUE,
                   show_delta = TRUE,
                   inset_plot = "auto",
                   inset_type = "boxplot")

# ============================================================
# Panel B: Standardized distance scatter
# ============================================================
z_ref <- (d_ref_all - mean(d_ref_all)) / sd(d_ref_all)
z_emb <- (d_emb_all - mean(d_emb_all)) / sd(d_emb_all)
diff_z <- z_emb - z_ref

df_scatter <- data.frame(
  PCA   = z_ref,
  DNN   = z_emb,
  group = ifelse(diff_z > 0, "DNN > PCA", "DNN <= PCA"),
  stringsAsFactors = FALSE
)
z_range  <- range(c(z_ref, z_emb))
z_lim    <- c(floor(z_range[1]), ceiling(z_range[2]))
z_breaks <- seq(z_lim[1], z_lim[2], by = 1)

p2 <- UtilsR::PlotScatter2(data            = df_scatter,
                   control_col     = "PCA",
                   treat_col       = "DNN",
                   group_col       = "group",
                   group_levels    = c("DNN <= PCA", "DNN > PCA"),
                   highlight_color = "#42B540FF",
                   point_shape     = 21,
                   point_size      = 1.5,
                   point_stroke    = 0.5,
                   x_label         = "Standardized Distance (PCA)",
                   y_label         = "Standardized Distance (DNN)",
                   axis_limits     = z_lim,
                   axis_breaks     = z_breaks,
                   show_ref_lines  = TRUE,
                   test_method     = "wilcox",
                   p_digits        = 2,
                   p_size          = 4,
                   p_x             = z_lim[1] + diff(z_lim) * 0.05,
                   p_y             = z_lim[2] * 0.95,
                   cor_method      = "pearson",
                   cor_size        = 4,
                   hist_bins       = 30,
                   hist_x_limits   = z_lim * 1.1,
                   hist_x_breaks   = z_breaks,
                   median_label_size = 5,
                   inset_left      = 0.6,
                   inset_bottom    = 0.6,
                   inset_right     = 1.1,
                   inset_top       = 1.1,
                   inset_angle     = -45,
                   inset_vp_width  = 0.85,
                   inset_vp_height = 0.75,
                   theme_use       = UtilsR::theme_my(border = FALSE),
                   plot_margin     = c(80, 60, 5, 5),
                   filename        = NULL)

# ============================================================
# Panel C: PCA vs DNN UMAP comparison
# ============================================================
celltype_colors <- UtilsR::pal_paraSC

df_pca     <- umap_data$df_pca |>
  dplyr::mutate(celltype = factor(celltype, levels = names(celltype_colors)))
df_dnn_emb <- umap_data$df_dnn_emb |>
  dplyr::mutate(celltype = factor(celltype, levels = names(celltype_colors)))

p_umap_pca <- scMMR::DimPlot2(df_pca,
                       group.by = "celltype",
                       palcolor = celltype_colors,
                       cells.highlight = TRUE,
                       sizes.highlight = 1,
                       label = TRUE, label_insitu = TRUE, label.size = 2.5,
                       title = "PCA-based UMAP (Seurat)")

p_umap_emb <- scMMR::DimPlot2(df_dnn_emb,
                       group.by = "celltype",
                       palcolor = celltype_colors,
                       cells.highlight = TRUE,
                       sizes.highlight = 1,
                       label = TRUE, label_insitu = TRUE, label.size = 2.5,
                       title = "DNN-based UMAP (scMMR)")

p3 <- p_umap_pca + p_umap_emb + patchwork::plot_layout(ncol = 2, guides = "collect")

# --- 组合与保存 ---
fig4 <- (p1 + p2) / p3
ggplot2::ggsave("./figures/Fig5_EvalEmbedding.pdf", fig4,
                width = 14, height = 10, device = cairo_pdf, dpi = 300)

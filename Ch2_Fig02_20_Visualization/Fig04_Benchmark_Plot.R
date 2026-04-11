# ============================================================
# Fig04_Benchmark_Plot.R
# Figure 4A: ComplexHeatmap (Accuracy + Macro F1)
# Figure 4B: funkyheatmap 综合评估
# 读取 Fig04_Benchmark.R 输出的汇总数据进行可视化
# ============================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

out_dir <- "./data/03_benchmark"

# --- 方法分类与配色 ---
type_map <- c(
  "scMMR" = "Deep Learning", "XGBoost" = "Machine Learning",
  "Random Forest" = "Machine Learning", "SVM" = "Machine Learning",
  "Elastic Net" = "Machine Learning", "LDA" = "Machine Learning",
  "Naive Bayes" = "Machine Learning", "KNN" = "Machine Learning",
  "NNET" = "Machine Learning", "Seurat LT" = "Single-cell",
  "SingleR" = "Single-cell", "CellTypist" = "Single-cell")
type_colors <- c("Deep Learning" = "#00468B",
                 "Machine Learning" = "#42B540",
                 "Single-cell" = "#0099B4")

# --- 读取汇总数据并整理 ---
results_mean <- read.csv(file.path(out_dir, "benchmark_10x10cv_mean.csv"))

# 方法/数据集重命名
results_mean <- results_mean %>%
  dplyr::mutate(Method = dplyr::recode(Method,
    "DNN" = "scMMR", "RF" = "Random Forest", "ElasticNet" = "Elastic Net",
    "NaiveBayes" = "Naive Bayes", "NNet" = "NNET", "Seurat_LT" = "Seurat LT"),
  Dataset = dplyr::recode(Dataset,
    "Baron" = "scRNAseq_Baron", "Muraro" = "scRNAseq_Muraro",
    "Segerstolpe" = "scRNAseq_Segerstolpe", "Xin" = "scRNAseq_Xin",
    "seu_para" = "ParaDataset", "seu_GSE190773" = "Para_GSE190773",
    "seu_GSE233962" = "Para_GSE233962"))

# 固定顺序
ds_order <- c("ParaDataset", "Para_GSE190773", "Para_GSE233962",
              "scRNAseq_Baron", "scRNAseq_Muraro",
              "scRNAseq_Segerstolpe", "scRNAseq_Xin")
method_order <- c("scMMR", "XGBoost", "Random Forest", "SVM",
                  "Elastic Net", "LDA", "Naive Bayes", "KNN",
                  "NNET", "Seurat LT", "SingleR", "CellTypist")

# --- Fig4A: ComplexHeatmap (Accuracy + Macro F1) ---
# 构建矩阵
make_matrix <- function(metric_name) {
  results_mean %>%
    dplyr::filter(Metric == metric_name) %>%
    dplyr::select(Dataset, Method, Mean) %>%
    tidyr::pivot_wider(names_from = Method, values_from = Mean) %>%
    column_to_rownames("Dataset") %>%
    as.matrix()
}

mat_acc <- make_matrix("Accuracy")[ds_order, method_order]
mat_f1  <- make_matrix("Macro_F1")[ds_order, method_order]

# 颜色映射
col_fun <- colorRamp2(c(0.84, 1), c("#F7F7F7", "#B2182B"))
type_vec <- type_map[colnames(mat_acc)]

# 顶部柱形图注释
ha_acc <- HeatmapAnnotation(
  Mean = anno_barplot(colMeans(mat_acc, na.rm = TRUE),
    gp = gpar(fill = type_colors[type_vec]),
    ylim = c(0.85, 0.98), baseline = 0.85,
    bar_width = 0.85, height = unit(4, "cm")),
  annotation_name_side = "left")

ha_f1 <- HeatmapAnnotation(
  Mean = anno_barplot(colMeans(mat_f1, na.rm = TRUE),
    gp = gpar(fill = type_colors[type_vec]),
    ylim = c(0.85, 0.98), baseline = 0.85,
    bar_width = 0.85, height = unit(4, "cm"), axis = FALSE),
  show_annotation_name = FALSE)

# 热图
ht_acc <- Heatmap(mat_acc,
  name = "Value", column_title = "Accuracy",
  top_annotation = ha_acc,
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = col_fun, row_names_side = "left", column_names_rot = 45,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.3f", mat_acc[i, j]), x, y, gp = gpar(fontsize = 8))
  })

ht_f1 <- Heatmap(mat_f1,
  name = "Value2", column_title = "Macro F1",
  top_annotation = ha_f1,
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = col_fun, show_heatmap_legend = FALSE, show_row_names = FALSE,
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.3f", mat_f1[i, j]), x, y, gp = gpar(fontsize = 8))
  })

lgd_type <- Legend(labels = names(type_colors),
                   legend_gp = gpar(fill = type_colors),
                   title = "Method Type", nrow = 3)

pdf("./figures/Fig4a_benchmark_heatmap.pdf", width = 14, height = 6.5)
draw(ht_acc + ht_f1, annotation_legend_list = list(lgd_type), merge_legend = TRUE)
dev.off()

# --- Fig4B: funkyheatmap 综合评估 ---
library(funkyheatmap)

funky_data <- read.csv(file.path(out_dir, "funkyheatmap_data_v4.csv"))
funky_data <- funky_data %>%
  dplyr::mutate(Overall = round(0.40*S_annot + 0.25*S_gen + 0.15*S_scale +
                           0.10*S_robust + 0.10*S_usab, 4))

# Top-3标签
label_top3 <- function(x) {
  r <- rank(-x, ties.method = "min", na.last = "keep")
  ifelse(r <= 3, as.character(r), "")
}
funky_data <- funky_data %>%
  dplyr::mutate(Rank = rank(-Overall, ties.method = "min"),
         Rank_Label = ifelse(Rank <= 3, as.character(Rank), ""),
         Acc_Label  = label_top3(Mean_Acc),
         F1_Label   = label_top3(Mean_F1),
         Gen_Label  = label_top3(Cross_F1))

plot_data <- funky_data %>%
  dplyr::select(id, Method, Model, Platform, GPU, Multi_Task, Interpretable,
         Overall, Rank_Label,
         Mean_Acc, Acc_Label, Mean_F1, F1_Label, Mean_BalAcc, Mean_Kappa,
         Cross_Acc, Cross_F1, Cross_BalAcc, Gen_Label,
         Log_Time, Log_Mem, CV_Acc, Min_BalAcc) %>%
  dplyr::filter(!is.na(id))

# 行排序
row_order <- c("DNN", "XGBoost", "RF", "SVM", "ElasticNet",
               "LDA", "NaiveBayes", "KNN", "NNet",
               "Seurat_LT", "SingleR", "CellTypist")
plot_data <- plot_data[match(row_order, plot_data$id), ]
row_info <- data.frame(id = row_order)

# 列信息
column_info <- tribble(
  ~id,              ~name,          ~geom,       ~group,           ~palette,      ~options,
  "Method",         "Method",       "text",      "info",           NA, list(width = 5, hjust = 0),
  "Model",          "Model",        "text",      "info",           NA, list(width = 3),
  "Platform",       "Platform",     "text",      "info",           NA, list(width = 2),
  "GPU",            "GPU",          "text",      "info",           NA, list(width = 1.2, fontface = "bold"),
  "Multi_Task",     "Multi-task",   "text",      "info",           NA, list(width = 2, fontface = "bold"),
  "Interpretable",  "Interpret.",   "text",      "info",           NA, list(width = 2, fontface = "bold"),
  "Overall",        "Overall",      "bar",       "overall",        "overall_pal", list(width = 4, draw_outline = FALSE),
  "Rank_Label",     NA,             "text",      "overall",        NA, list(hjust = 0.1, overlay = TRUE),
  "Mean_Acc",       "Accuracy",     "bar",       "classification", "blues", list(width = 3, draw_outline = FALSE),
  "Acc_Label",      NA,             "text",      "classification", NA, list(hjust = 0.1, overlay = TRUE),
  "Mean_F1",        "Macro F1",     "bar",       "classification", "blues", list(width = 3, draw_outline = FALSE),
  "F1_Label",       NA,             "text",      "classification", NA, list(hjust = 0.1, overlay = TRUE),
  "Mean_BalAcc",    "Bal.Acc",      "bar",       "classification", "blues", list(width = 3, draw_outline = FALSE),
  "Mean_Kappa",     "Kappa",        "bar",       "classification", "blues", list(width = 3, draw_outline = FALSE),
  "Cross_Acc",      "Acc",          "bar",       "generalization", "purples", list(width = 2.5, draw_outline = FALSE),
  "Cross_F1",       "F1",           "bar",       "generalization", "purples", list(width = 2.5, draw_outline = FALSE),
  "Cross_BalAcc",   "BAcc",         "bar",       "generalization", "purples", list(width = 2.5, draw_outline = FALSE),
  "Gen_Label",      NA,             "text",      "generalization", NA, list(hjust = 0.1, overlay = TRUE),
  "Log_Time",       "Time",         "bar",       "scalability",    "greens_rev", list(width = 3),
  "Log_Mem",        "Memory",       "bar",       "scalability",    "greens_rev", list(width = 3),
  "CV_Acc",         "CV",           "funkyrect", "robustness",     "greys", list(),
  "Min_BalAcc",     "Min BAcc",     "funkyrect", "robustness",     "greys", list())

column_groups <- tribble(
  ~group,           ~Category,             ~palette,
  "info",           "Method",              "method_header",
  "overall",        "Overall",             "overall_header",
  "classification", "Annotation Quality",  "class_header",
  "generalization", "Generalization",      "gen_header",
  "scalability",    "Scalability @100K",   "scale_header",
  "robustness",     "Robustness",          "robust_header")

palettes <- list(
  method_header  = c("grey", "grey"),
  overall_header = c("#FC9272", "#FC9272"),
  class_header   = c("#AED6F1", "#AED6F1"),
  gen_header     = c("#D7BDE2", "#D7BDE2"),
  scale_header   = c("#ABEBC6", "#ABEBC6"),
  robust_header  = c("#F5CBA7", "#F5CBA7"),
  overall_pal    = RColorBrewer::brewer.pal(9, "Reds")[1:7],
  blues          = RColorBrewer::brewer.pal(9, "Blues")[1:7],
  purples        = RColorBrewer::brewer.pal(9, "Purples")[1:7],
  greens_rev     = RColorBrewer::brewer.pal(9, "Greens")[1:7],
  greys          = RColorBrewer::brewer.pal(9, "Oranges")[1:7],
  category       = c("Deep Learning" = "#FFFFFF",
                     "Machine Learning" = "#FFFFFF",
                     "Single-cell" = "#FFFFFF"))

legends <- list(
  list(title = "Overall Score", palette = "overall_pal", geom = "bar",
       labels = c("0", "", "0.25", "", "0.5", "", "0.75", "", "1")),
  list(title = "Annotation Quality", palette = "blues", geom = "bar",
       labels = c("0.85", "", "0.90", "", "0.95", "", "1.0")),
  list(title = "Generalization", palette = "purples", geom = "bar",
       labels = c("0.85", "", "0.90", "", "0.95", "", "1.0")),
  list(title = "Scalability (log10)", palette = "greens_rev", geom = "bar",
       labels = c("1x", "", "10x", "", "100x", "", "500x")),
  list(palette = "method_header",  enabled = FALSE),
  list(palette = "overall_header", enabled = FALSE),
  list(palette = "class_header",   enabled = FALSE),
  list(palette = "gen_header",     enabled = FALSE),
  list(palette = "scale_header",   enabled = FALSE),
  list(palette = "robust_header",  enabled = FALSE),
  list(title = "Robustness", palette = "greys", geom = "funkyrect",
       labels = c("Low", "", "", "", "High")))

# 绘图
g <- funky_heatmap(
  data = plot_data, column_info = column_info,
  column_groups = column_groups, row_info = row_info,
  palettes = palettes, legends = legends, scale_column = TRUE,
  position_args = position_arguments(
    col_annot_offset = 3, expand_xmin = 1, expand_xmax = 4,
    expand_ymin = 0.5, expand_ymax = 0.5))

g <- g + ggplot2::theme(text = ggplot2::element_text(family = "Arial"))
ggplot2::ggsave("./figures/Fig4b_funkyheatmap.pdf", g,
                device = cairo_pdf, width = g$width * 0.9, height = g$height)



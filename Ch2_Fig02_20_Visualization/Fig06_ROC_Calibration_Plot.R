# ============================================================
# Fig06_ROC_Calibration.R
# Panel A/C/E: Per-class ROC curves (one-vs-rest)
# Panel B/D/F: Pooled calibration curves
# 3 datasets: Reference / Validation 1 / Validation 2
# ============================================================

library(ggplot2)
library(patchwork)
library(pROC)
library(UtilsR)
library(qs)

# --- 加载数据 (由 Figure5_step1 生成的中间结果) ---
dat <- readRDS("./out/Figure5_intermediate.rds")
list2env(dat, envir = environment())

ct_colors <- UtilsR::pal_paraSC

# --- 主题 ---
theme_fig5 <- ggplot2::theme_classic(base_size = 11) + ggplot2::theme(
  axis.text   = ggplot2::element_text(color = "black", size = 13, face = "bold"),
  axis.title  = ggplot2::element_text(color = "black", size = 14, face = "bold"),
  axis.line   = ggplot2::element_line(linewidth = 0.8, color = "black"),
  axis.ticks  = ggplot2::element_line(linewidth = 0.8, color = "black"),
  plot.title  = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
  plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "grey40"),
  legend.text   = ggplot2::element_text(size = 9),
  legend.title  = ggplot2::element_text(size = 8),
  legend.position = c(0.8, 0.3),
  plot.margin   = ggplot2::margin(10, 10, 10, 10))

# ============================================================
# 辅助函数: plot_roc (per-class one-vs-rest)
# ============================================================
plot_roc <- function(true_labels, probs, title = "ROC", colors = ct_colors) {
  class_names <- colnames(probs)
  present <- intersect(class_names, unique(true_labels))

  pred_labels <- class_names[apply(probs, 1, which.max)]
  acc <- mean(true_labels == pred_labels)
  f1_vals <- sapply(present, function(cls) {
    tp <- sum(pred_labels == cls & true_labels == cls)
    fp <- sum(pred_labels == cls & true_labels != cls)
    fn <- sum(pred_labels != cls & true_labels == cls)
    prec <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    rec  <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
  })
  macro_f1 <- mean(f1_vals, na.rm = TRUE)

  ct_order <- if (is.factor(true_labels)) levels(true_labels) else
    names(sort(table(true_labels), decreasing = TRUE))

  roc_list <- list(); auc_vals <- c()
  for (cls in present) {
    bt <- as.integer(true_labels == cls)
    if (sum(bt) < 2 || sum(bt) == length(bt)) next
    ci <- which(class_names == cls); if (length(ci) == 0) next
    ro <- tryCatch(pROC::roc(bt, probs[, ci], quiet = TRUE), error = function(e) NULL)
    if (!is.null(ro)) { roc_list[[cls]] <- ro; auc_vals[cls] <- as.numeric(pROC::auc(ro)) }
  }
  if (length(roc_list) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())
  macro_auc <- mean(auc_vals, na.rm = TRUE)

  roc_df <- do.call(rbind, lapply(names(roc_list), function(cls) {
    r <- roc_list[[cls]]
    data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities,
               CellType = cls, AUC = auc_vals[cls])
  }))
  roc_df$label <- paste0(roc_df$CellType, " (", sprintf("%.3f", roc_df$AUC), ")")
  uc <- ct_order[ct_order %in% names(auc_vals)]
  label_levels <- paste0(uc, " (", sprintf("%.3f", auc_vals[uc]), ")")
  roc_df$label <- factor(roc_df$label, levels = label_levels)
  pal <- colors[uc]; names(pal) <- label_levels

  ggplot2::ggplot(roc_df, ggplot2::aes(FPR, TPR, color = label)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.4) +
    ggplot2::scale_color_manual(values = pal, name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "FPR (1-Specificity)", y = "TPR (Sensitivity)", title = title,
         subtitle = sprintf("Macro AUC=%.4f | Acc=%.4f | Macro F1=%.4f",
                            macro_auc, acc, macro_f1)) +
    theme_fig5 +
    ggplot2::theme(legend.key.size = unit(0.35, "cm"),
          legend.spacing.y = unit(0.05, "cm"))
}

# ============================================================
# 辅助函数: plot_calibration (pooled)
# ============================================================
plot_calibration <- function(true_labels, probs, title = "Calibration",
                             n_bins = 10, colors = ct_colors) {
  class_names <- colnames(probs)
  present <- intersect(class_names, unique(true_labels))
  breaks <- seq(0, 1, length.out = n_bins + 1)

  pred_labels <- class_names[apply(probs, 1, which.max)]
  acc <- mean(true_labels == pred_labels)
  f1_vals <- sapply(present, function(cls) {
    tp <- sum(pred_labels == cls & true_labels == cls)
    fp <- sum(pred_labels == cls & true_labels != cls)
    fn <- sum(pred_labels != cls & true_labels == cls)
    prec <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    rec  <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
  })
  macro_f1 <- mean(f1_vals, na.rm = TRUE)

  # Pooled: 所有class合并计算Brier和校准曲线
  all_bt <- c(); all_p <- c()
  for (cls in present) {
    ci <- which(class_names == cls); if (length(ci) == 0) next
    all_bt <- c(all_bt, as.integer(true_labels == cls))
    all_p  <- c(all_p, probs[, ci])
  }
  brier_pooled <- mean((all_p - all_bt)^2)
  bins <- cut(all_p, breaks = breaks, include.lowest = TRUE)
  agg <- aggregate(cbind(p = all_p, y = all_bt) ~ bins, FUN = mean)
  agg$n <- as.integer(table(bins)[as.character(agg$bins)])

  ggplot2::ggplot(agg, ggplot2::aes(p, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.4) +
    ggplot2::geom_line(color = "#D55E00", linewidth = 0.8) +
    ggplot2::geom_point(color = "#D55E00", size = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "Mean Predicted Probability", y = "Observed Frequency",
         title = title,
         subtitle = sprintf("Brier=%.4f | Acc=%.4f | Macro F1=%.4f",
                            brier_pooled, acc, macro_f1)) +
    theme_fig5
}

# --- 绘图 ---
p_roc_ref <- plot_roc(true_ref, probs_ref, "Reference (n=170,852)")
p_cal_ref <- plot_calibration(true_ref, probs_ref, "Reference (n=170,852)")

p_roc_v1 <- plot_roc(true_v1, probs_v1, "Validation 1 (n=70,000, top-2 noise)")
p_cal_v1 <- plot_calibration(true_v1, probs_v1, "Validation 1 (n=70,000, top-2 noise)")

p_roc_v2 <- plot_roc(true_v2, probs_v2, "Validation 2 (n=40,000, top-2 noise)")
p_cal_v2 <- plot_calibration(true_v2, probs_v2, "Validation 2 (n=40,000, top-2 noise)")

# --- 组合与保存 ---
fig5 <- patchwork::wrap_plots(p_roc_ref, p_cal_ref,
                              p_roc_v1,  p_cal_v1,
                              p_roc_v2,  p_cal_v2, ncol = 2)
ggplot2::ggsave("./figures/Fig6_ROC_Calibration.pdf", fig5,
                width = 16, height = 16, device = pdf, limitsize = FALSE)

# ============================================================
# CH4_Fig29_Plot.R
# PHPT vs Normal TWEAK 通路通讯对比可视化
# 含自定义函数: PlotCellChatHeatmap / Bubble / Circle / Violin
# 输入: ./out/cellchat/cellchat_PH_subcluster.qs
#       ./out/cellchat/cellchat_Normal_subcluster.qs
# 输出: ./figures/CH4_Fig29_TWEAK_Circle_Comparison.pdf
# ============================================================

library(CellChat)
library(ComplexHeatmap)
library(circlize)
library(UtilsR)
library(qs)
library(ggplot2)
library(patchwork)

out_dir <- "./out/cellchat"
fig_dir <- "./figures"

# ============================================================
# PlotCellChatHeatmap: 通路通讯热力图 (单个或两组对比, 统一色阶)
# ============================================================

PlotCellChatHeatmap <- function(
    cellchat1, cellchat2 = NULL, signaling = "TWEAK",
    label1 = "Group1", label2 = "Group2",
    cell_colors = NULL, cell_order = NULL,
    color.heatmap = "Reds", row_names_side = "left",
    fontsize = 8, bar_height = 1.2, bar_width = 1.2,
    anno_size = 0.2, as_ggplot = TRUE,
    filename = NULL, width = 12, height = 9) {

  mat1 <- cellchat1@netP$prob[,, signaling]
  two_groups <- !is.null(cellchat2)
  if (two_groups) {
    mat2 <- cellchat2@netP$prob[,, signaling]
    all_cells <- sort(union(rownames(mat1), rownames(mat2)))
  } else {
    all_cells <- sort(rownames(mat1))
  }
  if (!is.null(cell_order)) {
    missing <- setdiff(all_cells, cell_order)
    all_cells <- c(cell_order[cell_order %in% all_cells], missing)
  }

  expand_mat <- function(mat, cells) {
    m <- matrix(0, length(cells), length(cells),
                dimnames = list(cells, cells))
    shared <- intersect(rownames(mat), cells)
    m[shared, shared] <- mat[shared, shared]; m
  }
  mat1 <- expand_mat(mat1, all_cells)
  if (two_groups) mat2 <- expand_mat(mat2, all_cells)

  global_max <- if (two_groups) max(mat1, mat2) else max(mat1)
  heatmap_cols <- switch(color.heatmap,
    "Reds" = c("white", "red"), "Blues" = c("white", "blue"),
    "Greens" = c("white", "darkgreen"), c("white", "red"))
  col_fun <- colorRamp2(c(0, global_max), heatmap_cols)

  bar_max_col <- if (two_groups) max(colSums(mat1), colSums(mat2)) else max(colSums(mat1))
  bar_max_row <- if (two_groups) max(rowSums(mat1), rowSums(mat2)) else max(rowSums(mat1))

  if (is.null(cell_colors)) {
    pal <- UtilsR::pal_lancet
    cell_colors <- setNames(rep_len(pal, length(all_cells)), all_cells)
  }

  build_heatmap <- function(mat, title, heatmap_name, show_legend = TRUE) {
    top_anno <- HeatmapAnnotation(
      bar = anno_barplot(colSums(mat), height = unit(bar_height, "cm"),
                         ylim = c(0, bar_max_col),
                         gp = gpar(fill = cell_colors[all_cells], col = NA)),
      show_annotation_name = FALSE)
    right_anno <- rowAnnotation(
      bar = anno_barplot(rowSums(mat), width = unit(bar_width, "cm"),
                         ylim = c(0, bar_max_row),
                         gp = gpar(fill = cell_colors[all_cells], col = NA)),
      show_annotation_name = FALSE)
    row_anno <- rowAnnotation(
      celltype = all_cells, col = list(celltype = cell_colors),
      show_legend = FALSE, show_annotation_name = FALSE,
      simple_anno_size = unit(anno_size, "cm"))
    bottom_anno <- HeatmapAnnotation(
      celltype = all_cells, col = list(celltype = cell_colors),
      show_legend = FALSE, show_annotation_name = FALSE,
      simple_anno_size = unit(anno_size, "cm"))
    Heatmap(mat, name = heatmap_name, col = col_fun,
            show_heatmap_legend = show_legend,
            column_title = paste0(title, " - ", signaling, " signaling network"),
            row_title = "Sources (Sender)",
            left_annotation = row_anno, right_annotation = right_anno,
            top_annotation = top_anno, bottom_annotation = bottom_anno,
            cluster_rows = FALSE, cluster_columns = FALSE,
            row_names_side = row_names_side,
            row_names_gp = gpar(fontsize = fontsize),
            column_names_gp = gpar(fontsize = fontsize),
            column_names_rot = 45,
            column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  }

  ht1 <- build_heatmap(mat1, label1, "Comm.\nProb.")
  if (two_groups) {
    ht2 <- build_heatmap(mat2, label2, "Comm.\nProb.2", show_legend = FALSE)
    ht_list <- ht1 + ht2
    draw_expr <- function() draw(ht_list, ht_gap = unit(1.5, "cm"),
      heatmap_legend_side = "right",
      column_title = paste0(signaling, " Signaling: ", label1, " vs ", label2,
                            " (unified color scale)"),
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      merge_legend = TRUE)
  } else {
    draw_expr <- function() draw(ht1)
  }
  if (!is.null(filename)) { pdf(filename, width = width, height = height); draw_expr(); dev.off() }
  if (as_ggplot) ggplotify::as.ggplot(function() draw_expr()) else ht1
}

# ============================================================
# PlotCellChatBubble: 通讯气泡图 (两组对比, 支持 facet by target)
# ============================================================

PlotCellChatBubble <- function(
    cellchat1, cellchat2, signaling = "TWEAK",
    label1 = "Group1", label2 = "Group2",
    target_pattern = NULL, target_levels = NULL,
    exclude_source = NULL, colors = NULL,
    point_size = c(2, 9), base_size = 11,
    title = NULL, filename = NULL, width = 8, height = 14) {

  df1 <- as.data.frame(subsetCommunication(cellchat1, signaling = signaling))
  df2 <- as.data.frame(subsetCommunication(cellchat2, signaling = signaling))
  df1$group <- label1; df2$group <- label2
  df_all <- rbind(df1, df2)
  df_all$pair <- paste0(df_all$source, " -> ", df_all$target)
  if (is.null(colors)) colors <- c("#2166AC", "#67A9CF", "#FDDBC7", "#EF8A62", "#B2182B")

  if (is.null(target_pattern)) {
    df_all$group <- factor(df_all$group, levels = c(label1, label2))
    df_all <- df_all[order(df_all$target, -df_all$prob), ]
    df_all$pair <- factor(df_all$pair, levels = rev(unique(df_all$pair)))
    if (is.null(title)) title <- paste0(signaling, ": ", label1, " vs ", label2)
    p <- ggplot(df_all, aes(x = group, y = pair)) +
      geom_point(aes(size = prob, color = prob)) +
      scale_color_gradientn(colors = colors, name = "Comm. Prob.") +
      scale_size_continuous(range = point_size, name = "Comm. Prob.") +
      theme_bw(base_size = base_size) + labs(title = title, x = "", y = "")
  } else {
    df_sub <- df_all[grepl(target_pattern, df_all$target), ]
    if (!is.null(exclude_source)) df_sub <- df_sub[!df_sub$source %in% exclude_source, ]
    df_sub$group <- factor(df_sub$group, levels = c(label1, label2))
    if (!is.null(target_levels)) df_sub$target <- factor(df_sub$target, levels = target_levels)
    df_sub <- df_sub[order(df_sub$target, -df_sub$prob), ]
    df_sub$source <- factor(df_sub$source, levels = rev(unique(df_sub$source)))
    if (is.null(title)) title <- paste0(signaling, " -> ", target_pattern, ": ", label1, " vs ", label2)
    p <- ggplot(df_sub, aes(x = group, y = source)) +
      geom_point(aes(size = prob, color = prob)) +
      scale_color_gradientn(colors = colors, name = "Comm. Prob.") +
      scale_size_continuous(range = point_size, name = "Comm. Prob.") +
      facet_grid(target ~ ., scales = "free_y", space = "free_y") +
      theme_bw(base_size = base_size) +
      theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) +
      labs(title = title, x = "", y = "Sender cell type")
  }
  if (!is.null(filename)) ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
# PlotCellChatCircle: 圆形图 (统一线宽映射)
# ============================================================

PlotCellChatCircle <- function(
    cellchat1, cellchat2 = NULL, signaling = "TWEAK",
    label1 = "Group1", label2 = "Group2",
    cell_colors = NULL, cell_order = NULL,
    edge.width.max = 10, vertex.label.cex = 0.8,
    arrow.size = 0.3, layout = "circle",
    filename = NULL, width = 14, height = 7) {

  two_groups <- !is.null(cellchat2)
  mat1 <- cellchat1@netP$prob[,, signaling]
  global_max <- if (two_groups) max(mat1, cellchat2@netP$prob[,, signaling]) else max(mat1)

  all_idents <- levels(cellchat1@idents)
  if (two_groups) all_idents <- union(all_idents, levels(cellchat2@idents))
  if (is.null(cell_colors)) {
    pal <- UtilsR::pal_lancet
    cell_colors <- setNames(rep_len(pal, length(all_idents)), all_idents)
  }

  reorder_idents <- function(cc, ord) {
    current <- levels(cc@idents)
    new_levels <- c(ord[ord %in% current], setdiff(current, ord))
    cc@idents <- factor(cc@idents, levels = new_levels)
    for (sn in c("prob", "pval")) {
      arr <- cc@net[[sn]]; common <- intersect(new_levels, dimnames(arr)[[1]])
      cc@net[[sn]] <- arr[common, common, , drop = FALSE]
    }; cc
  }
  if (!is.null(cell_order)) {
    cellchat1 <- reorder_idents(cellchat1, cell_order)
    if (two_groups) cellchat2 <- reorder_idents(cellchat2, cell_order)
  }

  if (!is.null(filename)) pdf(filename, width = width, height = height)
  if (two_groups) par(mfrow = c(1, 2), xpd = TRUE) else par(mfrow = c(1, 1), xpd = TRUE)

  cols1 <- cell_colors[levels(cellchat1@idents)]
  netVisual_aggregate(cellchat1, signaling = signaling, layout = layout,
    edge.weight.max = global_max, edge.width.max = edge.width.max,
    vertex.label.cex = vertex.label.cex, arrow.size = arrow.size, color.use = cols1)
  title(paste0(label1, " - ", signaling), cex.main = 1.5)

  if (two_groups) {
    cols2 <- cell_colors[levels(cellchat2@idents)]
    netVisual_aggregate(cellchat2, signaling = signaling, layout = layout,
      edge.weight.max = global_max, edge.width.max = edge.width.max,
      vertex.label.cex = vertex.label.cex, arrow.size = arrow.size, color.use = cols2)
    title(paste0(label2, " - ", signaling), cex.main = 1.5)
  }
  if (!is.null(filename)) dev.off()
}

# ============================================================
# PlotCellChatViolin: 配受体基因表达小提琴图 (统一 x 轴)
# ============================================================

PlotCellChatViolin <- function(
    cellchat1, cellchat2 = NULL, signaling = "TWEAK",
    label1 = "Group1", label2 = "Group2",
    cell_order = NULL, cell_colors = NULL,
    type = "violin", x_angle = 45, shared_x = FALSE,
    show_color_bar = FALSE, color_bar_height = 0.1,
    filename = NULL, width = 12, height = 10) {

  two_groups <- !is.null(cellchat2)
  levels1 <- levels(cellchat1@idents)
  all_cells <- if (two_groups) sort(union(levels1, levels(cellchat2@idents))) else levels1
  if (!is.null(cell_order)) {
    all_cells <- c(cell_order[cell_order %in% all_cells], setdiff(all_cells, cell_order))
  }
  cellchat1@idents <- factor(cellchat1@idents, levels = all_cells)
  if (two_groups) cellchat2@idents <- factor(cellchat2@idents, levels = all_cells)

  if (is.null(cell_colors)) {
    pal <- UtilsR::pal_lancet
    cell_colors <- setNames(rep_len(pal, length(all_cells)), all_cells)
  }
  color_use <- cell_colors[all_cells]

  p1 <- plotGeneExpression(cellchat1, signaling = signaling, type = type,
                           color.use = color_use)
  p1 <- p1 & scale_x_discrete(drop = FALSE, limits = all_cells)

  if (two_groups) {
    p2 <- plotGeneExpression(cellchat2, signaling = signaling, type = type,
                             color.use = color_use)
    p2 <- p2 & scale_x_discrete(drop = FALSE, limits = all_cells)
  }

  if (shared_x && two_groups) {
    plots1 <- p1$patches$plots; plots2 <- p2$patches$plots
    all_plots <- c(plots1, plots2); n <- length(all_plots)
    all_plots[[1]] <- all_plots[[1]] + ggtitle(label1) +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
    all_plots[[length(plots1) + 1]] <- all_plots[[length(plots1) + 1]] + ggtitle(label2) +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
    for (i in seq_len(n)) {
      all_plots[[i]] <- all_plots[[i]] +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    if (show_color_bar) {
      df_bar <- data.frame(x = factor(all_cells, levels = all_cells), y = 1)
      hjust <- if (!is.null(x_angle) && x_angle > 0) 1 else 0
      p_bar <- ggplot(df_bar, aes(x = x, y = y, fill = x)) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_manual(values = color_use, guide = "none") +
        scale_x_discrete(drop = FALSE, limits = all_cells) +
        theme_void() +
        theme(axis.text.x = element_text(angle = x_angle, hjust = hjust, vjust = 1, size = 9))
      p_out <- wrap_plots(c(all_plots, list(p_bar)), ncol = 1,
                          heights = c(rep(1, n), color_bar_height))
    } else {
      if (!is.null(x_angle) && x_angle != 0) {
        hjust <- if (x_angle > 0) 1 else 0
        all_plots[[n]] <- all_plots[[n]] +
          theme(axis.text.x = element_text(angle = x_angle, hjust = hjust, vjust = 1))
      }
      p_out <- wrap_plots(all_plots, ncol = 1)
    }
  } else {
    p1 <- p1 + plot_annotation(title = label1)
    if (two_groups) {
      p2 <- p2 + plot_annotation(title = label2)
      p_out <- wrap_elements(p1) / wrap_elements(p2)
    } else { p_out <- p1 }
  }

  if (!is.null(filename)) ggsave(filename, p_out, width = width, height = height)
  p_out
}

# ============================================================
# 5. 加载数据 & 调色板
# ============================================================

cellchat_PH     <- qs::qread(file.path(out_dir, "cellchat_PH_subcluster.qs"))
cellchat_Normal <- qs::qread(file.path(out_dir, "cellchat_Normal_subcluster.qs"))

# --- 项目专用调色板 (Para subclusters 版) ---
pal_subcluster <- c(
  "Para_cluster0"          = "#E8725C",
  "Para_cluster1"          = "#C43C3C",
  "Para_cluster3"          = "#8B1A1A",
  "T cells"                = "#56B4E9",
  "Cycling T cells"        = "#8B5E3C",
  "NK cells"               = "#2E8B45",
  "B cells"                = "#1A7B7B",
  "Monocytes"              = "#7B68AA",
  "M1-like Macrophages"    = "#D4919A",
  "M2/M3-like Macrophages" = "#F0C8A0",
  "Neutrophils"            = "#C0C0C0",
  "cDC2s"                  = "#282828",
  "Mast cells"             = "#F0E442",
  "iTAFs"                  = "#9AB83C",
  "mTAFs"                  = "#E07B1A",
  "Pericytes"              = "#6B3FA0",
  "Capillary ECs"          = "#D4C8E8",
  "Venous ECs"             = "#8AAAC8"
)

# ============================================================
# 6. Circle: 统一线宽对比 → p1 (base graphics → ggplot)
# ============================================================

p1 <- ggplotify::as.ggplot(function() {
  par(mfrow = c(1, 2), xpd = TRUE)
  mat1 <- cellchat_PH@netP$prob[,, "TWEAK"]
  mat2 <- cellchat_Normal@netP$prob[,, "TWEAK"]
  global_max <- max(mat1, mat2)

  reorder_idents <- function(cc, ord) {
    current <- levels(cc@idents)
    new_levels <- c(ord[ord %in% current], setdiff(current, ord))
    cc@idents <- factor(cc@idents, levels = new_levels)
    for (sn in c("prob", "pval")) {
      arr <- cc@net[[sn]]; common <- intersect(new_levels, dimnames(arr)[[1]])
      cc@net[[sn]] <- arr[common, common, , drop = FALSE]
    }; cc
  }
  cc_ph <- reorder_idents(cellchat_PH, names(pal_subcluster))
  cc_nm <- reorder_idents(cellchat_Normal, names(pal_subcluster))

  cols1 <- pal_subcluster[levels(cc_ph@idents)]
  netVisual_aggregate(cc_ph, signaling = "TWEAK", layout = "circle",
    edge.weight.max = global_max, edge.width.max = 10,
    vertex.label.cex = 0.8, arrow.size = 0.3, color.use = cols1)
  title("PHPT - TWEAK", cex.main = 1.5)

  cols2 <- pal_subcluster[levels(cc_nm@idents)]
  netVisual_aggregate(cc_nm, signaling = "TWEAK", layout = "circle",
    edge.weight.max = global_max, edge.width.max = 10,
    vertex.label.cex = 0.8, arrow.size = 0.3, color.use = cols2)
  title("Normal - TWEAK", cex.main = 1.5)
})

# ============================================================
# 7. Bubble: 仅 -> Para subclusters (facet) → p2
# ============================================================

p2 <- PlotCellChatBubble(
  cellchat_PH, cellchat_Normal,
  signaling      = "TWEAK",
  label1         = "PHPT",
  label2         = "Normal",
  target_pattern = "Para_",
  target_levels  = c("Para_cluster3", "Para_cluster1", "Para_cluster0"),
  exclude_source = c("Para_cluster0", "Para_cluster1", "Para_cluster3"),
  point_size     = c(3, 11),
  base_size      = 12,
  width          = 12,
  height         = 5
)

# ============================================================
# 9. Violin: 统一 x 轴基因表达 → p3
# ============================================================

p3 <- PlotCellChatViolin(
  cellchat_PH, cellchat_Normal,
  signaling      = "TWEAK",
  label1         = "PHPT",
  label2         = "Normal",
  cell_colors    = pal_subcluster,
  cell_order     = names(pal_subcluster),
  shared_x       = TRUE,
  show_color_bar = TRUE,
  width          = 14,
  height         = 10
)

# ============================================================
# 10. 组合 Circle + Bubble + Violin → 保存
# ============================================================

p_combined <- wrap_plots(p1, p2, p3, ncol = 1,
                         heights = c(1, 0.6, 1.2))
ggsave(file.path(fig_dir, "CH4_Fig29_TWEAK_Comparison.pdf"),
       p_combined, width = 14, height = 22)


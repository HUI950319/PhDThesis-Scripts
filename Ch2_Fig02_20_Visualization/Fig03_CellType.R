# ============================================================
# Fig03_CellType.R
# Harmony整合 → 聚类 → 细胞类型注释 → UMAP/热图
# Figure 3: A-分组UMAP / B-Marker热图
# ============================================================

library(scMMR)
library(UtilsR)
library(Seurat)
library(tidyverse)

# 细胞类型配色 (16 types, Capillary ECs I/II合并)
my_color <- UtilsR::pal_paraSC

# --- 加载质控后数据 ---
seu <- qs::qread("./data/seu_QC.qs")

seu@meta.data <- seu@meta.data %>%
  dplyr::mutate(sample = orig.ident,
         group = dplyr::case_when(
           sample %in% paste0("SH", 1:6) ~ "SH",
           sample %in% paste0("PT", 1:3) ~ "PT",
           sample %in% paste0("PH", 1:6) ~ "PH"))

# --- 标准预处理 + Harmony整合 ---
seu <- Seurat::NormalizeData(seu)
seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000)
seu <- Seurat::ScaleData(seu)
seu <- Seurat::RunPCA(seu)
seu <- harmony::RunHarmony(seu, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:50)
seu <- Seurat::RunUMAP(seu, reduction = "harmony", dims = 1:30)

# 聚类 (测试多个resolution，最终选择0.2)
seu <- Seurat::FindNeighbors(seu, reduction = "harmony", dims = 1:30, k.param = 20)
resolutions <- seq(0.1, 0.6, 0.1)
seu <- Seurat::FindClusters(seu, resolution = resolutions)
scMMR::DimPlot2(seu, reduction = "umap",
                  group.by = paste0("RNA_snn_res.", resolutions),
                  ncol = 2, label = TRUE) & NoLegend()
ggplot2::ggsave("./figures/Fig3_resolution_test.pdf", width = 10, height = 12)

# 选择 resolution = 0.2
seu$seurat_clusters <- seu$RNA_snn_res.0.2
Seurat::Idents(seu) <- seu$seurat_clusters

# --- 细胞类型注释 (基于7篇甲状旁腺文献marker汇总) ---
seu <- RenameIdents(seu,
                    "0"  = "Parathyroid cells",
                    "1"  = "M1-like Macrophages",
                    "2"  = "Neutrophils",
                    "3"  = "Cycling T cells",
                    "4"  = "Monocytes",
                    "5"  = "T cells",
                    "6"  = "NK cells",
                    "7"  = "iTAFs",
                    "8"  = "Capillary ECs I",
                    "9"  = "M2/M3-like Macrophages",
                    "10" = "cDC2s",
                    "11" = "Mast cells",
                    "12" = "mTAFs",
                    "13" = "Venous ECs",
                    "14" = "B cells",
                    "15" = "Pericytes",
                    "16" = "Capillary ECs II")
seu$celltype <- Seurat::Idents(seu)

# 合并 Capillary ECs I/II
seu$celltype <- dplyr::recode(seu$celltype,
                              "Capillary ECs I"  = "Capillary ECs",
                              "Capillary ECs II" = "Capillary ECs")

seu$celltype <- factor(seu$celltype, levels = c(
  "Parathyroid cells",
  "T cells", "Cycling T cells", "NK cells", "B cells",
  "Monocytes", "M1-like Macrophages", "M2/M3-like Macrophages",
  "Neutrophils", "cDC2s", "Mast cells",
  "iTAFs", "mTAFs", "Pericytes",
  "Capillary ECs", "Venous ECs"))

seu$group <- dplyr::recode(seu$group, "PH" = "PHPT", "PT" = "Normal", "SH" = "SHPT")
seu$group <- factor(seu$group, levels = c("PHPT", "Normal", "SHPT"))
Seurat::Idents(seu) <- seu$celltype

qs::qsave(seu, "./data/seu_annotated.qs")

# --- Fig3A: UMAP按分组拆分 (主面板) ---
scMMR::DimPlot2(seu, reduction = "umap", group.by = "celltype", split.by = "group",
                  palcolor = my_color, theme_use = "theme_blank",
                  cells.highlight = TRUE, sizes.highlight = 2) +
  patchwork::plot_layout(guides = "collect")
ggplot2::ggsave("./figures/Fig3a_umap_by_group.pdf", width = 12, height = 6)

# --- Fig3B: GroupHeatmap (marker基因热图) ---
markers_df <- data.frame(
  marker = c(
    "EPCAM",   "CASR",
    "CD3D",    "CD3E",
    "STMN1",   "HMGB2",
    "GNLY",    "NKG7",
    "MS4A1",   "CD79A",
    "LYZ",     "FCN1",
    "CCL3",    "RGS1",
    "APOE",    "C1QA",
    "S100A12", "NAMPT",
    "CLEC10A", "CD1C",
    "TPSAB1",  "CPA3",
    "PDGFRA",  "CXCL12",
    "MYH11",   "TAGLN",
    "RGS5",    "PDGFRB",
    "PLVAP",   "CA4",
    "ACKR1",   "SELP"),
  celltype = rep(levels(seu$celltype), each = 2))

seu_sub <- subset(seu, downsample = 500)

scop::GroupHeatmap(
  seu_sub,
  features = markers_df$marker,
  feature_split = markers_df$celltype,
  group.by = "celltype",
  group_palcolor = my_color,
  feature_split_palcolor = my_color,
  heatmap_palette = "YlOrRd",
  cell_annotation = c("group", "Phase", "percent.mt"),
  cell_annotation_palette = c("Dark2", "Chinese", "simspec"),
  cell_annotation_params = list(height = grid::unit(10, "mm")),
  add_dot = TRUE, dot_size = unit(6, "mm"),
  add_bg = TRUE, nlabel = 0, show_row_names = TRUE)
ggplot2::ggsave("./figures/Fig3b_heatmap.pdf", width = 14.5, height = 8)


# ============================================================
# Fig02_QC.R
# 单细胞质量控制: QC指标 → 过滤 → 双细胞检测 → 游离RNA去污
# Figure 2: A-QC小提琴图 / B-Doublet UMAP / C-DecontX UMAP
# ============================================================

library(scMMR)
library(UtilsR)
library(Seurat)
library(tidyverse)

pal_lan <- as.character(UtilsR::pal_lancet)

# --- 导入15个样本原始矩阵并合并 ---
samples <- list.dirs(path = "./data/raw_mat", full.names = FALSE, recursive = FALSE)
seu.list <- pbapply::pblapply(samples, function(sn) {
  counts <- Seurat::Read10X(file.path("./data/raw_mat", sn))
  sn <- gsub("_", "-", sn)
  colnames(counts) <- paste(sn, colnames(counts), sep = "_")
  Seurat::CreateSeuratObject(counts = counts)
}, cl = length(samples))

seu <- merge(seu.list[[1]], y = seu.list[-1])
seu <- Seurat::JoinLayers(seu, assay = "RNA")
rm(seu.list); gc()

# --- 计算QC指标 ---
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
ribo.genes <- intersect(scMMR::ribo.genes, rownames(seu))
seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, features = ribo.genes)

# 添加分组信息
seu@meta.data <- seu@meta.data %>%
  dplyr::mutate(sample = orig.ident,
         group = dplyr::case_when(
           sample %in% paste0("SH", 1:6) ~ "SH",
           sample %in% paste0("PT", 1:3) ~ "PT",
           sample %in% paste0("PH", 1:6) ~ "PH"))

seu$sample <- factor(seu$sample,
                     levels = c(paste0("PH", 1:6), paste0("PT", 1:3), paste0("SH", 1:6)),
                     labels = c(paste0("PHPT", 1:6), paste0("Normal", 1:3), paste0("SHPT", 1:6)))

# --- Fig2A: QC小提琴图 ---
pal_grp <- c(rep(pal_lan[1], 6), rep(pal_lan[2], 3), rep(pal_lan[3], 6))

p1 <- scMMR::FeaturePlot3(
  seu, group.by = "sample", bg.by = "group", features = "nFeature_RNA",
  add_box = TRUE, stack = TRUE, palcolor = pal_grp,
  bg_palette = "lancet", y.min = 0, y.max = 14000, y.nbreaks = 3,
  legend.position = "none")

p2 <- scMMR::FeaturePlot3(
  seu, group.by = "sample", bg.by = "group", features = "nCount_RNA",
  add_box = TRUE, stack = TRUE, palcolor = pal_grp,
  bg_palette = "lancet", y.min = 0, y.max = 60000, y.nbreaks = 3,
  legend.position = "none")

p3 <- scMMR::FeaturePlot3(
  seu, group.by = "sample", bg.by = "group", features = "percent.mt",
  add_box = TRUE, stack = TRUE, palcolor = pal_grp,
  bg_palette = "lancet", y.min = 0, y.max = 100, y.nbreaks = 3,
  legend.position = "none")

pA <- (p1 %>% UtilsR::fmt_ref(y = c(200, 10000), color = "red")) /
  (p2 %>% UtilsR::fmt_ref(y = c(1000, 50000), color = "red")) /
  (p3 %>% UtilsR::fmt_ref(y = 40, color = "red"))

# --- 质控过滤 ---
seu$QC <- ifelse(seu$nFeature_RNA < 200 | seu$nFeature_RNA > 10000 | seu$nCount_RNA < 1000 | seu$nCount_RNA > 50000 | seu$percent.mt > 40, "low", "high")
seu <- subset(seu, QC == "high")

# --- 细胞周期评分 ---
seu <- Seurat::JoinLayers(seu)
seu <- Seurat::NormalizeData(seu)
seu <- Seurat::CellCycleScoring(seu,
                        s.features   = Seurat::cc.genes.updated.2019$s.genes,
                        g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)

# --- 双细胞检测 ---
seu <- scMMR::ComputeDoublets(seu, split.by = "orig.ident", num.cores = 10)

# Fig2B: Doublet分类UMAP
pB <- scMMR::DimPlot2(seu, reduction = "umap", group.by = "DF.classifications_2",
                        palcolor = c("grey", pal_lan),
                        theme_use = "theme_blank",
                        cells.highlight = TRUE, sizes.highlight = 2,
                        label = FALSE,
                        raster = TRUE, raster.dpi = c(512, 512))

seu <- subset(seu, DF.classifications == "Singlet")

# --- 游离RNA污染估计 (DecontX) ---
QuickCluster <- function(object) {
  object <- Seurat::NormalizeData(object)
  object <- Seurat::FindVariableFeatures(object, nfeatures = 2000)
  object <- Seurat::ScaleData(object)
  object <- Seurat::RunPCA(object)
  object <- Seurat::FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- Seurat::FindClusters(object, resolution = 0.1)
  return(object)
}

seu.list <- Seurat::SplitObject(seu, split.by = "orig.ident")
seu.list <- lapply(seu.list, QuickCluster)
clusters <- lapply(seu.list, function(x) x$seurat_clusters) %>% Reduce(c, .)
seu$quick_clusters <- clusters[rownames(seu@meta.data)]
rm(seu.list); gc()

seu <- Seurat::JoinLayers(seu)
seu <- scMMR::ComputeAmbientRNA(seu, split.by = "orig.ident", cluster.name = "quick_clusters")

# Fig2C: DecontX污染程度UMAP热图
pC <- scMMR::FeaturePlot2(seu, reduction = "umap",
                           features = "decontX_contamination",
                           theme_use = "theme_blank",
                           raster = TRUE, raster.dpi = c(512, 512))

# --- 组合Figure 2并保存 ---
pBC <- pB | pC
p_fig1 <- pA / pBC + patchwork::plot_layout(heights = c(3, 1)) +
  patchwork::plot_annotation(tag_levels = "A")
ggplot2::ggsave("./figures/Fig2_QC.pdf", p_fig1, width = 10, height = 10)

# 过滤高污染细胞
seu <- subset(seu, subset = decontX_contamination < 0.5)
seu[["decontX"]] <- NULL

# --- 保存质控后数据 ---
qs::qsave(seu, "./data/seu_QC.qs")

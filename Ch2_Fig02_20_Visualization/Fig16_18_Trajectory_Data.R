# ============================================================
# Fig16_18_Trajectory_Data.R
# Parathyroid cells 轨迹分析数据生成
# Harmony → CytoTRACE2 → Monocle3 (3 lineages) → Velocity
# 输出: ./out/trajectory/seu_with_pseudotime.qs
# 注意: Monocle3 lineage 选择需要交互式操作 (CellSelector)
# ============================================================

library(scMMR)
library(Seurat)
library(qs)
library(monocle3)

# --- 加载并提取 Parathyroid cells ---
seu <- qs::qread("./data/seu.qs")
seu <- subset(seu, subset = celltype %in% "Parathyroid cells")

# --- 标准流程: Normalize → PCA → Harmony → UMAP → Cluster ---
seu <- Seurat::NormalizeData(seu) %>%
  Seurat::FindVariableFeatures(nfeatures = 2000) %>%
  Seurat::ScaleData() %>%
  Seurat::RunPCA() %>%
  harmony::RunHarmony(group.by.vars = "sample",
                      reduction.use = "pca", dims.use = 1:50) %>%
  Seurat::RunUMAP(reduction = "harmony", dims = 1:30) %>%
  Seurat::FindNeighbors(reduction = "harmony", dims = 1:30, k.param = 20) %>%
  Seurat::FindClusters(resolution = 0.1)
# --- CytoTRACE2 ---
seu <- scMMR::RunCytoTRACE2(seu)

# --- Monocle3 轨迹分析 (交互式选择 root cells) ---
p <- Seurat::DimPlot(seu, reduction = "umap", group.by = "celltype")
selected_cells <- Seurat::CellSelector(plot = p)
root_barcodes <- selected_cells

seu <- scMMR::RunMonocle3(
  seu,
  reduction     = "umap",
  root_cells    = root_barcodes,
  use_partition = FALSE
)

# --- 提取 lineages (交互式选择节点) ---
cds <- seu@tools$Monocle3$cds

cds_sub1 <- monocle3::choose_graph_segments(
  cds, starting_pr_node = "Y_251", ending_pr_nodes = "Y_12")
seu$Lineage1 <- NA
seu$Lineage1[colnames(cds_sub1)] <- seu$Monocle3_Pseudotime[colnames(cds_sub1)]

cds_sub2 <- monocle3::choose_graph_segments(
  cds, starting_pr_node = "Y_251", ending_pr_nodes = "Y_341")
seu$Lineage2 <- NA
seu$Lineage2[colnames(cds_sub2)] <- seu$Monocle3_Pseudotime[colnames(cds_sub2)]

cds_sub3 <- monocle3::choose_graph_segments(
  cds, starting_pr_node = "Y_251", ending_pr_nodes = "Y_501")
seu$Lineage3 <- NA
seu$Lineage3[colnames(cds_sub3)] <- seu$Monocle3_Pseudotime[colnames(cds_sub3)]

# --- 加入 RNA velocity metadata ---
meta_velo <- read.csv("./out/trajectory/velocity_metadata.csv", row.names = 1)
cols_to_add <- c("velocity_pseudotime", "velocity_confidence", "velocity_length")
seu <- Seurat::AddMetaData(seu, metadata = meta_velo[Seurat::Cells(seu), cols_to_add, drop = FALSE])

# --- 保存 ---
qs::qsave(seu, "./out/trajectory/seu_with_pseudotime.qs")
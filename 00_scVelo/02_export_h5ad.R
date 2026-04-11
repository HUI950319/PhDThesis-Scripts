# ============================================================
# Fig16_18_Velocity_Data.R
# RNA Velocity 数据准备: Seurat → add loom → export h5ad
# 输出: ./out/velocyto/para_velocity.h5ad
# 注意: 需要先运行 velocyto CLI 生成 loom 文件并合并为 qs
# ============================================================

library(Seurat)
library(qs)

# --- 加载含拟时序信息的 Seurat 对象 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime_clean.qs")

# --- 添加 spliced/unspliced assay ---
source("./scripts/scVelo/01_velocity_prep.R")
seu <- add_loom_to_seurat(
  seu          = seu,
  loom_qs_file = "./data/loom_merged.qs"
)
# --- 提取 spliced/unspliced 矩阵 ---
emat <- Seurat::GetAssayData(seu, assay = "spliced", layer = "counts")
nmat <- Seurat::GetAssayData(seu, assay = "unspliced", layer = "counts")

# --- 使用 anndataR 导出 h5ad ---
outdir <- "./out/velocyto"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

obs_df <- seu@meta.data

obsm <- list()
if ("umap" %in% Seurat::Reductions(seu))
  obsm[["X_umap"]] <- as.matrix(Seurat::Embeddings(seu, "umap"))
if ("pca" %in% Seurat::Reductions(seu))
  obsm[["X_pca"]] <- as.matrix(Seurat::Embeddings(seu, "pca"))
if ("harmony" %in% Seurat::Reductions(seu))
  obsm[["X_harmony"]] <- as.matrix(Seurat::Embeddings(seu, "harmony"))
layers <- list(
  spliced   = Matrix::t(emat),
  unspliced = Matrix::t(nmat)
)

adata <- anndataR::AnnData(
  X      = Matrix::t(emat),
  obs    = obs_df,
  var    = data.frame(row.names = rownames(emat)),
  obsm   = obsm,
  layers = layers
)

h5ad_file <- file.path(outdir, "para_velocity.h5ad")
anndataR::write_h5ad(adata, h5ad_file, mode = "w")
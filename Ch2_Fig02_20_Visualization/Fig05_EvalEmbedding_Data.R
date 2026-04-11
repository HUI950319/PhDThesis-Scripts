# ============================================================
# Fig05_EvalEmbedding_Data.R
# 生成 Fig05_EvalEmbedding.R 所需的画图数据
# 输入: seu.qs + scMMR model.pt
# 输出: out/evaluate_embedding/
#   - silhouette_per_type.csv  (Panel A)
#   - eval_results.qs          (Panel B: dist_cor)
#   - umap_comparison_data.qs  (Panel C: UMAP坐标)
# ============================================================

library(scMMR)
library(Seurat)
library(qs)

use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# --- 参数配置 ---
seu_path      <- "./data/seu.qs"
model_path    <- system.file("extdata", "model.pt", package = "scMMR")
out_dir       <- "./out/evaluate_embedding"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

celltype_col <- "celltype"
k_values     <- c(10L, 20L, 50L, 100L, 200L)
n_pairs      <- 10000L
max_cells    <- 20000L
seed         <- 42L

# --- 1. 加载数据 ---
seu <- qs::qread(seu_path)

if (is.null(Seurat::Reductions(seu, "pca"))) {
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = 50, verbose = FALSE)
}

# --- 2. DNN预测 (提取512维embedding) ---
pred <- DNN_predict(
  query            = seu,
  model_path       = model_path,
  true_label_col   = celltype_col,
  return_embedding = TRUE,
  device           = "auto"
)

# --- 3. Embedding质量评估 (silhouette + dist_cor) ---
eval_res <- EvaluateEmbedding(
  embedding     = pred$shared_embedding,
  seurat_obj    = seu,
  cell_type_col = celltype_col,
  n_pcs         = 50L,
  k_values      = k_values,
  n_pairs       = n_pairs,
  max_cells     = max_cells,
  seed          = seed,
  verbose       = TRUE
)

# --- 4. UMAP对比数据 (PCA vs DNN) ---
seu <- Seurat::RunUMAP(seu, reduction = "pca", dims = 1:30,
                       reduction.name = "umap_pca", verbose = FALSE)

seu[["dnn_emb"]] <- Seurat::CreateDimReducObject(
  embeddings = pred$shared_embedding[colnames(seu), ],
  key = "DNN_", assay = Seurat::DefaultAssay(seu)
)
seu <- Seurat::RunUMAP(seu, reduction = "dnn_emb",
                       dims = 1:ncol(pred$shared_embedding),
                       reduction.name = "umap_dnn", verbose = FALSE)

df_pca <- data.frame(
  Seurat::Embeddings(seu, "umap_pca"),
  celltype = seu@meta.data[[celltype_col]]
)
colnames(df_pca)[1:2] <- c("UMAP_1", "UMAP_2")

df_dnn_emb <- data.frame(
  Seurat::Embeddings(seu, "umap_dnn"),
  celltype = seu@meta.data[[celltype_col]]
)
colnames(df_dnn_emb)[1:2] <- c("UMAP_1", "UMAP_2")

# --- 5. 保存 ---
write.csv(eval_res$consistency$silhouette$per_type,
          file.path(out_dir, "silhouette_per_type.csv"), row.names = FALSE)

qs::qsave(list(eval_res = eval_res),
           file.path(out_dir, "eval_results.qs"))

qs::qsave(list(df_pca = df_pca, df_dnn_emb = df_dnn_emb),
           file.path(out_dir, "umap_comparison_data.qs"))

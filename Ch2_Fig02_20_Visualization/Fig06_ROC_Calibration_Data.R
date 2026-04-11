# ============================================================
# Fig06_ROC_Calibration_Data.R
# 生成 Figure 6 所需的预测概率矩阵
# 输入: seu.qs, seu_GSE190773.qs, seu_GSE233962.qs + model.pt
# 输出: out/Figure5_intermediate.rds
#   - true_ref / probs_ref      (Reference)
#   - true_v1  / probs_v1       (Validation 1: GSE190773)
#   - true_v2  / probs_v2       (Validation 2: GSE233962)
# ============================================================

library(scMMR)
library(Seurat)
library(reticulate)

use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# --- 路径 ---
model_path <- system.file("extdata", "model.pt", package = "scMMR")
out_dir <- "./out"

# --- 辅助函数: 获取预测标签 + 全类别概率矩阵 ---
get_pred_with_probs <- function(seu, model_path, true_label_col = NULL) {
  pred <- DNN_predict(
    query          = seu,
    model_path     = model_path,
    true_label_col = true_label_col,
    device         = "auto"
  )

  # 通过Python获取完整softmax概率矩阵
  py_env <- scMMR:::.scMMR_env
  model  <- py_env$mt_load_model(normalizePath(model_path))

  counts_mat <- Seurat::GetAssayData(seu, layer = "counts")
  obs_df     <- reticulate::r_to_py(seu@meta.data)
  var_names  <- rownames(seu)
  obs_names  <- colnames(seu)

  obsm_dict <- list()
  if ("umap" %in% Seurat::Reductions(seu)) {
    obsm_dict[["X_umap"]] <- Seurat::Embeddings(seu, "umap")
  }

  adata <- py_env$mt_srt_to_adata(
    counts_matrix = counts_mat,
    obs_df        = obs_df,
    var_names     = var_names,
    obs_names     = obs_names,
    obsm_dict     = if (length(obsm_dict) > 0) obsm_dict else NULL
  )

  X <- py_env$mt_prepare_query(adata, model)
  py_result <- model$predict(X)
  probs     <- reticulate::py_to_r(py_result[["probabilities"]])

  class_names <- reticulate::py_to_r(model$label_encoder$classes_)
  colnames(probs) <- class_names
  rownames(probs) <- colnames(seu)

  list(predictions = pred$predictions,
       probs       = probs,
       class_names = class_names)
}

# --- 1. Reference数据集 (已知celltype) ---
seu_ref  <- qs::qread("./data/seu.qs")
pred_ref <- get_pred_with_probs(seu_ref, model_path, true_label_col = "celltype")

true_ref  <- factor(seu_ref$celltype)
probs_ref <- pred_ref$probs
rm(seu_ref); gc()

# --- 2. GSE190773 (外部验证, 已有celltype注释) ---
seu_v1  <- qs::qread("./data/seu_GSE190773.qs")
pred_v1 <- get_pred_with_probs(seu_v1, model_path, true_label_col = "celltype")

true_v1  <- factor(seu_v1$celltype)
probs_v1 <- pred_v1$probs
rm(seu_v1); gc()

# --- 3. GSE233962 (外部验证, 已有celltype注释) ---
seu_v2  <- qs::qread("./data/seu_GSE233962.qs")
pred_v2 <- get_pred_with_probs(seu_v2, model_path, true_label_col = "celltype")

true_v2  <- factor(seu_v2$celltype)
probs_v2 <- pred_v2$probs
rm(seu_v2); gc()

# --- 保存 ---
saveRDS(list(
  true_ref = true_ref, probs_ref = probs_ref,
  true_v1  = true_v1,  probs_v1  = probs_v1,
  true_v2  = true_v2,  probs_v2  = probs_v2
), file.path(out_dir, "Figure5_intermediate.rds"))

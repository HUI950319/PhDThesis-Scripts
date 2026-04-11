# ============================================================
# Fig07_08_Data.R
# 生成 Fig07 / Fig08 所需的预测数据
# DNN_predict → 添加预测标签与UMAP坐标 → 保存 predicted.qs
# ============================================================

library(scMMR)
library(Seurat)
library(qs)

use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# --- 路径 ---
model_path <- system.file("extdata", "model.pt", package = "scMMR")
out_dir    <- "./out/06_GSE_prediction"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
datasets <- list(
  GSE190773 = "./data/seu_GSE190773.qs",
  GSE233962 = "./data/seu_GSE233962.qs"
)

# --- 预测并保存 ---
for (ds_name in names(datasets)) {
  seu <- qs::qread(datasets[[ds_name]])

  pred <- scMMR::DNN_predict(
    query      = seu,
    model_path = model_path,
    explain    = FALSE,
    device     = "gpu"
  )

  seu <- Seurat::AddMetaData(seu, pred$predictions)
  qs::qsave(seu, file.path(out_dir, paste0(ds_name, "_predicted.qs")))
}
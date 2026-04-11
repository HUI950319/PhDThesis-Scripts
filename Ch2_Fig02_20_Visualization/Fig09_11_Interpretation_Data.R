# ============================================================
# Fig09_11_Interpretation_Data.R
# 模型可解释性分析：Gene / KEGG / Regulon Importance (IG)
# DNN_predict(explain=TRUE) × 3 GMT → fig9_intermediate.qs
# ============================================================

library(scMMR)
library(Seurat)
library(qs)
library(dplyr)

use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# --- 路径 ---
model_path   <- system.file("extdata", "model.pt", package = "scMMR")
hallmark_gmt <- system.file("extdata/gmt", "h.all.v2022.1.Hs.symbols.gmt",
                            package = "scMMR")
kegg_gmt     <- system.file("extdata/gmt", "c2.cp.kegg.v2022.1.Hs.symbols.gmt",
                            package = "scMMR")
regulon_gmt  <- "./data/regulons.gmt"
out_dir      <- "./out/09_model_interpretation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

celltype_col    <- "celltype"
n_cells_explain <- 200L
top_k_global    <- 50L
top_k_class     <- 50L

# --- 加载数据 ---
seu <- qs::qread("./data/seu.qs")

# ============================================================
# DNN_predict × 3 GMT
# ============================================================

# --- Hallmark ---
pred_hallmark <- scMMR::DNN_predict(
  query            = seu,
  model_path       = model_path,
  true_label_col   = celltype_col,
  explain          = TRUE,  top_k_global     = top_k_global,
  top_k_class      = top_k_class,
  n_cells_explain  = n_cells_explain,
  pathway_gmt      = hallmark_gmt,
  return_embedding = FALSE,
  device           = "auto"
)

# --- KEGG ---
pred_kegg <- scMMR::DNN_predict(
  query            = seu,
  model_path       = model_path,
  true_label_col   = celltype_col,
  explain          = TRUE,
  top_k_global     = top_k_global,
  top_k_class      = top_k_class,
  n_cells_explain  = n_cells_explain,
  pathway_gmt      = kegg_gmt,
  return_embedding = FALSE,
  device           = "auto"
)

# --- Regulon (SCENIC) ---
pred_regulon <- scMMR::DNN_predict(
  query            = seu,
  model_path       = model_path,
  true_label_col   = celltype_col,  explain          = TRUE,
  top_k_global     = top_k_global,
  top_k_class      = top_k_class,
  n_cells_explain  = n_cells_explain,
  pathway_gmt      = regulon_gmt,
  return_embedding = FALSE,
  device           = "auto"
)

# ============================================================
# 整理 per-class importance tables
# ============================================================

gene_pc  <- pred_hallmark$imp_per_class
ct_levels <- sort(unique(seu@meta.data[[celltype_col]]))

# --- Hallmark per-class ---
pw_hallmark <- pred_hallmark$pathway_scores
pw_ct_h     <- seu@meta.data[rownames(pw_hallmark), celltype_col]

hallmark_long <- do.call(rbind, lapply(ct_levels, function(ct) {
  idx <- which(pw_ct_h == ct)
  scores <- if (length(idx) > 1) {
    colMeans(as.matrix(pw_hallmark[idx, , drop = FALSE]))
  } else {
    as.numeric(pw_hallmark[idx, ])
  }
  data.frame(cell_type = ct, name = colnames(pw_hallmark),
             importance = scores, stringsAsFactors = FALSE)
}))
hallmark_long$name_clean <- gsub("^HALLMARK_", "", hallmark_long$name)
hallmark_long$name_clean <- gsub("_", " ", hallmark_long$name_clean)
hallmark_long$name_clean <- paste0(
  toupper(substr(hallmark_long$name_clean, 1, 1)),
  tolower(substr(hallmark_long$name_clean, 2, nchar(hallmark_long$name_clean)))
)
hallmark_top <- hallmark_long %>%
  dplyr::group_by(cell_type) %>%
  dplyr::slice_max(importance, n = top_k_class) %>%
  dplyr::ungroup()

# --- KEGG per-class ---
pw_kegg <- pred_kegg$pathway_scores
pw_ct_k <- seu@meta.data[rownames(pw_kegg), celltype_col]

kegg_long <- do.call(rbind, lapply(ct_levels, function(ct) {
  idx <- which(pw_ct_k == ct)
  scores <- if (length(idx) > 1) {
    colMeans(as.matrix(pw_kegg[idx, , drop = FALSE]))
  } else {
    as.numeric(pw_kegg[idx, ])
  }
  data.frame(cell_type = ct, name = colnames(pw_kegg),             importance = scores, stringsAsFactors = FALSE)
}))
kegg_long$name_clean <- gsub("^KEGG_", "", kegg_long$name)
kegg_long$name_clean <- gsub("_", " ", kegg_long$name_clean)
kegg_long$name_clean <- paste0(
  toupper(substr(kegg_long$name_clean, 1, 1)),
  tolower(substr(kegg_long$name_clean, 2, nchar(kegg_long$name_clean)))
)
kegg_top <- kegg_long %>%
  dplyr::group_by(cell_type) %>%
  dplyr::slice_max(importance, n = top_k_class) %>%
  dplyr::ungroup()

# --- Regulon per-class ---
pw_regulon <- pred_regulon$pathway_scores
pw_ct_r    <- seu@meta.data[rownames(pw_regulon), celltype_col]

regulon_long <- do.call(rbind, lapply(ct_levels, function(ct) {
  idx <- which(pw_ct_r == ct)
  scores <- if (length(idx) > 1) {
    colMeans(as.matrix(pw_regulon[idx, , drop = FALSE]))
  } else {
    as.numeric(pw_regulon[idx, ])
  }
  data.frame(cell_type = ct, name = colnames(pw_regulon),
             importance = scores, stringsAsFactors = FALSE)}))
regulon_long$name_clean <- gsub("\\(\\d+g\\)$", "", regulon_long$name)
regulon_top <- regulon_long %>%
  dplyr::group_by(cell_type) %>%
  dplyr::slice_max(importance, n = top_k_class) %>%
  dplyr::ungroup()

# --- 保存中间结果 ---
qs::qsave(list(
  gene_pc      = gene_pc,
  hallmark_top = hallmark_top,
  kegg_top     = kegg_top,
  regulon_top  = regulon_top,
  params = list(
    n_cells_explain = n_cells_explain,
    top_k_global    = top_k_global,
    top_k_class     = top_k_class
  )
), file.path(out_dir, "fig9_intermediate.qs"))
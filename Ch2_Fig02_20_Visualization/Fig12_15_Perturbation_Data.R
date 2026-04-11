# ============================================================
# Fig12_15_Perturbation_Data.R
# 细胞比例与表达扰动分析 (使用全部数据，无subsample)
# 4 methods: miloR, scMMR RankPercent (edgeR GLM),
#            Augur, scMMR RankPerturbation
# 输出: ./out/perturbation/
# ============================================================

library(scMMR)
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(Augur)
library(qs)
library(ggplot2)
library(dplyr)
library(S4Vectors)
library(FNN)

use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# --- 路径与参数 ---
model_path   <- system.file("extdata", "model.pt", package = "scMMR")
out_root     <- "./out/perturbation"
celltype_col <- "celltype"
group_col    <- "group"
sample_col   <- "sample"
seed         <- 42L

k_graph    <- 30L
d_pca      <- 30L
prop_nhood <- 0.2

k_percent      <- 30L
prop_percent   <- 0.2
min_cells_pct  <- 10L

n_threads_augur <- 8L
augur_subsample <- 20L
augur_folds     <- 3L

n_permutations <- 1000L
min_cells_pert <- 5L
n_pcs_pert     <- 20L

comparisons <- list(
  "PHPT_vs_Normal" = c("Normal", "PHPT"),
  "SHPT_vs_Normal" = c("Normal", "SHPT")
)

milo_dir <- file.path(out_root, "miloR")
pct_dir  <- file.path(out_root, "scMMR_percent")
augur_dir <- file.path(out_root, "augur")
pert_dir  <- file.path(out_root, "scMMR_perturbation")
for (d in c(out_root, milo_dir, pct_dir, augur_dir, pert_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. 加载全部数据
# ============================================================

seu <- qs::qread("./data/seu.qs")

if (is.null(Reductions(seu, "pca"))) {
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = 50, verbose = FALSE)
}
if (is.null(Reductions(seu, "umap"))) {
  seu <- Seurat::RunUMAP(seu, dims = 1:30, verbose = FALSE)
}

qs::qsave(seu@meta.data, file.path(out_root, "seu_meta.data.qs"))

# ============================================================
# 2. DNN embedding
# ============================================================

pred <- scMMR::DNN_predict(
  query            = seu,
  model_path       = model_path,
  true_label_col   = celltype_col,
  return_embedding = TRUE,
  device           = "auto"
)
emb <- pred$shared_embedding
seu <- Seurat::AddMetaData(seu, pred$predictions)

# ============================================================
# 3. miloR
# ============================================================

milo_results   <- list()
milo_summaries <- list()

for (comp_name in names(comparisons)) {
  pair <- comparisons[[comp_name]]
  keep <- seu@meta.data[[group_col]] %in% pair
  seu_pair <- seu[, keep]

  meta_pair <- seu_pair@meta.data[, c(sample_col, group_col, celltype_col)]
  colnames(meta_pair) <- c("sample", "group", "celltype")

  sce <- SingleCellExperiment(
    assays      = list(counts = Seurat::GetAssayData(seu_pair, layer = "counts")),
    reducedDims = SimpleList(PCA  = Seurat::Embeddings(seu_pair, "pca"),
                             UMAP = Seurat::Embeddings(seu_pair, "umap"))
  )
  milo <- Milo(sce)
  colData(milo) <- DataFrame(meta_pair)

  k_use <- ifelse(ncol(seu_pair) < 500, 15L, k_graph)
  milo <- buildGraph(milo, k = k_use, d = d_pca, reduced.dim = "PCA")
  milo <- makeNhoods(milo, prop = prop_nhood, k = k_use, d = d_pca,
                     refined = TRUE, reduced_dims = "PCA")
  milo <- countCells(milo, meta.data = as.data.frame(colData(milo)),
                     sample = "sample")

  milo_design <- as.data.frame(colData(milo)) %>%
    dplyr::select(sample, group) %>% dplyr::distinct()
  milo_design$group <- factor(milo_design$group, levels = pair)
  rownames(milo_design) <- milo_design$sample

  milo <- calcNhoodDistance(milo, d = d_pca, reduced.dim = "PCA")
  da_res <- testNhoods(milo, design = ~ group, design.df = milo_design,
                       reduced.dim = "PCA")
  da_res <- annotateNhoods(milo, da_res, coldata_col = "celltype")
  da_res$celltype <- ifelse(da_res$celltype_fraction < 0.7, "Mixed", da_res$celltype)
  da_res$comparison <- comp_name

  da_summary <- da_res %>%
    dplyr::filter(celltype != "Mixed") %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(
      mean_logFC = mean(logFC, na.rm = TRUE),
      sd_logFC   = sd(logFC, na.rm = TRUE),
      n_nhoods   = n(),
      n_sig      = sum(SpatialFDR < 0.1, na.rm = TRUE),
      frac_sig   = mean(SpatialFDR < 0.1, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::arrange(desc(abs(mean_logFC))) %>%
    dplyr::mutate(comparison = comp_name)
  milo_results[[comp_name]]   <- da_res
  milo_summaries[[comp_name]] <- da_summary

  write.csv(da_res, file.path(milo_dir,
    sprintf("miloR_%s_da_results.csv", comp_name)), row.names = FALSE)
  write.csv(da_summary, file.path(milo_dir,
    sprintf("miloR_%s_summary.csv", comp_name)), row.names = FALSE)
}

write.csv(dplyr::bind_rows(milo_summaries),
          file.path(milo_dir, "miloR_combined_summary.csv"), row.names = FALSE)

# ============================================================
# 4. scMMR RankPercent (edgeR GLM mode, with sample_col)
# ============================================================

pct_results <- list()

for (comp_name in names(comparisons)) {
  pair <- comparisons[[comp_name]]
  keep <- seu@meta.data[[group_col]] %in% pair
  emb_pair  <- emb[keep, ]
  meta_pair <- seu@meta.data[keep, ]

  res <- scMMR::RankPercent(
    embedding     = emb_pair,
    cell_meta     = meta_pair,
    cell_type_col = celltype_col,    condition_col = group_col,
    sample_col    = sample_col,
    conditions    = pair,
    k             = k_percent,
    prop_sample   = prop_percent,
    test          = "fisher",
    min_cells     = min_cells_pct,
    seed          = seed,
    verbose       = TRUE
  )

  res$cell_type_summary$comparison <- comp_name
  pct_results[[comp_name]] <- res

  write.csv(res$da_results, file.path(pct_dir,
    sprintf("percent_%s_da_results.csv", comp_name)), row.names = FALSE)
  write.csv(res$cell_type_summary, file.path(pct_dir,
    sprintf("percent_%s_summary.csv", comp_name)), row.names = FALSE)
}

write.csv(dplyr::bind_rows(lapply(pct_results, function(x) x$cell_type_summary)),
          file.path(pct_dir, "percent_combined.csv"), row.names = FALSE)

# ============================================================
# 5. Augur
# ============================================================

all_aucs <- list()

for (comp_name in names(comparisons)) {  pair <- comparisons[[comp_name]]
  keep <- seu@meta.data[[group_col]] %in% pair
  seu_pair <- seu[, keep]
  seu_pair@meta.data[["condition"]] <- factor(
    seu_pair@meta.data[[group_col]], levels = pair)

  augur_res <- Augur::calculate_auc(
    seu_pair,
    cell_type_col  = celltype_col,
    label_col      = "condition",
    n_threads      = n_threads_augur,
    folds          = augur_folds,
    subsample_size = augur_subsample
  )

  augur_auc <- augur_res$AUC
  augur_auc$comparison <- comp_name
  all_aucs[[comp_name]] <- augur_auc

  write.csv(augur_auc, file.path(augur_dir,
    sprintf("augur_%s_auc.csv", comp_name)), row.names = FALSE)
}

write.csv(dplyr::bind_rows(all_aucs),
          file.path(augur_dir, "augur_combined_auc.csv"), row.names = FALSE)

# ============================================================
# 6. scMMR RankPerturbation
# ============================================================
pert_results <- list()

for (comp_name in names(comparisons)) {
  pair <- comparisons[[comp_name]]
  keep <- seu@meta.data[[group_col]] %in% pair
  emb_pair  <- emb[keep, ]
  meta_pair <- seu@meta.data[keep, ]

  res <- scMMR::RankPerturbation(
    embedding      = emb_pair,
    cell_meta      = meta_pair,
    cell_type_col  = celltype_col,
    condition_col  = group_col,
    conditions     = pair,
    method         = "wasserstein",
    n_permutations = n_permutations,
    min_cells      = min_cells_pert,
    n_pcs          = n_pcs_pert,
    seed           = seed,
    verbose        = TRUE
  )

  res$results$comparison <- comp_name
  pert_results[[comp_name]] <- res

  write.csv(res$results, file.path(pert_dir,
    sprintf("perturbation_%s_results.csv", comp_name)), row.names = FALSE)
}

write.csv(dplyr::bind_rows(lapply(pert_results, function(x) x$results)),
          file.path(pert_dir, "perturbation_combined.csv"), row.names = FALSE)
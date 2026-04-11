# ============================================================
# Fig04_Benchmark.R
# scMMR DNN vs 11 methods: 10次重复10折交叉验证
# 数据集: 4个胰腺公开数据 + 3个甲状旁腺子集
# Figure 4A: ComplexHeatmap / Figure 4B: funkyheatmap
# ============================================================

library(scMMR)
library(Seurat)
library(tidyverse)
library(mlr3verse)
library(data.table)
library(scRNAseq)
library(SingleCellExperiment)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

out_dir <- "./data/03_benchmark"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 参数配置 ---
N_REPEATS  <- 10L
N_FOLDS    <- 10L
N_EPOCHS   <- 50L

method_names <- c("DNN", "RF", "XGBoost", "SVM", "ElasticNet",
                  "LDA", "NaiveBayes", "KNN", "NNet",
                  "SingleR", "Seurat_LT", "CellTypist")

metric_names <- c("Accuracy", "Balanced_Acc", "Macro_F1",
                  "Weighted_F1", "Kappa", "Min_Recall")

type_map <- c(DNN = "Deep Learning",
              RF = "Machine Learning", XGBoost = "Machine Learning",
              SVM = "Machine Learning", ElasticNet = "Machine Learning",
              LDA = "Machine Learning", NaiveBayes = "Machine Learning",
              KNN = "Machine Learning", NNet = "Machine Learning",
              SingleR = "Single-cell", Seurat_LT = "Single-cell",
              CellTypist = "Single-cell")
type_colors <- c("Deep Learning" = "#00468B",
                 "Machine Learning" = "#42B540",
                 "Single-cell" = "#0099B4")

# --- 辅助函数: 计算6个分类指标 ---
compute_metrics <- function(true_vec, pred_vec) {
  na_result <- setNames(rep(NA_real_, length(metric_names)), metric_names)
  valid <- !is.na(true_vec) & !is.na(pred_vec) & true_vec != "" & pred_vec != ""
  if (sum(valid) < 2) return(na_result)
  tru <- true_vec[valid]; prd <- pred_vec[valid]
  lvls <- sort(unique(tru)); n <- length(tru)

  pc <- vapply(lvls, function(ct) {
    tp <- sum(tru == ct & prd == ct)
    fp <- sum(tru != ct & prd == ct)
    fn <- sum(tru == ct & prd != ct)
    pr <- if (tp + fp > 0) tp / (tp + fp) else 0
    rc <- if (tp + fn > 0) tp / (tp + fn) else 0
    f1 <- if (pr + rc > 0) 2 * pr * rc / (pr + rc) else 0
    c(pr, rc, f1, sum(tru == ct))
  }, numeric(4))
  rownames(pc) <- c("prec", "recall", "f1", "n")

  acc <- mean(tru == prd)
  bacc <- mean(pc["recall", ])
  macro_f1 <- mean(pc["f1", ])
  weighted_f1 <- sum(pc["f1", ] * pc["n", ] / n)
  min_recall <- min(pc["recall", ])
  ct <- table(factor(tru, lvls), factor(prd, lvls))
  po <- sum(diag(ct)) / n; pe <- sum(rowSums(ct) * colSums(ct)) / n^2
  kappa <- if (pe < 1) (po - pe) / (1 - pe) else 1

  c(Accuracy = acc, Balanced_Acc = bacc, Macro_F1 = macro_f1,
    Weighted_F1 = weighted_f1, Kappa = kappa, Min_Recall = min_recall)
}

# --- 辅助函数: 从mlr3 benchmark提取预测 ---
extract_bmr_preds <- function(bmr, cell_names) {
  scores_dt <- as.data.table(bmr$score(msr("classif.acc")))
  n <- length(cell_names)
  out <- list()
  for (lid in unique(scores_dt$learner_id)) {
    sub <- scores_dt[learner_id == lid]
    p <- character(n)
    for (j in seq_len(nrow(sub))) {
      idx <- sub$resampling[[j]]$test_set(sub$iteration[j])
      p[idx] <- as.character(sub$prediction[[j]]$response)
    }
    out[[lid]] <- setNames(p, cell_names)
  }
  out
}

# --- 辅助函数: SCE转Seurat (胰腺数据集) ---
sce_to_seurat <- function(sce, dataset_name) {
  cts <- if ("counts" %in% assayNames(sce)) counts(sce) else assay(sce, 1)
  if (dataset_name %in% c("Muraro", "Xin")) {
    syms <- gsub("_", "-", as.character(rowData(sce)$symbol))
    ok <- !is.na(syms) & syms != "" & !duplicated(syms)
    cts <- cts[ok, ]; rownames(cts) <- syms[ok]
  }

  label_map <- list(
    Baron = list(col = "label", map = c(
      acinar = "acinar", alpha = "alpha", beta = "beta", delta = "delta",
      ductal = "ductal", endothelial = "endothelial", epsilon = "epsilon",
      gamma = "gamma", activated_stellate = "stellate",
      quiescent_stellate = "stellate", macrophage = "macrophage",
      mast = "mast", schwann = "schwann", t_cell = "t_cell")),
    Muraro = list(col = "label", map = c(
      acinar = "acinar", alpha = "alpha", beta = "beta", delta = "delta",
      duct = "ductal", endothelial = "endothelial", epsilon = "epsilon",
      pp = "gamma", mesenchymal = "mesenchyme")),
    Segerstolpe = list(col = "cell type", map = c(
      "acinar cell" = "acinar", "alpha cell" = "alpha", "beta cell" = "beta",
      "delta cell" = "delta", "ductal cell" = "ductal",
      "endothelial cell" = "endothelial", "epsilon cell" = "epsilon",
      "gamma cell" = "gamma", "PSC cell" = "stellate", "mast cell" = "mast")),
    Xin = list(col = "cell.type", map = c(
      alpha = "alpha", beta = "beta", delta = "delta", PP = "gamma")))

  info <- label_map[[dataset_name]]
  mapped <- info$map[as.character(colData(sce)[[info$col]])]
  keep <- !is.na(mapped)
  cts <- cts[, keep]; mapped <- mapped[keep]
  rn <- gsub("_", "-", rownames(cts))
  dup <- duplicated(rn)
  if (any(dup)) { cts <- cts[!dup, ]; rn <- rn[!dup] }
  rownames(cts) <- rn

  seu <- Seurat::CreateSeuratObject(counts = cts, min.cells = 0, min.features = 0)
  seu$cell_type <- as.character(mapped)
  seu
}

# --- 1. 加载数据集 ---

# 甲状旁腺数据 (3个独立数据集)
seu_para <- qs::qread("./data/seu_annotated.qs")
seu_para$cell_type <- seu_para$celltype

seu_GSE190773 <- qs::qread("./data/seu_GSE190773.qs")
seu_GSE190773$cell_type <- seu_GSE190773$celltype

seu_GSE233962 <- qs::qread("./data/seu_GSE233962.qs")
seu_GSE233962$cell_type <- seu_GSE233962$celltype

# 胰腺公开数据 (4个)
datasets <- list(
  Baron       = sce_to_seurat(scRNAseq::BaronPancreasData("human"), "Baron"),
  Muraro      = sce_to_seurat(scRNAseq::MuraroPancreasData(), "Muraro"),
  Segerstolpe = sce_to_seurat(scRNAseq::SegerstolpePancreasData(), "Segerstolpe"),
  Xin         = sce_to_seurat(scRNAseq::XinPancreasData(), "Xin"),
  seu_para      = seu_para,
  seu_GSE190773 = seu_GSE190773,
  seu_GSE233962 = seu_GSE233962)
dataset_names <- names(datasets)

# --- 2. 10次重复10折交叉验证 ---
celltypist_ok <- tryCatch({
  reticulate::import("celltypist"); TRUE
}, error = function(e) FALSE)

# 存储每次repeat的指标: all_metrics[[repeat]][[dataset]][[method]] = named numeric
all_metrics <- list()

for (rep_i in seq_len(N_REPEATS)) {
  set.seed(rep_i * 42)
  rep_metrics <- list()

  for (ds_name in dataset_names) {
    seu <- datasets[[ds_name]]
    true_labels <- setNames(as.character(seu$cell_type), colnames(seu))

    # 去除稀有类型 (需 >= N_FOLDS*2 个细胞)
    ct_tab <- table(true_labels)
    rare <- names(ct_tab[ct_tab < N_FOLDS * 2])
    if (length(rare) > 0) {
      keep <- !true_labels %in% rare
      seu <- seu[, keep]; true_labels <- true_labels[keep]
    }

    # 预处理
    seu <- Seurat::NormalizeData(seu, verbose = FALSE) |>
      Seurat::FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
      Seurat::ScaleData(verbose = FALSE) |>
      Seurat::RunPCA(npcs = 50, verbose = FALSE) |>
      Seurat::RunUMAP(dims = 1:30, verbose = FALSE)

    cell_names  <- colnames(seu)
    n_cells     <- length(cell_names)
    true_labels <- true_labels[cell_names]
    pca_mat     <- Seurat::Embeddings(seu, "pca")
    norm_data   <- Seurat::GetAssayData(seu, layer = "data")
    preds       <- list()

    # 构建分层10折
    task_df <- as.data.frame(pca_mat)
    task_df$cell_type <- factor(true_labels)
    task <- TaskClassif$new(id = "ct", backend = task_df, target = "cell_type")
    task$col_roles$stratum <- "cell_type"
    cv <- rsmp("cv", folds = N_FOLDS)
    cv$instantiate(task)

    # DNN (10折CV)
    tryCatch({
      dnn_p <- setNames(character(n_cells), cell_names)
      for (fold in seq_len(N_FOLDS)) {
        tri <- cv$train_set(fold); ti <- cv$test_set(fold)
        mp <- tempfile(fileext = ".pt")
        suppressMessages({
          DNN_train(input = seu[, tri], label_column = "cell_type",
                    embedding_key = "umap", save_path = mp,
                    n_top_genes = 6000L, num_epochs = N_EPOCHS,
                    batch_size = 256L, device = "auto")
          pred <- DNN_predict(query = seu[, ti], model_path = mp,
                              true_label_col = "cell_type", device = "auto")
        })
        dnn_p[colnames(seu[, ti])] <- as.character(pred$predictions$cell_type_pred)
        file.remove(mp)
      }
      preds[["DNN"]] <- dnn_p
    }, error = function(e) NULL)

    # 8个ML方法 (共享同一套fold)
    tryCatch({
      learners <- list(
        lrn("classif.ranger",      id = "RF",         num.trees = 500),
        lrn("classif.xgboost",     id = "XGBoost",    nrounds = 100, max_depth = 6, verbose = 0),
        lrn("classif.svm",         id = "SVM",        type = "C-classification", kernel = "radial"),
        lrn("classif.glmnet",      id = "ElasticNet",  alpha = 0.5),
        lrn("classif.lda",         id = "LDA"),
        lrn("classif.naive_bayes", id = "NaiveBayes"),
        lrn("classif.kknn",        id = "KNN",        k = 15),
        lrn("classif.nnet",        id = "NNet",       size = 20, MaxNWts = 50000, maxit = 200, trace = FALSE))
      bmr <- suppressWarnings(suppressMessages(
        benchmark(benchmark_grid(task, learners, cv), store_models = FALSE)))
      ml_preds <- extract_bmr_preds(bmr, cell_names)
      for (lid in names(ml_preds)) preds[[lid]] <- ml_preds[[lid]]
    }, error = function(e) NULL)

    # SingleR (10折CV)
    tryCatch({
      sr_p <- setNames(character(n_cells), cell_names)
      for (fold in seq_len(N_FOLDS)) {
        tri <- cv$train_set(fold); ti <- cv$test_set(fold)
        ref <- SingleCellExperiment(assays = list(logcounts = norm_data[, tri]))
        ref$label <- as.character(true_labels[tri])
        qry <- SingleCellExperiment(assays = list(logcounts = norm_data[, ti]))
        sr <- SingleR::SingleR(test = qry, ref = ref, labels = ref$label, de.method = "wilcox")
        sr_p[cell_names[ti]] <- sr$labels
      }
      preds[["SingleR"]] <- sr_p
    }, error = function(e) NULL)

    # Seurat Label Transfer (10折CV)
    tryCatch({
      lt_p <- setNames(character(n_cells), cell_names)
      for (fold in seq_len(N_FOLDS)) {
        tri <- cv$train_set(fold); ti <- cv$test_set(fold)
        anchors <- FindTransferAnchors(reference = seu[, tri], query = seu[, ti],
                                       dims = 1:30, verbose = FALSE)
        transferred <- TransferData(anchorset = anchors,
                                    refdata = as.character(true_labels[tri]),
                                    dims = 1:30, verbose = FALSE)
        lt_p[colnames(seu[, ti])] <- as.character(transferred$predicted.id)
      }
      preds[["Seurat_LT"]] <- lt_p
    }, error = function(e) NULL)

    # CellTypist (10折CV)
    if (celltypist_ok) {
      tryCatch({
        ct_mod <- reticulate::import("celltypist")
        sc_py  <- reticulate::import("scanpy")
        pd     <- reticulate::import("pandas")
        ct_p <- setNames(character(n_cells), cell_names)
        for (fold in seq_len(N_FOLDS)) {
          tri <- cv$train_set(fold); ti <- cv$test_set(fold)
          ref_ad <- sc_py$AnnData(X = as.matrix(t(norm_data[, tri])),
                                  obs = pd$DataFrame(list(cell_type = as.character(true_labels[tri]))))
          ref_ad$var_names <- rownames(norm_data)
          qry_ad <- sc_py$AnnData(X = as.matrix(t(norm_data[, ti])))
          qry_ad$var_names <- rownames(norm_data)
          m <- ct_mod$train(ref_ad, labels = "cell_type", use_SGD = FALSE, n_jobs = 4L, max_iter = 200L)
          pr <- ct_mod$annotate(qry_ad, model = m, majority_voting = FALSE)
          ct_p[cell_names[ti]] <- as.character(reticulate::py_to_r(pr$predicted_labels)[["predicted_labels"]])
        }
        preds[["CellTypist"]] <- ct_p
      }, error = function(e) NULL)
    }

    # 计算本次repeat的指标
    ds_metrics <- list()
    for (method in intersect(method_names, names(preds))) {
      ds_metrics[[method]] <- compute_metrics(true_labels, preds[[method]])
    }
    rep_metrics[[ds_name]] <- ds_metrics
  }
  all_metrics[[rep_i]] <- rep_metrics
}

# --- 3. 汇总结果: 10次repeat取均值±SD ---
# 转为长表: repeat / dataset / method / metric / value
results_long <- do.call(rbind, lapply(seq_len(N_REPEATS), function(r) {
  do.call(rbind, lapply(dataset_names, function(ds) {
    do.call(rbind, lapply(method_names, function(m) {
      vals <- all_metrics[[r]][[ds]][[m]]
      if (is.null(vals)) vals <- setNames(rep(NA_real_, length(metric_names)), metric_names)
      data.frame(Repeat = r, Dataset = ds, Method = m,
                 Metric = names(vals), Value = as.numeric(vals),
                 stringsAsFactors = FALSE)
    }))
  }))
}))

# 均值矩阵 (用于热图)
results_mean <- results_long %>%
  dplyr::group_by(Dataset, Method, Metric) %>%
  dplyr::summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE), .groups = "drop")

# 保存完整结果
write.csv(results_long, file.path(out_dir, "benchmark_10x10cv_full.csv"), row.names = FALSE)
write.csv(results_mean, file.path(out_dir, "benchmark_10x10cv_mean.csv"), row.names = FALSE)
qs::qsave(all_metrics, file.path(out_dir, "benchmark_10x10cv_metrics.qs"))

# --- 4. 排名汇总 ---
mean_acc <- results_mean %>%
  dplyr::filter(Metric == "Accuracy") %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(Avg = mean(Mean, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(Avg))

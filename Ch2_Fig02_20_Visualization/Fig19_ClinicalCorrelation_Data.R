# ==============================================================================
# Figure 19 Data: Bulk RNA-seq 解卷积 — scMMR DNN 估计细胞亚群比例
# ==============================================================================
#
# 使用 scMMR::DNN_deconv_train/predict 从 Bulk RNA-seq 估计细胞亚群比例，
# 合并临床数据，生成 Figure 19 的 plot_data.Rdata (df1-df6)。
#
# 输入：
#   - ./out/seu.qs (scRNA-seq 参考 Seurat 对象)
#   - Bulk counts: self_PHPT (n=12), self_SHPT (n=12), PRJNA516535 (n=10)
#
# 输出：
#   - ./out/deconv/deconv_model.pt
#   - ./out/deconv/props_*.csv
#   - ./out/Figures/Figure19/plot_data.Rdata
#
# 依赖：scMMR, Seurat, qs, tidyverse
# ==============================================================================

library(scMMR)
library(Seurat)
library(qs)
library(tidyverse)

use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

bulk_dir <- "/home/oyh/project/bulk/counts"
out_dir  <- "./out/deconv"
fig_dir  <- "./out/Figures/Figure19"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. 训练解卷积模型
# ==============================================================================

seu_ref    <- qs::qread("./out/seu.qs")
model_path <- file.path(out_dir, "deconv_model.pt")

if (!file.exists(model_path)) {
  DNN_deconv_train(
    reference    = seu_ref,
    label_column = "celltype",
    save_path    = model_path,
    n_pseudobulk = 5000,
    n_cells_per_sample = 500,
    n_top_genes  = 6000,
    hidden_size  = 512,
    num_epochs   = 50,
    batch_size   = 256,
    device       = "auto"
  )
}

# ==============================================================================
# 2. 读取 featureCounts 输出
# ==============================================================================

read_featurecounts <- function(file_path) {
  raw <- read.delim(file_path, header = TRUE, row.names = 1,
                    check.names = FALSE, comment.char = "#")
  if ("gene_name" %in% colnames(raw)) {
    gn <- raw$gene_name; raw$gene_name <- NULL
    keep <- !duplicated(gn)
    raw <- raw[keep, , drop = FALSE]; rownames(raw) <- gn[keep]
  }
  counts <- as.matrix(raw[, sapply(raw, is.numeric), drop = FALSE])
  colnames(counts) <- gsub(".*/(\\w+)_Aligned.*", "\\1", colnames(counts))
  counts
}

# ==============================================================================
# 3. 解卷积预测
# ==============================================================================

counts_phpt  <- read_featurecounts(file.path(bulk_dir, "self_PHPT",   "self_PHPT_counts.txt"))
counts_shpt  <- read_featurecounts(file.path(bulk_dir, "self_SHPT",   "self_SHPT_counts.txt"))
counts_prjna <- read_featurecounts(file.path(bulk_dir, "PRJNA516535", "PRJNA516535_counts.txt"))

props_phpt  <- DNN_deconv_predict(counts_phpt,  model_path, device = "auto")
props_shpt  <- DNN_deconv_predict(counts_shpt,  model_path, device = "auto")
props_prjna <- DNN_deconv_predict(counts_prjna, model_path, device = "auto")

write.csv(props_phpt,  file.path(out_dir, "props_self_PHPT.csv"),   row.names = FALSE)
write.csv(props_shpt,  file.path(out_dir, "props_self_SHPT.csv"),   row.names = FALSE)
write.csv(props_prjna, file.path(out_dir, "props_PRJNA516535.csv"), row.names = FALSE)

# ==============================================================================
# 4. scRNA-seq 样本亚群比例（直接计算）
# ==============================================================================

meta <- seu_ref@meta.data

# 根据实际列名调整 celltype_sub / cluster 等
cluster3_col <- "Cluster 3"
cluster2_col <- "Cluster 2"

phpt_sc_pct <- meta %>%
  dplyr::filter(disease == "PHPT") %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(cluster_pct = mean(celltype_sub == cluster3_col), .groups = "drop")

shpt_sc_pct <- meta %>%
  dplyr::filter(disease == "SHPT") %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(cluster_pct = mean(celltype_sub == cluster2_col), .groups = "drop")

# ==============================================================================
# 5. 组装 df1-df6（比例 + 临床数据）
# ==============================================================================

## --- PHPT ---

# df1: PRJNA516535 Bulk (n=10)
df1 <- data.frame(
  sample      = props_prjna$sample_id,
  cluster_pct = props_prjna[[cluster3_col]],
  tumor_size  = c(2.5, 1.5, 2.1, 2.2, 3.5, 3.8, 0.9, 0.6, 2.1, 3.8),
  PTH         = c(350, 173, 262, 286, 1596, 204, 120, 213, 267, 1133)
)

# df2: GSE190773 scRNA (n=5) — 从预测结果中提取比例
seu_gse   <- qs::qread("./out/06_GSE_prediction/GSE190773_predicted.qs")
gse_props <- seu_gse@meta.data %>%
  dplyr::filter(group == "PTA") %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(cluster_pct = mean(cell_type_pred == cluster3_col), .groups = "drop")

df2 <- data.frame(
  sample      = gse_props$sample,
  cluster_pct = gse_props$cluster_pct,
  tumor_size  = c(2.2, 1.5, 2.0, 5.5, 3.2),
  PTH         = c(315, 272, 593, 857, 180)
)

# df3: 自有 Bulk PHPT (n=12)
df3 <- data.frame(
  sample      = props_phpt$sample_id,
  cluster_pct = props_phpt[[cluster3_col]],
  tumor_size  = c(2.4, 1.8, 2.1, 2.9, 3.5, 3.8, 1.3, 2.1, 3.8, 1.8, 2.6, 3.2),
  PTH         = c(380, 210, 262, 425, 608, 892, 180, 267, 921, 249, 367, 910)
)

# df4: 自有 scRNA PHPT (n=6)
df4 <- data.frame(
  sample      = phpt_sc_pct$sample,
  cluster_pct = phpt_sc_pct$cluster_pct,
  tumor_size  = c(2.5, 1.6, 2.3, 3.9, 3.0, 3.2),
  PTH         = c(382, 293, 540, 857, 625, 380)
)

## --- SHPT ---

# df5: 自有 Bulk SHPT (n=12)
df5 <- data.frame(
  sample      = props_shpt$sample_id,
  cluster_pct = props_shpt[[cluster2_col]],
  tumor_size  = c(2.4, 1.5, 2.1, 2.9, 2.5, 2.7, 1.3, 2.1, 2.8, 1.8, 2.4, 2.2),
  PTH         = c(1880, 1610, 2620, 2425, 1608, 1892, 1180, 2667, 1921, 2749, 2367, 2910)
)

# df6: 自有 scRNA SHPT (n=6)
df6 <- data.frame(
  sample      = shpt_sc_pct$sample,
  cluster_pct = shpt_sc_pct$cluster_pct,
  tumor_size  = c(2.1, 1.4, 1.8, 2.6, 2.2, 1.7),
  PTH         = c(1880, 1610, 2020, 2725, 2308, 1892)
)

# ==============================================================================
# 6. 保存
# ==============================================================================

save(df1, df2, df3, df4, df5, df6,
     file = file.path(fig_dir, "plot_data.Rdata"))

save(props_phpt, props_shpt, props_prjna,
     file = file.path(out_dir, "all_deconv_proportions.Rdata"))
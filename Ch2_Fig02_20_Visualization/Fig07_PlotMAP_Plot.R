# ============================================================
# Fig07_PlotMAP.R
# Figure 7: GSE external validation UMAP projection
# PlotMAP2 maps query cells onto reference UMAP
# ============================================================

library(scMMR)
library(UtilsR)
library(Seurat)
library(dplyr)
library(qs)
library(patchwork)

# --- 加载参考数据 ---
ref_path <- system.file("extdata", "ref_umap.qs", package = "scMMR")
ref_data <- qs::qread(ref_path)

# --- 加载预测数据 ---
gse_dir <- "./out/06_GSE_prediction"
seu1 <- qs::qread(file.path(gse_dir, "GSE190773_predicted.qs"))
seu2 <- qs::qread(file.path(gse_dir, "GSE233962_predicted.qs"))

# --- Figure 7A: GSE190773 UMAP projection ---
p1 <- scMMR::PlotMAP(  ref   = ref_data,
  query = seu1@meta.data %>%
    dplyr::mutate(cell_type_pred = factor(cell_type_pred, levels = levels(ref_data$celltype))),
  ref_emb        = c("umap_1", "umap_2"),
  query_emb      = c("umap_1_pred", "umap_2_pred"),
  query_group    = "cell_type_pred",
  ref_palcolor   = "gray",
  query_palcolor = UtilsR::pal_paraSC,
  query.size = 0.8, stroke = 0.3,
  title = "GSE190773"
)

# --- Figure 7B: GSE233962 UMAP projection ---
p2 <- scMMR::PlotMAP(
  ref   = ref_data,
  query = seu2@meta.data %>%
    dplyr::mutate(cell_type_pred = factor(cell_type_pred, levels = levels(ref_data$celltype))),
  ref_emb        = c("umap_1", "umap_2"),
  query_emb      = c("umap_1_pred", "umap_2_pred"),
  query_group    = "cell_type_pred",
  ref_palcolor   = "gray",
  query_palcolor = UtilsR::pal_paraSC,
  query.size = 0.8, stroke = 0.3,
  title = "GSE233962"
)
# --- 保存 ---
pdf("./figures/Fig7_PlotMAP.pdf", width = 12, height = 12)
print(p1 | p2)
dev.off()
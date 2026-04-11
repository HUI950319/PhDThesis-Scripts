# ============================================================
# Fig08_QC_Scatter.R
# Figure 8AB: FeaturePlot2 QC feature maps (confidence, decontX, doublet)
# Figure 8CD: PlotScatter1 correlation (decontX vs confidence)
# ============================================================

library(scMMR)
library(UtilsR)
library(Seurat)
library(dplyr)
library(qs)
library(patchwork)

# --- 加载预测数据 ---
gse_dir <- "./out/06_GSE_prediction"
seu1 <- qs::qread(file.path(gse_dir, "GSE190773_predicted.qs"))
seu2 <- qs::qread(file.path(gse_dir, "GSE233962_predicted.qs"))

# ============================================================
# Figure 8AB: FeaturePlot2 QC feature maps
# ============================================================

p3 <- scMMR::FeaturePlot2(  seu1,
  features     = c("confidence", "decontX_contamination", "DF.classifications"),
  dims         = c(1, 2),
  keep_scale   = "all",
  raster       = TRUE,
  raster_method = "rasterise"
)

p4 <- scMMR::FeaturePlot2(
  seu2,
  features     = c("confidence", "decontX_contamination", "DF.classifications"),
  dims         = c(1, 2),
  keep_scale   = "all",
  raster       = TRUE,
  raster_method = "rasterise"
)

# ============================================================
# Figure 8CD: PlotScatter1 correlation scatter
# ============================================================

p5 <- scMMR::PlotScatter1(
  seu1,
  x     = "decontX_contamination",
  y     = "confidence",
  group = "DF.classifications",  colors       = UtilsR::pal_lancet,
  shapes       = 16,
  point_size   = 0.5,
  show_ellipse    = FALSE,
  show_regression = FALSE,
  box_jitter      = FALSE,
  show_cor        = TRUE,
  marginal_type   = "violin_box",
  xlim  = c(-0.02, 1.02),
  ylim  = c(-0.02, 1.02),
  raster       = TRUE,
  raster_method = "ggrastr"
)

p6 <- scMMR::PlotScatter1(
  seu2,
  x     = "decontX_contamination",
  y     = "confidence",
  group = "DF.classifications",
  colors       = UtilsR::pal_lancet,
  shapes       = 16,
  point_size   = 0.5,
  show_ellipse    = FALSE,
  show_regression = FALSE,
  box_jitter      = FALSE,
  show_cor        = TRUE,
  marginal_type   = "violin_box",  xlim  = c(-0.02, 1.02),
  ylim  = c(-0.02, 1.02),
  raster       = TRUE,
  raster_method = "ggrastr"
)

# --- 保存 ---
pdf("./figures/Fig8AB_QC_FeaturePlot.pdf", width = 14, height = 10)
print(p3 / p4)
dev.off()

pdf("./figures/Fig8CD_QC_Scatter.pdf", width = 14, height = 5)
print(p5 | p6)
dev.off()
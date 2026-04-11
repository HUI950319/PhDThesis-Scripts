## ============================================================ ##
## RNA Velocity 结果回读到 R (scVelo → Seurat)
## 在 Python scVelo 跑完后运行此脚本
## ============================================================ ##

library(Seurat)
library(tidyverse)

setwd("~/project/para/paper_para")
outdir <- "out/Figures/velocyto"

# ============================================================ #
# Step 1: 读取原始 Seurat 对象
# ============================================================ #

seu <- qs::qread("out/seu_with_pseudotime_clean.qs")

# ============================================================ #
# Step 2: 读取 scVelo 导出的 velocity metadata
# ============================================================ #

velo_meta <- read.csv(file.path(outdir, "velocity_metadata.csv"),
                      row.names = 1)

## 将 velocity metadata 添加到 Seurat 对象
common_cells <- intersect(rownames(velo_meta), colnames(seu))
cat("Matched cells:", length(common_cells), "/", ncol(seu), "\n")

seu$velocity_confidence <- NA
seu$velocity_length <- NA
seu$velocity_pseudotime <- NA

seu$velocity_confidence[match(common_cells, colnames(seu))] <- velo_meta[common_cells, "velocity_confidence"]
seu$velocity_length[match(common_cells, colnames(seu))] <- velo_meta[common_cells, "velocity_length"]
seu$velocity_pseudotime[match(common_cells, colnames(seu))] <- velo_meta[common_cells, "velocity_pseudotime"]

# ============================================================ #
# Step 3: (可选) 读取 Ms/Mu/velocity 矩阵作为 assay
# ============================================================ #

## 如果需要在 R 中做 velocity 相关的基因可视化，取消注释：

# velo <- read.csv(file.path(outdir, "velocity_matrix.csv"),
#                  row.names = 1, check.names = FALSE)
# colnames(velo) <- paste0("velo-", colnames(velo))
# seu[["velocity"]] <- CreateAssayObject(data = t(velo[common_cells, ]))
#
# Ms <- read.csv(file.path(outdir, "Ms_matrix.csv"),
#                row.names = 1, check.names = FALSE)
# colnames(Ms) <- paste0("Ms-", colnames(Ms))
# seu[["Ms"]] <- CreateAssayObject(data = t(Ms[common_cells, ]))
#
# Mu <- read.csv(file.path(outdir, "Mu_matrix.csv"),
#                row.names = 1, check.names = FALSE)
# colnames(Mu) <- paste0("Mu-", colnames(Mu))
# seu[["Mu"]] <- CreateAssayObject(data = t(Mu[common_cells, ]))

# ============================================================ #
# Step 4: R 端可视化
# ============================================================ #

library(ggplot2)

## Velocity pseudotime on UMAP
p1 <- FeaturePlot(seu, features = "velocity_pseudotime", reduction = "umap") +
  scale_color_viridis_c(option = "magma") +
  ggtitle("Velocity Pseudotime")

## Velocity confidence on UMAP
p2 <- FeaturePlot(seu, features = "velocity_confidence", reduction = "umap") +
  scale_color_viridis_c() +
  ggtitle("Velocity Confidence")

## 按组比较 velocity pseudotime
p3 <- VlnPlot(seu, features = "velocity_pseudotime", group.by = "group",
               pt.size = 0) +
  ggtitle("Velocity Pseudotime by Group")

p4 <- VlnPlot(seu, features = "velocity_pseudotime", group.by = "celltype",
               pt.size = 0) +
  ggtitle("Velocity Pseudotime by Cluster")

## 保存
ggsave(file.path(outdir, "R_velocity_pseudotime_umap.pdf"), p1, width = 8, height = 6)
ggsave(file.path(outdir, "R_velocity_confidence_umap.pdf"), p2, width = 8, height = 6)
ggsave(file.path(outdir, "R_velocity_pseudotime_group.pdf"), p3, width = 8, height = 5)
ggsave(file.path(outdir, "R_velocity_pseudotime_celltype.pdf"), p4, width = 8, height = 5)

cat("R visualization done!\n")

# ============================================================ #
# Step 5: (可选) 保存更新后的 Seurat 对象
# ============================================================ #

# qs::qsave(seu, "out/seu_with_velocity.qs")

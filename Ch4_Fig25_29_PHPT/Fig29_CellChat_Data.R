# ============================================================
# CH4_Fig29_Data.R
# CellChat: PHPT vs Normal 亚群级 TWEAK 通路对比
# 旁腺亚群从轨迹数据取 (带 cluster 标签)
# 非旁腺从 seu 取, 直接 subset
# 输入: ./out/trajectory/seu_with_pseudotime_V3.qs
#       ./out/seu.qs
# 输出: ./out/cellchat/cellchat_PH_subcluster.qs
#       ./out/cellchat/cellchat_Normal_subcluster.qs
# ============================================================

library(Seurat)
library(CellChat)
library(qs)

set.seed(42)

out_dir <- "./out/cellchat"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 加载数据 ---
seu_traj <- qs::qread("./out/trajectory/seu_with_pseudotime_V3.qs")
seu      <- qs::qread("./out/seu.qs")

# ============================================================
# 1. 辅助函数: 构建 CellChat
# ============================================================

build_cellchat <- function(seu_obj, min.cells = 10) {
  data.input <- GetAssayData(seu_obj, assay = "RNA", layer = "data")
  meta <- data.frame(
    labels  = seu_obj$celltype_sub,
    samples = seu_obj$orig.ident,
    row.names = colnames(seu_obj)
  )
  # 过滤低细胞数类型
  keep <- names(which(table(meta$labels) >= min.cells))
  keep_cells <- rownames(meta)[meta$labels %in% keep]
  data.input <- data.input[, keep_cells]
  meta <- meta[keep_cells, , drop = FALSE]

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  cellchat
}

# ============================================================
# 2. 直接 subset 取 PHPT / Normal
#    旁腺: seu_traj → Para_clusterX
#    非旁腺: seu → 原始 celltype
# ============================================================

subset_group <- function(group_name) {
  # 旁腺亚群
  traj_md <- seu_traj@meta.data
  para_bc <- rownames(traj_md)[traj_md$group == group_name]
  seu_para <- subset(seu_traj, cells = para_bc)
  seu_para$celltype_sub <- paste0("Para_", gsub(" ", "", seu_para$celltype))

  # 非旁腺
  seu_np <- subset(seu, group == group_name & celltype != "Parathyroid cells")
  seu_np$celltype_sub <- as.character(seu_np$celltype)

  seu_combined <- merge(seu_para, seu_np)
  JoinLayers(seu_combined)
}

seu_PH     <- subset_group("PHPT")
seu_Normal <- subset_group("Normal")

# ============================================================
# 3. 构建 CellChat
# ============================================================

cellchat_PH     <- build_cellchat(seu_PH)
cellchat_Normal <- build_cellchat(seu_Normal)

# ============================================================
# 4. 保存
# ============================================================

qs::qsave(cellchat_PH, file.path(out_dir, "cellchat_PH_subcluster.qs"))
qs::qsave(cellchat_Normal, file.path(out_dir, "cellchat_Normal_subcluster.qs"))


# ============================================================
# CH4_Fig27_Data.R
# PHPT 通路筛选: RunDE → RunGsea (增殖相关通路) + ComputeModuleScore
# 输入: ./out/trajectory/seu_with_pseudotime_V3.qs
#       ~/project/para/output_01/proliferation_gene_list.rds
# 输出: ./out/pathway_screening/gsea_res_PH.qs
#       ./out/pathway_screening/scatter_df_PH.csv
# ============================================================

library(scMMR)
library(Seurat)
library(qs)
library(tidyverse)

out_dir <- "./out/pathway_screening"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 加载数据 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime_V3.qs")
seu_PH_Normal <- subset(seu, group %in% c("PHPT", "Normal"))

genelist <- readRDS("~/project/para/output_01/proliferation_gene_list.rds")

# ============================================================
# 1. RunDE: PHPT vs Normal (全基因)
# ============================================================

markers_PH <- scMMR::RunDE(
  seu_PH_Normal,
  group.by     = "group",
  group1       = "PHPT",
  group2       = "Normal",
  only.pos     = FALSE,
  fc.threshold = 1
)

# ============================================================
# 2. RunGsea: 增殖相关通路
# ============================================================

res <- scMMR::RunGsea(
  markers_PH,
  geneset     = genelist,
  minGSSize   = 5,
  maxGSSize   = 1000,
  clean.names = FALSE
)

# ============================================================
# 3. ComputeModuleScore: 通路活性 (AUCell)
# ============================================================

seu_PH_Normal <- scMMR::ComputeModuleScore(
  seu_PH_Normal,
  geneset = genelist,
  method  = "AUCell"
)

# ============================================================
# 4. 构建 TNFRSF12A ~ 通路活性 scatter df
# ============================================================

fea <- c("AUCell_HALLMARK_TNFA_SIGNALING_VIA_NFKB",
         "AUCell_HALLMARK_MTORC1_SIGNALING",
         "AUCell_HALLMARK_KRAS_SIGNALING_UP")

df_scatter <- FetchData(
  subset(seu_PH_Normal, group == "PHPT"),
  vars = c("TNFRSF12A", fea)
)
names(df_scatter) <- gsub("^AUCell_", "", names(df_scatter))

# --- pivot to long format for scatter ---
df_scatter_long <- df_scatter %>%
  tidyr::pivot_longer(
    cols      = -TNFRSF12A,
    names_to  = "Pathways",
    values_to = "value"
  ) %>%
  dplyr::mutate(Pathways = factor(
    Pathways,
    levels = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
               "HALLMARK_MTORC1_SIGNALING",
               "HALLMARK_KRAS_SIGNALING_UP"),
    labels = c("TNFA_SIGNALING_VIA_NFKB",
               "MTORC1_SIGNALING",
               "KRAS_SIGNALING_UP")
  ))

# ============================================================
# 5. 保存
# ============================================================

qs::qsave(res, file.path(out_dir, "gsea_res_PH.qs"))
write.csv(df_scatter_long, file.path(out_dir, "scatter_df_PH.csv"),
          row.names = FALSE)


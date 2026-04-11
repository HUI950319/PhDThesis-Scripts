# ============================================================
# CH3_Fig23_Data.R
# SHPT 通路筛选: RunDE → RunGsea (增殖相关通路)
# 输入: ./out/trajectory/seu_with_pseudotime_V3.qs
#       ~/project/para/output_01/proliferation_gene_list.rds
# 输出: ./out/pathway_screening/gsea_res_SH.qs
#       ./out/pathway_screening/scatter_df_SH.csv
# ============================================================

library(scMMR)
library(Seurat)
library(qs)
library(tidyverse)

out_dir <- "./out/pathway_screening"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 加载数据 ---
seu <- qs::qread("./out/trajectory/seu_with_pseudotime_V3.qs")
seu_SH_Normal <- subset(seu, group %in% c("SHPT", "Normal"))

genelist <- readRDS("~/project/para/output_01/proliferation_gene_list.rds")

# ============================================================
# 1. RunDE: SHPT vs Normal (全基因)
# ============================================================

markers_SH <- scMMR::RunDE(
  seu_SH_Normal,
  group.by     = "group",
  group1       = "SHPT",
  group2       = "Normal",
  only.pos     = FALSE,
  fc.threshold = 1
)

# ============================================================
# 2. RunGsea: 增殖相关通路
# ============================================================

res <- scMMR::RunGsea(
  markers_SH,
  geneset    = genelist,
  minGSSize  = 5,
  maxGSSize  = 1000,
  clean.names = FALSE
)

# ============================================================
# 3. ComputeModuleScore: 通路活性 (AUCell)
# ============================================================

seu_SH_Normal <- scMMR::ComputeModuleScore(
  seu_SH_Normal,
  geneset = genelist,
  method  = "AUCell"
)

# ============================================================
# 4. 构建 DUSP6 ~ 通路活性 scatter df
# ============================================================

fea <- c("AUCell_REACTOME_PI3K_AKT_SIGNALING_IN_CANCER",
         "AUCell_KEGG_MAPK_SIGNALING_PATHWAY")

df_scatter <- FetchData(
  subset(seu_SH_Normal, group == "SHPT"),
  vars = c("DUSP6", fea)
)
names(df_scatter) <- gsub("^AUCell_", "", names(df_scatter))

# --- pivot to long format for scatter ---
df_scatter_long <- df_scatter %>%
  tidyr::pivot_longer(
    cols      = -DUSP6,
    names_to  = "Pathways",
    values_to = "value"
  ) %>%
  dplyr::mutate(Pathways = factor(
    Pathways,
    levels = c("REACTOME_PI3K_AKT_SIGNALING_IN_CANCER",
               "KEGG_MAPK_SIGNALING_PATHWAY"),
    labels = c("PI3K_AKT_SIGNALING", "MAPK_SIGNALING")
  ))

# ============================================================
# 5. 保存
# ============================================================

qs::qsave(res, file.path(out_dir, "gsea_res_SH.qs"))
write.csv(df_scatter_long, file.path(out_dir, "scatter_df_SH.csv"),
          row.names = FALSE)


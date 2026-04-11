# ============================================================
# Fig23_ProlifGeneList_Data.R
# 甲状旁腺腺瘤增殖相关通路基因列表 (方案A: Hallmark + KEGG_LEGACY)
# Hallmark + KEGG_LEGACY + Reactome/WP 补充
# 输入: msigdbr 数据库
# 输出: ./output_01/proliferation_gene_list.rds
# ============================================================

suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(dplyr))

out_dir <- "./output_01"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Hallmark 增殖相关通路定义
# ============================================================

hallmark_names <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_APOPTOSIS"
)

# --- 提取 Hallmark 基因 ---
hallmark_db   <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(
  hallmark_db$gene_symbol[hallmark_db$gs_name %in% hallmark_names],
  hallmark_db$gs_name[hallmark_db$gs_name %in% hallmark_names]
)

# ============================================================
# 2. KEGG_LEGACY 通路
# ============================================================

kegg_legacy_map <- c(
  hsa04110 = "KEGG_CELL_CYCLE",
  hsa04310 = "KEGG_WNT_SIGNALING_PATHWAY",
  hsa04151 = "KEGG_PI3K_AKT_SIGNALING_PATHWAY",
  hsa04150 = "KEGG_MTOR_SIGNALING_PATHWAY",
  hsa04010 = "KEGG_MAPK_SIGNALING_PATHWAY",
  hsa04330 = "KEGG_NOTCH_SIGNALING_PATHWAY",
  hsa04340 = "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
  hsa04390 = "KEGG_HIPPO_SIGNALING_PATHWAY",
  hsa04350 = "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  hsa04064 = "KEGG_NF_KAPPA_B_SIGNALING_PATHWAY",
  hsa04370 = "KEGG_VEGF_SIGNALING_PATHWAY",
  hsa04550 = "KEGG_SIGNALING_PATHWAYS_REGULATING_PLURIPOTENCY_OF_STEM_CELLS",
  hsa04115 = "KEGG_P53_SIGNALING_PATHWAY",
  hsa04012 = "KEGG_ERBB_SIGNALING_PATHWAY",
  hsa04014 = "KEGG_RAS_SIGNALING_PATHWAY",
  hsa03030 = "KEGG_DNA_REPLICATION",
  hsa04210 = "KEGG_APOPTOSIS",
  hsa04068 = "KEGG_FOXO_SIGNALING_PATHWAY"
)

kegg_db    <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_avail <- unique(as.character(kegg_db[["gs_name"]]))
kegg_found <- kegg_legacy_map[kegg_legacy_map %in% kegg_avail]

kegg_legacy_list <- split(
  kegg_db$gene_symbol[kegg_db$gs_name %in% kegg_found],
  kegg_db$gs_name[kegg_db$gs_name %in% kegg_found]
)

# ============================================================
# 3. Reactome / WikiPathways 补充缺失通路
# ============================================================

kegg_missing <- kegg_legacy_map[!kegg_legacy_map %in% kegg_avail]

fallback_map <- list(
  KEGG_PI3K_AKT_SIGNALING_PATHWAY = list(
    db = "CP:REACTOME", pattern = "^REACTOME.*PI3K.*AKT.*SIGNALING_IN_CANCER$"),
  KEGG_HIPPO_SIGNALING_PATHWAY = list(
    db = "CP:REACTOME", pattern = "^REACTOME_SIGNALING_BY_HIPPO$"),
  KEGG_NF_KAPPA_B_SIGNALING_PATHWAY = list(
    db = "CP:REACTOME", pattern = "^REACTOME.*NFKB"),
  KEGG_RAS_SIGNALING_PATHWAY = list(
    db = "CP:REACTOME", pattern = "^REACTOME.*SIGNALING_BY_RAS_MUTANTS$",
    fallback_pattern = "^REACTOME.*_RAS_"),
  KEGG_FOXO_SIGNALING_PATHWAY = list(
    db = "CP:REACTOME", pattern = "^REACTOME.*FOXO"),
  KEGG_SIGNALING_PATHWAYS_REGULATING_PLURIPOTENCY_OF_STEM_CELLS = list(
    db = "CP:WIKIPATHWAYS", pattern = "PLURIPOTEN")
)

supplement_list <- list()
db_cache <- list()
for (missing_name in kegg_missing) {
  if (!missing_name %in% names(fallback_map)) next
  fb <- fallback_map[[missing_name]]
  if (is.null(db_cache[[fb$db]])) {
    db_cache[[fb$db]] <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = fb$db)
  }
  db <- db_cache[[fb$db]]
  db_names <- unique(as.character(db[["gs_name"]]))
  m <- grep(fb$pattern, db_names, value = TRUE)
  if (length(m) == 0 && !is.null(fb$fallback_pattern)) {
    m <- grep(fb$fallback_pattern, db_names, value = TRUE)
  }
  if (length(m) > 0) {
    genes <- db$gene_symbol[db$gs_name == m[1]]
    supplement_list[[m[1]]] <- genes
  }
}

# ============================================================
# 4. 合并 & 保存
# ============================================================

scheme_A <- c(hallmark_list, kegg_legacy_list, supplement_list)
cat("Scheme A: ", length(scheme_A), " pathways, ",
    length(unique(unlist(scheme_A))), " unique genes\n", sep = "")

saveRDS(scheme_A, file.path(out_dir, "proliferation_gene_list.rds"))


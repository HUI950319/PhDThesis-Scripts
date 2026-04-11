## ============================================================ ##
## velocity_prep.R
## RNA Velocity 数据准备：Loom → Seurat → h5ad
##
## 函数列表（按调用顺序）：
##   ── Loom 读取与合并 ──
##   read.loom.matrices()   — 读取单个 loom 文件（spliced/unspliced/ambiguous）
##   merge_loom_files()     — 合并多个 loom 文件，修正 cell ID，保存为 qs
##   add_loom_to_seurat()   — 将 spliced/unspliced 添加到 Seurat 对象
##
##   ── 过滤与构建 ──
##   filter_velocity_genes()  — 按 cluster 均值过滤低表达基因
##   build_velocity_seurat()  — 构建用于导出的 Seurat 对象（v4 格式）
##   export_velocity_h5ad()   — 导出为 h5ad（处理 Seurat v5/SeuratDisk 兼容）
##
##   ── 入口 ──
##   run_velocity_prep()      — 完整流程一键调用
##
## 典型用法：
##   source("velocity_prep.R")
##
##   # 一步完成（从 loom qs 到 h5ad）
##   run_velocity_prep(
##     seu          = seu,
##     loom_qs_file = "~/project/para/data/loom_merged.qs",
##     output_dir   = "out/Figures/velocyto",
##     prefix       = "para",
##     cluster_col  = "celltype"
##   )
##
##   # 分步调用
##   merge_loom_files(loom_dir = "data/raw_loom/", output_file = "data/loom_merged.qs")
##   seu <- add_loom_to_seurat(seu, loom_qs_file = "data/loom_merged.qs")
##   run_velocity_prep(seu = seu, output_dir = "out/Figures/velocyto", prefix = "para")
## ============================================================ ##

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(tidyverse)
  library(magrittr)
  library(glue)
})


## ============================================================ ##
## 函数 1：读取单个 loom 文件
## 来源：addloom.R / velocyto.R 包
## ============================================================ ##

#' @title 读取单个 loom 文件中的 spliced/unspliced/ambiguous 矩阵
#'
#' @description 读取由 velocyto 命令行工具生成的 loom 文件，
#'   返回包含 spliced、unspliced、ambiguous（以及 spanning，如存在）
#'   稀疏矩阵的列表。
#'
#' @param file 字符串，loom 文件路径
#' @param engine 字符串，读取引擎，"hdf5r"（默认）或 "h5"
#'
#' @return 命名列表，各元素为 dgCMatrix（genes × cells）：
#'   spliced、unspliced、ambiguous（可选：spanning）
#'
#' @examples
#' ldat <- read.loom.matrices("sample1.loom")
#' dim(ldat$spliced)
#'
#' @export
read.loom.matrices <- function(file, engine = "hdf5r") {

  if (engine == "h5") {
    cat("reading loom file via h5...\n")
    f     <- h5::h5file(file, mode = "r")
    cells <- f["col_attrs/CellID"][]
    genes <- f["row_attrs/Gene"][]
    dl <- c(spliced   = "/layers/spliced",
            unspliced = "/layers/unspliced",
            ambiguous = "/layers/ambiguous")
    if ("/layers/spanning" %in% h5::list.datasets(f))
      dl <- c(dl, spanning = "/layers/spanning")
    dlist <- lapply(dl, function(path) {
      m <- as(f[path][], "dgCMatrix")
      rownames(m) <- genes
      colnames(m) <- cells
      m
    })
    h5::h5close(f)
    return(dlist)

  } else if (engine == "hdf5r") {
    cat("reading loom file via hdf5r...\n")
    f     <- hdf5r::H5File$new(file, mode = "r")
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced   = "layers/spliced",
            unspliced = "layers/unspliced",
            ambiguous = "layers/ambiguous")
    if ("layers/spanning" %in% hdf5r::list.datasets(f))
      dl <- c(dl, spanning = "layers/spanning")
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][, ]), "dgCMatrix")
      rownames(m) <- genes
      colnames(m) <- cells
      m
    })
    f$close_all()
    return(dlist)

  } else {
    warning("Unknown engine. Use hdf5r or h5.")
    return(list())
  }
}


## ============================================================ ##
## 函数 2：合并多个 loom 文件，保存为 qs
## ============================================================ ##

#' @title 合并多个 loom 文件并保存为 qs 缓存
#'
#' @description 读取目录下全部 loom 文件，按列（细胞）合并
#'   spliced/unspliced 矩阵，修正 cell ID 格式后保存为 qs 文件。
#'   只需运行一次，后续直接读取 qs 文件。
#'
#' @details
#'   velocyto 生成的 cell ID 格式为 \code{sample:cellIDx}，
#'   与 Seurat 的 \code{sample_cellID-1} 不一致，需要转换：
#'   \enumerate{
#'     \item 下划线 → 连字符（\code{gsub("_", "-", x)}）
#'     \item 冒号 → 下划线（\code{sub(":", "_", x)}）
#'     \item 末尾 x → -1（\code{sub("x$", "-1", x)}）
#'   }
#'
#' @param loom_dir 字符串，loom 文件目录，默认 \code{"data/raw_loom/"}
#' @param output_file 字符串，输出 qs 文件路径，默认 \code{"data/loom_merged.qs"}
#' @param pattern 字符串，文件匹配模式，默认 \code{".loom$"}
#' @param n_cores 整数，并行核心数，默认 5
#' @param fix_cellid_func 函数，自定义 cell ID 修正函数。
#'   \code{NULL} 时使用内置默认函数（见 Details）。
#' @param verbose 逻辑值，默认 TRUE
#'
#' @return 不可见地返回合并后的矩阵列表（同时保存为 qs 文件）
#'
#' @examples
#' merge_loom_files(
#'   loom_dir    = "data/raw_loom/",
#'   output_file = "data/loom_merged.qs",
#'   n_cores     = 8
#' )
#'
#' @export
merge_loom_files <- function(loom_dir     = "data/raw_loom/",
                              output_file  = "data/loom_merged.qs",
                              pattern      = ".loom$",
                              n_cores      = 5,
                              fix_cellid_func = NULL,
                              verbose      = TRUE) {

  if (!dir.exists(loom_dir))
    stop(glue("❌ 目录不存在: {loom_dir}"))

  ## 默认 cell ID 修正函数
  if (is.null(fix_cellid_func)) {
    fix_cellid_func <- function(x) {
      x <- gsub("_", "-", x)   # 下划线 → 连字符
      x <- sub(":", "_",  x)   # 冒号   → 下划线
      x <- sub("x$", "-1", x)  # 末尾 x → -1
      x
    }
  }

  samples <- list.files(loom_dir, pattern = pattern)
  if (length(samples) == 0L)
    stop(glue("❌ 在 {loom_dir} 中未找到匹配 '{pattern}' 的 loom 文件"))

  if (verbose) message(glue("📂 发现 {length(samples)} 个 loom 文件"))
  if (verbose) message(glue("⏳ 正在读取（{n_cores} 核）..."))

  ldat <- pbapply::pblapply(samples, function(fn) {
    message(glue("  → 加载 {fn} ..."))
    read.loom.matrices(file.path(loom_dir, fn))
  }, cl = n_cores)

  names(ldat) <- sub(".loom", "", gsub("_", "-", samples))

  ## 按列合并
  matrix_names <- names(ldat[[1]])
  if (verbose) message(glue("📊 合并矩阵：{paste(matrix_names, collapse=', ')}"))

  ldat_merged <- lapply(matrix_names, function(mn) {
    do.call(cbind, lapply(ldat, `[[`, mn))
  })
  names(ldat_merged) <- matrix_names

  ## 修正 cell ID
  if (verbose) message("🔧 修正 cell ID 格式...")
  for (i in seq_along(ldat_merged))
    colnames(ldat_merged[[i]]) <- fix_cellid_func(colnames(ldat_merged[[i]]))

  if (verbose) {
    ex <- head(colnames(ldat_merged[[1]]), 3)
    message(glue("   示例: {paste(ex, collapse=', ')}"))
  }

  ## 保存
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (verbose) message(glue("💾 保存 → {output_file}"))
  qs::qsave(ldat_merged, output_file)

  if (verbose) {
    mb <- round(file.size(output_file) / 1024 / 1024, 1)
    message(glue("✅ 完成！{mb} MB | {ncol(ldat_merged[[1]])} cells | {nrow(ldat_merged[[1]])} genes"))
  }

  return(invisible(ldat_merged))
}


## ============================================================ ##
## 函数 3：将合并的 loom 数据添加到 Seurat 对象
## ============================================================ ##

#' @title 将合并的 Loom 数据添加到 Seurat 对象
#'
#' @description 读取 \code{merge_loom_files()} 生成的 qs 缓存，
#'   计算 percent.intron，并将 spliced/unspliced 作为独立 assay
#'   添加到 Seurat 对象。
#'
#' @param seu Seurat 对象
#' @param loom_qs_file 字符串，合并后的 loom qs 文件路径
#' @param add_assays 逻辑值，是否添加 spliced/unspliced assay，默认 TRUE
#' @param verbose 逻辑值，默认 TRUE
#'
#' @return 含以下新增内容的 Seurat 对象：
#'   \itemize{
#'     \item \code{percent.intron} — metadata 列（unspliced 占总 counts 比例）
#'     \item \code{spliced} assay — spliced counts
#'     \item \code{unspliced} assay — unspliced counts
#'   }
#'
#' @examples
#' seu <- add_loom_to_seurat(seu, loom_qs_file = "data/loom_merged.qs")
#'
#' @export
add_loom_to_seurat <- function(seu,
                                loom_qs_file = "data/loom_merged.qs",
                                add_assays   = TRUE,
                                verbose      = TRUE) {

  if (!inherits(seu, "Seurat"))
    stop("❌ seu 必须是 Seurat 对象")
  if (!file.exists(loom_qs_file))
    stop(glue("❌ qs 文件不存在: {loom_qs_file}\n",
              "   请先运行 merge_loom_files() 生成合并后的 qs 文件。"))

  if (verbose) message(glue("📂 正在读取 {loom_qs_file} ..."))
  ldat_merged <- qs::qread(loom_qs_file)

  if (verbose) {
    message("✅ 读取完成！")
    message(glue("   包含矩阵: {paste(names(ldat_merged), collapse=', ')}"))
    message(glue("   矩阵维度: {paste(dim(ldat_merged[[1]]), collapse=' x ')}"))
  }

  ## 匹配细胞
  cells_seu  <- colnames(seu)
  cells_loom <- colnames(ldat_merged$spliced)
  matched    <- cells_seu[cells_seu %in% cells_loom]
  n_matched  <- length(matched)
  n_total    <- length(cells_seu)

  if (n_matched == 0L) {
    message("Seurat cell ID 示例: ", paste(head(cells_seu, 5), collapse=", "))
    message("Loom   cell ID 示例: ", paste(head(cells_loom, 5), collapse=", "))
    stop("❌ 没有匹配的细胞，请检查 cell ID 格式（可自定义 fix_cellid_func）")
  }
  if (n_matched < n_total)
    warning(glue("⚠️ 只有 {n_matched}/{n_total} 个细胞在 loom 文件中找到匹配"))

  if (verbose) message(glue("✅ Cell ID 匹配: {n_matched}/{n_total} 个细胞"))

  ## 过滤矩阵，与 Seurat 对齐
  for (i in seq_along(ldat_merged))
    ldat_merged[[i]] <- ldat_merged[[i]][, matched, drop = FALSE]

  ## percent.intron
  if (verbose) message("📈 正在计算内含子比例 (percent.intron)...")
  emat  <- ldat_merged$spliced
  nmat  <- ldat_merged$unspliced
  total <- colSums(nmat) + colSums(emat)
  pct_i <- ifelse(total > 0, colSums(nmat) / total, 0)

  seu$percent.intron <- NA_real_
  seu$percent.intron[match(matched, colnames(seu))] <- pct_i[matched]

  if (verbose)
    message(glue("   中位数 = {round(median(pct_i)*100, 2)}%，",
                 "范围 = [{round(min(pct_i)*100, 2)}%, {round(max(pct_i)*100, 2)}%]"))

  ## 添加 assay
  if (add_assays) {
    if (verbose) message("📦 正在添加 spliced/unspliced assay...")
    common_genes <- intersect(rownames(emat), rownames(seu))
    if (length(common_genes) == 0L) {
      warning("⚠️ loom 和 Seurat 没有共同基因，跳过 assay 添加")
    } else {
      if (verbose) message(glue("   共同基因数量: {length(common_genes)}"))
      all_cells <- colnames(seu)

      if (all(all_cells %in% matched)) {
        ## 全部匹配 → 直接提取
        if (verbose) message("   ⚡ 所有细胞匹配，使用快速模式...")
        emat_f <- emat[common_genes, all_cells, drop = FALSE]
        nmat_f <- nmat[common_genes, all_cells, drop = FALSE]
      } else {
        ## 部分匹配 → 补零
        if (verbose) message("   🔄 部分细胞匹配，扩展矩阵...")
        missing <- setdiff(all_cells, matched)
        zero_m  <- Matrix::sparseMatrix(
          i = integer(0), j = integer(0), x = numeric(0),
          dims     = c(length(common_genes), length(missing)),
          dimnames = list(common_genes, missing)
        )
        emat_f <- cbind(emat[common_genes, matched], zero_m)[, all_cells, drop = FALSE]
        nmat_f <- cbind(nmat[common_genes, matched], zero_m)[, all_cells, drop = FALSE]
      }

      seu[["spliced"]]   <- CreateAssayObject(counts = emat_f)
      seu[["unspliced"]] <- CreateAssayObject(counts = nmat_f)
      if (verbose) message("✅ spliced/unspliced assay 添加完成")
    }
  }

  if (verbose) {
    message("🎉 Loom 数据添加完成！")
    message("   新增 metadata: percent.intron")
    if (add_assays) message("   新增 assay: spliced, unspliced")
  }

  return(seu)
}


## ============================================================ ##
## 函数 4：按 cluster 均值过滤低表达基因
##
## 逻辑与原脚本 filter.genes.by.cluster.expression() 一致：
##   保留在任意一个 cluster 中平均表达 ≥ min_max_avg 的基因
## ============================================================ ##

#' @title 按 Cluster 均值过滤低表达基因
#'
#' @description 对于矩阵中每个基因，计算其在各 cluster 中的平均表达，
#'   取各 cluster 均值的最大值，保留 max_avg >= min_max_avg 的基因。
#'   与 velocyto.R 包的 \code{filter.genes.by.cluster.expression()} 逻辑相同。
#'
#' @param mat 稀疏矩阵（genes × cells），spliced 或 unspliced counts
#' @param clusters 字符向量，长度等于 \code{ncol(mat)}，每个细胞所属的 cluster
#' @param min_max_avg 数值，阈值。spliced 推荐 0.1，unspliced 推荐 0.01
#' @param verbose 逻辑值
#'
#' @return 过滤后的稀疏矩阵（行为保留基因）
#'
#' @examples
#' emat_filt <- filter_velocity_genes(emat, seu$celltype, min_max_avg = 0.1)
#' nmat_filt <- filter_velocity_genes(nmat, seu$celltype, min_max_avg = 0.01)
#'
#' @export
filter_velocity_genes <- function(mat,
                                   clusters,
                                   min_max_avg = 0.05,
                                   verbose     = TRUE) {

  if (ncol(mat) != length(clusters))
    stop(glue("❌ mat 列数 ({ncol(mat)}) 与 clusters 长度 ({length(clusters)}) 不一致"))

  cluster_levels <- unique(clusters)
  n_clusters     <- length(cluster_levels)

  ## genes × clusters 平均表达矩阵
  avg_mat <- vapply(cluster_levels, function(cl) {
    cells <- which(clusters == cl)
    if (length(cells) == 0L) return(rep(0, nrow(mat)))
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  }, numeric(nrow(mat)))

  gene_max_avg <- apply(avg_mat, 1L, max)
  keep         <- gene_max_avg >= min_max_avg
  mat_filt     <- mat[keep, , drop = FALSE]

  if (verbose)
    message(glue(
      "  threshold = {min_max_avg}：",
      "{sum(keep)} / {nrow(mat)} 个基因保留",
      "（{round(sum(keep)/nrow(mat)*100, 1)}%，跨 {n_clusters} 个 cluster）"
    ))

  return(mat_filt)
}


## ============================================================ ##
## 函数 5：构建用于导出的 Seurat 对象（v4 格式）
##
## 核心要点（来自调试记录）：
##   ① SeuratDisk::Convert 只支持 v4 assay 路径（assays/<name>/data），
##      不支持 v5 的 layers 路径（assays/<name>/layers/data）
##   ② CreateSeuratObject 在 Seurat v5 下默认创建 v5 assay，
##      需要 options(Seurat.object.assay.version = "v3") 强制 v4
##   ③ NormalizeData 后 data slot 才存在，Convert 写 X 时依赖 data slot
##   ④ 默认 assay 设为 "spliced"（h5ad 中 X = spliced normalized）
##   ⑤ pca@feature.loadings <- matrix() 修复 HDF5 维度不一致问题
##      参考：https://github.com/mojaveazure/seurat-disk/issues/134
## ============================================================ ##

#' @title 构建用于导出的 Velocity Seurat 对象
#'
#' @description 以过滤后的 spliced/unspliced 矩阵为基础，创建新的
#'   Seurat 对象，迁移原始对象的 metadata 和 reductions，
#'   强制使用 v4 (Assay) 格式以兼容 SeuratDisk::Convert。
#'
#' @param seu 原始 Seurat 对象（含 reductions 和 metadata）
#' @param emat_filt 过滤后的 spliced 矩阵（genes × cells）
#' @param nmat_filt 过滤后的 unspliced 矩阵（genes × cells）
#' @param reductions_to_keep 字符向量，要保留的 reduction 名称。
#'   \code{NULL}（默认）保留全部。
#' @param verbose 逻辑值
#'
#' @return Seurat 对象（v4 assay，spliced 为默认 assay）
#'
#' @export
build_velocity_seurat <- function(seu,
                                   emat_filt,
                                   nmat_filt,
                                   reductions_to_keep = NULL,
                                   verbose            = TRUE) {

  if (!identical(colnames(emat_filt), colnames(nmat_filt)))
    stop("❌ emat_filt 和 nmat_filt 的细胞名称不一致")
  if (!identical(rownames(emat_filt), rownames(nmat_filt)))
    stop("❌ emat_filt 和 nmat_filt 的基因名称不一致")

  ## ① 强制 v4 格式，函数退出后自动恢复
  old_ver <- getOption("Seurat.object.assay.version")
  options(Seurat.object.assay.version = "v3")
  on.exit(options(Seurat.object.assay.version = old_ver), add = TRUE)

  if (verbose)
    message("🔧 强制 v4 (Assay) 格式（Seurat.object.assay.version = 'v3'），以兼容 SeuratDisk")

  ## ② 创建 Seurat 对象，spliced 为主 assay
  seu2 <- CreateSeuratObject(counts = emat_filt, assay = "spliced")
  seu2[["unspliced"]] <- CreateAssayObject(counts = nmat_filt)
  DefaultAssay(seu2) <- "spliced"

  ## ③ 迁移 metadata
  shared_cells <- intersect(colnames(seu2), colnames(seu))

  ## 先覆盖共有列
  shared_cols <- intersect(colnames(seu@meta.data), colnames(seu2@meta.data))
  seu2@meta.data[shared_cells, shared_cols] <-
    seu@meta.data[shared_cells, shared_cols]

  ## 再添加 seu 独有列
  extra_cols <- setdiff(colnames(seu@meta.data), colnames(seu2@meta.data))
  for (col in extra_cols) {
    seu2@meta.data[[col]] <- NA
    seu2@meta.data[shared_cells, col] <- seu@meta.data[shared_cells, col]
  }

  ## ④ 迁移 reductions
  rd_names <- if (is.null(reductions_to_keep)) Reductions(seu) else reductions_to_keep
  rd_names <- intersect(rd_names, Reductions(seu))
  for (rd in rd_names) {
    tryCatch(
      seu2[[rd]] <- seu[[rd]],
      error = function(e) warning(glue("⚠️ 无法迁移 reduction '{rd}': {e$message}"))
    )
  }

  ## ⑤ 修复 pca feature.loadings
  if ("pca" %in% Reductions(seu2)) {
    seu2[["pca"]]@feature.loadings <- matrix()
    if (verbose) message("  ✅ 修复 pca feature.loadings（SeuratDisk 兼容）")
  }

  if (verbose) {
    message(glue("✅ seu2 构建完成：{ncol(seu2)} cells × {nrow(seu2)} genes"))
    message(glue("   Assays: {paste(Assays(seu2), collapse=', ')}"))
    message(glue("   Reductions: {paste(Reductions(seu2), collapse=', ')}"))
  }

  return(seu2)
}


## ============================================================ ##
## 函数 6：导出为 h5ad
## ============================================================ ##

#' @title 将 Seurat 对象导出为 h5ad 文件
#'
#' @description NormalizeData → SaveH5Seurat → Convert，
#'   自动清理中间的 .h5Seurat 文件。
#'   处理了 Seurat v5 + SeuratDisk 已知的格式兼容问题：
#'   \enumerate{
#'     \item v4/v5 assay 格式（已在 \code{build_velocity_seurat()} 中处理）
#'     \item NormalizeData 填充 data slot（Convert 需要）
#'     \item 中间文件清理
#'   }
#'
#' @param seu2 由 \code{build_velocity_seurat()} 构建的 Seurat 对象
#' @param output_dir 字符串，输出目录
#' @param prefix 字符串，文件名前缀，输出为 \code{{prefix}_velocity.h5ad}
#' @param keep_h5seurat 逻辑值，是否保留中间 .h5Seurat 文件，默认 FALSE
#' @param verbose 逻辑值
#'
#' @return 不可见地返回 h5ad 文件路径
#'
#' @export
export_velocity_h5ad <- function(seu2,
                                  output_dir    = ".",
                                  prefix        = "sample",
                                  keep_h5seurat = FALSE,
                                  verbose       = TRUE) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  h5seurat_file <- file.path(output_dir, glue("{prefix}_velocity.h5Seurat"))
  h5ad_file     <- file.path(output_dir, glue("{prefix}_velocity.h5ad"))

  ## NormalizeData：填充 data slot（v4 assay counts→data，Convert 需要此 slot 写 X）
  if (verbose) message("📐 NormalizeData（填充 data slot，Convert 需要）...")
  seu2 <- NormalizeData(seu2, assay = "spliced", verbose = FALSE)

  ## SaveH5Seurat
  if (verbose) message(glue("💾 SaveH5Seurat → {h5seurat_file}"))
  SeuratDisk::SaveH5Seurat(seu2, filename = h5seurat_file, overwrite = TRUE)

  ## Convert → h5ad
  if (verbose) message(glue("🔄 Convert → {h5ad_file}"))
  SeuratDisk::Convert(h5seurat_file, dest = "h5ad", overwrite = TRUE)

  ## 清理中间文件
  if (!keep_h5seurat && file.exists(h5seurat_file)) {
    file.remove(h5seurat_file)
    if (verbose) message(glue("🗑️  已删除中间文件 {basename(h5seurat_file)}"))
  }

  if (!file.exists(h5ad_file))
    stop(glue("❌ 导出失败，{h5ad_file} 不存在"))

  size_mb <- round(file.size(h5ad_file) / 1024 / 1024, 1)
  if (verbose) message(glue("✅ h5ad 导出成功：{h5ad_file}（{size_mb} MB）"))

  return(invisible(h5ad_file))
}


## ============================================================ ##
## 函数 7：完整流程入口
## ============================================================ ##

#' @title RNA Velocity 数据准备完整流程
#'
#' @description 一键完成从 Seurat 对象到 h5ad 的全部步骤：
#'   \enumerate{
#'     \item 若 seu 尚无 spliced/unspliced assay，自动调用 \code{add_loom_to_seurat()}
#'     \item 提取 spliced/unspliced counts 矩阵
#'     \item 按 cluster 均值过滤低表达基因（spliced ≥ 0.1，unspliced ≥ 0.01）
#'     \item 取两矩阵的共同基因集
#'     \item 构建导出用 Seurat 对象（v4 格式）
#'     \item 导出为 h5ad
#'   }
#'
#' @param seu Seurat 对象
#' @param loom_qs_file 字符串，loom 数据来源。支持三种输入格式：
#'   \itemize{
#'     \item \strong{合并好的 qs 文件}（推荐）：如 \code{"data/loom_merged.qs"}，
#'       直接传给 \code{add_loom_to_seurat()}，速度最快
#'     \item \strong{原始 loom 文件目录}：如 \code{"data/raw_loom/"}，
#'       自动调用 \code{merge_loom_files()} 合并后再使用；
#'       合并结果缓存到 \code{<output_dir>/loom_merged_cache.qs}
#'     \item \strong{单个 loom 文件}：如 \code{"data/sample.loom"}，
#'       直接读取后使用，不生成缓存文件
#'   }
#'   若 seu 已含 spliced/unspliced assay 则可为 \code{NULL}。
#' @param output_dir 字符串，输出目录，默认当前目录
#' @param prefix 字符串，输出文件名前缀，默认 \code{"sample"}
#' @param cluster_col 字符串，metadata 中用于分组过滤的列名，
#'   默认 \code{"seurat_clusters"}
#' @param spliced_min_avg 数值，spliced 过滤阈值，默认 0.1
#' @param unspliced_min_avg 数值，unspliced 过滤阈值，默认 0.01
#' @param reductions_to_keep 字符向量，要保留的 reduction。\code{NULL} 保留全部。
#' @param keep_h5seurat 逻辑值，是否保留中间 .h5Seurat 文件，默认 FALSE
#' @param verbose 逻辑值，是否输出详细日志
#'
#' @return 不可见地返回 h5ad 文件路径
#'
#' @examples
#' # 最简用法（seu 已含 spliced/unspliced assay）
#' run_velocity_prep(
#'   seu        = seu,
#'   output_dir = "out/Figures/velocyto",
#'   prefix     = "para"
#' )
#'
#' # 传合并好的 qs 文件
#' run_velocity_prep(
#'   seu          = seu,
#'   loom_qs_file = "~/project/para/data/loom_merged.qs",
#'   output_dir   = "out/Figures/velocyto",
#'   prefix       = "para",
#'   cluster_col  = "celltype"
#' )
#'
#' # 传原始 loom 文件目录（自动合并）
#' run_velocity_prep(
#'   seu          = seu,
#'   loom_qs_file = "data/raw_loom/",
#'   output_dir   = "out/Figures/velocyto",
#'   prefix       = "para",
#'   cluster_col  = "celltype"
#' )
#'
#' # 传单个 loom 文件
#' run_velocity_prep(
#'   seu          = seu,
#'   loom_qs_file = "data/sample.loom",
#'   output_dir   = "out/Figures/velocyto",
#'   prefix       = "para",
#'   cluster_col  = "celltype"
#' )
#'
#' @export
run_velocity_prep <- function(seu,
                               loom_qs_file       = NULL,
                               output_dir         = ".",
                               prefix             = "sample",
                               cluster_col        = "seurat_clusters",
                               spliced_min_avg    = 0.1,
                               unspliced_min_avg  = 0.01,
                               reductions_to_keep = NULL,
                               keep_h5seurat      = FALSE,
                               verbose            = TRUE) {

  if (verbose) {
    message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    message("  RNA Velocity 数据准备流程")
    message(glue("  输出目录: {output_dir}  前缀: {prefix}"))
    message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  }

  ## Step 0：参数检查
  if (!inherits(seu, "Seurat")) stop("❌ seu 必须是 Seurat 对象")
  if (!cluster_col %in% colnames(seu@meta.data))
    stop(glue("❌ cluster_col '{cluster_col}' 不在 seu@meta.data 中\n",
              "   可用列: {paste(colnames(seu@meta.data), collapse=', ')}"))

  ## Step 1：确保 spliced/unspliced assay 存在
  has_loom <- all(c("spliced", "unspliced") %in% Assays(seu))
  if (!has_loom) {
    if (is.null(loom_qs_file))
      stop(glue(
        "❌ seu 中缺少 spliced/unspliced assay 且未提供 loom_qs_file\n",
        "   请提供 loom_qs_file（qs 文件 / loom 目录 / 单个 loom 文件），",
        "或先运行 merge_loom_files() + add_loom_to_seurat()"
      ))

    if (verbose) message("\n── Step 1：添加 spliced/unspliced assay ──")

    ## 解析 loom_qs_file 的三种输入格式
    resolved_qs <- .resolve_loom_input(
      loom_input = loom_qs_file,
      cache_dir  = output_dir,
      verbose    = verbose
    )

    seu <- add_loom_to_seurat(seu, loom_qs_file = resolved_qs, verbose = verbose)

  } else {
    if (verbose) message("\n── Step 1：spliced/unspliced assay 已存在，跳过 ──")
  }

  ## Step 2：提取矩阵
  if (verbose) message("\n── Step 2：提取 spliced/unspliced 矩阵 ──")
  emat <- GetAssayData(seu, assay = "spliced",   layer = "counts")
  nmat <- GetAssayData(seu, assay = "unspliced", layer = "counts")
  n_groups <- length(unique(seu[[cluster_col, drop = TRUE]]))
  if (verbose) {
    message(glue("  spliced:   {nrow(emat)} genes × {ncol(emat)} cells"))
    message(glue("  unspliced: {nrow(nmat)} genes × {ncol(nmat)} cells"))
    message(glue("  cluster_col: '{cluster_col}' ({n_groups} groups)"))
  }

  ## Step 3：过滤低表达基因
  if (verbose) message("\n── Step 3：过滤低表达基因 ──")
  clusters  <- seu[[cluster_col, drop = TRUE]]

  if (verbose) message(glue("  [spliced]   threshold = {spliced_min_avg}"))
  emat_filt <- filter_velocity_genes(emat, clusters,
                                      min_max_avg = spliced_min_avg,
                                      verbose = verbose)
  if (verbose) message(glue("  [unspliced] threshold = {unspliced_min_avg}"))
  nmat_filt <- filter_velocity_genes(nmat, clusters,
                                      min_max_avg = unspliced_min_avg,
                                      verbose = verbose)

  ## Step 4：共同基因集
  if (verbose) message("\n── Step 4：取共同基因集 ──")
  genes_use <- intersect(rownames(emat_filt), rownames(nmat_filt))
  if (length(genes_use) == 0L) stop("❌ 过滤后没有共同基因，请降低阈值")
  emat_filt <- emat_filt[genes_use, ]
  nmat_filt <- nmat_filt[genes_use, ]
  if (verbose) {
    message(glue("  共同基因数：{length(genes_use)}"))
    message(glue("  最终矩阵：{length(genes_use)} genes × {ncol(emat_filt)} cells"))
  }

  ## Step 5：构建导出用 Seurat 对象
  if (verbose) message("\n── Step 5：构建导出用 Seurat 对象（v4 格式）──")
  seu2 <- build_velocity_seurat(
    seu                = seu,
    emat_filt          = emat_filt,
    nmat_filt          = nmat_filt,
    reductions_to_keep = reductions_to_keep,
    verbose            = verbose
  )

  ## Step 6：导出为 h5ad
  if (verbose) message("\n── Step 6：导出为 h5ad ──")
  h5ad_path <- export_velocity_h5ad(
    seu2          = seu2,
    output_dir    = output_dir,
    prefix        = prefix,
    keep_h5seurat = keep_h5seurat,
    verbose       = verbose
  )

  if (verbose) {
    message("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    message("  ✅ 全部完成！")
    message(glue("  h5ad：{h5ad_path}"))
    message(glue("  {ncol(seu2)} cells × {nrow(seu2)} genes（过滤后）"))
    message("  下一步：运行 scVelo Python 脚本")
    message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  }

  return(invisible(h5ad_path))
}


## ============================================================ ##
## 内部辅助函数：解析 loom_qs_file 的三种输入格式
## 不对外暴露（以 . 开头）
## ============================================================ ##

#' @keywords internal
.resolve_loom_input <- function(loom_input, cache_dir = ".", verbose = TRUE) {

  loom_input <- path.expand(loom_input)

  ## ── 情况 1：合并好的 qs 文件 ────────────────────────────────
  if (grepl("\\.qs$", loom_input, ignore.case = TRUE)) {
    if (!file.exists(loom_input))
      stop(glue("❌ qs 文件不存在: {loom_input}"))
    if (verbose) message(glue("  📦 检测到 qs 文件，直接使用: {loom_input}"))
    return(loom_input)
  }

  ## ── 情况 2：原始 loom 文件目录 ──────────────────────────────
  if (dir.exists(loom_input)) {
    loom_files <- list.files(loom_input, pattern = "\\.loom$",
                              full.names = FALSE, ignore.case = TRUE)
    if (length(loom_files) == 0L)
      stop(glue("❌ 目录 {loom_input} 中未找到 .loom 文件"))

    if (verbose) {
      message(glue("  📂 检测到 loom 目录，发现 {length(loom_files)} 个文件"))
      message(glue("  ⏳ 自动合并 loom 文件..."))
    }

    ## 合并结果缓存到 output_dir，避免重复合并
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    cache_file <- file.path(cache_dir, "loom_merged_cache.qs")

    merge_loom_files(
      loom_dir    = loom_input,
      output_file = cache_file,
      verbose     = verbose
    )

    if (verbose) message(glue("  💾 合并缓存已保存: {cache_file}"))
    return(cache_file)
  }

  ## ── 情况 3：单个 loom 文件 ──────────────────────────────────
  if (grepl("\\.loom$", loom_input, ignore.case = TRUE)) {
    if (!file.exists(loom_input))
      stop(glue("❌ loom 文件不存在: {loom_input}"))

    if (verbose)
      message(glue("  🗂️  检测到单个 loom 文件，读取并转为 qs 缓存..."))

    ## 读取单个 loom 并保存为 qs
    ldat <- read.loom.matrices(loom_input)

    ## cell ID 修正（与 merge_loom_files 默认逻辑一致）
    fix_id <- function(x) {
      x <- gsub("_", "-", x)
      x <- sub(":", "_",  x)
      x <- sub("x$", "-1", x)
      x
    }
    for (i in seq_along(ldat))
      colnames(ldat[[i]]) <- fix_id(colnames(ldat[[i]]))

    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    cache_file <- file.path(cache_dir, "loom_merged_cache.qs")
    qs::qsave(ldat, cache_file)

    if (verbose) {
      mb <- round(file.size(cache_file) / 1024 / 1024, 1)
      message(glue("  💾 单个 loom 转换完成: {cache_file}（{mb} MB）"))
    }
    return(cache_file)
  }

  ## ── 无法识别的格式 ──────────────────────────────────────────
  stop(glue(
    "❌ 无法识别 loom_qs_file: {loom_input}\n",
    "   支持格式：\n",
    "   ① 合并好的 qs 文件（*.qs）\n",
    "   ② 包含 loom 文件的目录\n",
    "   ③ 单个 loom 文件（*.loom）"
  ))
}

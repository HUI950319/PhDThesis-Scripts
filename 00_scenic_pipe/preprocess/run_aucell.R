#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)
library(tidyverse)
library(argparse)

# 创建参数解析器
parser <- ArgumentParser(description = "Calculate gene module scores using AUCell or UCell")

# 必须参数
parser$add_argument("--input_file", required = TRUE, 
                   help = "Input Seurat object file (*.rds or *.qs)")
parser$add_argument("--regulon_file", required = TRUE,
                   help = "Input regulon file (*.gmt)")
parser$add_argument("--output_path", required = TRUE,
                   help = "Output directory path for results")

# 可选参数（带默认值）
parser$add_argument("--method", choices = c("AUCell", "UCell"), default = "AUCell",
                   help = "Method for calculating signature score (default: AUCell)")
parser$add_argument("--min_size", type = "integer", default = 10,
                   help = "Minimum number of genes in gene sets (default: 10)")
parser$add_argument("--batch_size", type = "integer", default = 500,
                   help = "Number of cells per batch for AUCell (default: 500)")
parser$add_argument("--cores", type = "integer", default = 10,
                   help = "Number of threads for parallel computing (default: 10)")
parser$add_argument("--assay", default = NULL,
                   help = "Name of the Seurat object assay (default: DefaultAssay)")

# 解析参数
args <- parser$parse_args()

# 参数验证函数
validate_parameters <- function(args) {
  # 检查输入文件是否存在
  if (!file.exists(args$input_file)) {
    stop("Error: Input file '", args$input_file, "' does not exist!")
  }
  
  if (!file.exists(args$regulon_file)) {
    stop("Error: Regulon file '", args$regulon_file, "' does not exist!")
  }
  
  # 检查输入文件扩展名
  if (!grepl("\\.(rds|qs)$", args$input_file, ignore.case = TRUE)) {
    stop("Error: Input file must have .rds or .qs extension!")
  }
  
  if (!grepl("\\.gmt$", args$regulon_file, ignore.case = TRUE)) {
    stop("Error: Regulon file must have .gmt extension!")
  }
  
  # 检查输出目录是否存在，如果不存在则创建
  if (!dir.exists(args$output_path)) {
    cat("Creating output directory:", args$output_path, "\n")
    dir.create(args$output_path, recursive = TRUE)
  }
  
  # 检查数值参数的有效性
  if (args$min_size <= 0) {
    stop("Error: min_size must be positive!")
  }
  
  if (args$batch_size <= 0) {
    stop("Error: batch_size must be positive!")
  }
  
  if (args$cores <= 0) {
    stop("Error: cores must be positive!")
  }
  
  cat("Parameter validation passed!\n")
}

# 执行参数验证
validate_parameters(args)

# 打印参数信息
cat("Running gene module score calculation with parameters:\n")
cat("Input file:", args$input_file, "\n")
cat("Regulon file:", args$regulon_file, "\n")
cat("Output path:", args$output_path, "\n")
cat("Method:", args$method, "\n")
cat("Min size:", args$min_size, "\n")
cat("Batch size:", args$batch_size, "\n")
cat("Cores:", args$cores, "\n")
if (!is.null(args$assay)) {
  cat("Assay:", args$assay, "\n")
}
cat("---\n")

#' Calculate gene module scores
#'
#' @param x gene expression matrix, rows are genes, columns are cells. Can be any
#' format, UMI, CPM, TPM, etc.
#' @param ... Arguments passed to other methods.
#' @return A signature score matrix or Seurat object.
#' @concept compute_module_score
#' @export
ComputeModuleScore <- function(x, ...) UseMethod('ComputeModuleScore')


#' @param gene.sets a list of gene sets in data.frame or named list.
#' @param bg.genes background genes for calculating signature score, NULL means use
#' all genes. Default: NULL
#' @param method method for calculating signature score, support "AUCell" or "UCell".
#' Default: UCell
#' @param min.size The minimal genes of the gene sets. The size of gene sets less
#' than this value were ignored.  Default: 20
#' @param batch.size The number of cells were calculated for each batch.
#' This parameter is for the parallel calculation in adopted in 'AUCell' method.
#' Default: 500
#' @param cores number of threads for parallel computing. Default: 1
#'
#' @rdname ComputeModuleScore
#' @concept compute_module_score
#' @export
ComputeModuleScore.default <- function(x, gene.sets, bg.genes=NULL, method="UCell",
                                       min.size=20, batch.size=500, cores=1, ...) {
  # Check if UCell and AUCell packages are installed
  if (method == "AUCell" && !requireNamespace("AUCell", quietly = TRUE)) {
    stop("AUCell package is not installed. Please install the package using BiocManager::install('AUCell').")
  } else if (method == "UCell" && !requireNamespace("UCell", quietly = TRUE)) {
    stop("UCell package is not installed. Please install the package using remotes::install_github('carmonalab/UCell', ref='v1.3').")
  }
  if (!is.list(gene.sets)) {
    stop("'gene.sets' should be a list or data.frame!")
  }
  # filter the gene sets
  gene.sets <- gene.sets[sapply(gene.sets, length) >= min.size]
  # check the bg.genes
  if (!is.null(bg.genes)) {
    genes.use <- intersect(bg.genes, rownames(x))
    matched.ratio <- length(genes.use) / length(bg.genes)
    matched.ratio <- round(matched.ratio, 4)*100
    if (matched.ratio < 50) {
      stop(sprintf("Too less background genes (%s%%) were matched between query and reference. Please check the gene names in 'x'", matched.ratio))
    }
    if (matched.ratio < 70) {
      warning(sprintf("Only %s%% background genes were matched between query and reference", matched.ratio))
    } else {
      message(sprintf("%s%% background genes were matched between query and reference", matched.ratio))
    }
    x <- x[genes.use, ]
  }
  # calculate signature scores
  if (method == "AUCell") {
    n.cells <- ncol(x)
    batches <- floor((1:n.cells-1) / batch.size)
    batch.levels <- unique(batches)
    aucell <- function(i) {
      dge.tmp <- x[, batches == i]
      cr <- AUCell::AUCell_buildRankings(dge.tmp, nCores=1, plotStats=F, verbose = F)
      auc <- AUCell::AUCell_calcAUC(gene.sets, cr, nCores=1, verbose = F)
      AUCell::getAUC(auc)
    }
    auc_scores <- parallel::mclapply(batch.levels, aucell, mc.cores = cores)
    auc_scores <- do.call(cbind, auc_scores)
    auc_scores <- as.data.frame(t(auc_scores))
  } else if (method == "UCell") {
    auc_scores <- UCell::ScoreSignatures_UCell(x, features = gene.sets, ncores = cores)
  } else {
    stop(sprintf("%s is not supported! Please set method to 'AUCell' or 'UCell'.", method))
  }
  colnames(auc_scores) <- names(gene.sets)
  return(auc_scores)
}


#' @rdname ComputeModuleScore
#' @param assay Name of the seurat object assay.
#' @concept compute_module_score
#' @export
ComputeModuleScore.Seurat <- function(x, gene.sets, bg.genes=NULL, method="UCell",
                                      min.size=20, batch.size=500, cores=1,
                                      assay = Seurat::DefaultAssay(x), ...) {
  dge <- Seurat::GetAssayData(x, slot = "counts", assay = assay)
  if (ncol(dge) != ncol(x)) {
    dge <- Seurat::GetAssayData(x, slot = "data", assay = assay)
    message("Using 'data' slot for computing gene module score.")
  } else {
    message("Using 'counts' slot for computing gene module score.")
  }
  auc_scores <- ComputeModuleScore.default(x = dge, gene.sets, bg.genes, method,
                                           min.size, batch.size, cores)
  x[[method]] <- Seurat::CreateAssayObject(data = t(auc_scores))
  return(x)
}

#### 读入数据 ####
cat("Reading Seurat object...\n")
if (grepl("\\.qs$", args$input_file, ignore.case = TRUE)) {
  seu <- qs::qread(args$input_file)
} else if (grepl("\\.rds$", args$input_file, ignore.case = TRUE)) {
  seu <- readRDS(args$input_file)
} else {
  stop("Error: Input file must have .rds or .qs extension!")
}

cat("Reading regulon file...\n")
regulons <- clusterProfiler::read.gmt(args$regulon_file)

## data.frame -> list, list中的每个元素为一个gene set
rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)

cat("Regulon summary:\n")
cat("Number of regulons:", length(regulon.list), "\n")
cat("Regulon size summary:\n")
print(summary(sapply(regulon.list, length)))
cat("First regulon example:\n")
print(regulon.list[1])

# 保存regulon list
regulon_list_file <- file.path(args$output_path, "regulons.rds")
cat("Saving regulon list to:", regulon_list_file, "\n")
saveRDS(regulon.list, regulon_list_file)

#### 计算模块得分 ####
cat("Calculating gene module scores using", args$method, "method...\n")
seu <- ComputeModuleScore(seu, gene.sets = regulon.list, 
                         min.size = args$min_size, 
                         method = args$method, 
                         cores = args$cores,
                         batch.size = args$batch_size)

# 提取得分矩阵
gs.mat <- seu[[args$method]]@data %>% t()

# 保存结果
regulon_activity_file <- file.path(args$output_path, "regulon_activity.qs")
cat("Saving regulon activity matrix to:", regulon_activity_file, "\n")
qs::qsave(gs.mat, regulon_activity_file)

cat("Gene module score calculation completed successfully!\n")

#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)
library(tidyverse)
library(argparse)

# 创建参数解析器
parser <- ArgumentParser(description = "SCENIC preprocessing script for gene expression data")

# 必须参数
parser$add_argument("--input_file", required = TRUE, 
                   help = "Input Seurat object file (*.rds or *.qs)")
parser$add_argument("--output_path", required = TRUE,
                   help = "Output directory path for results")
parser$add_argument("--species", required = TRUE, choices = c("human", "mouse"),
                   help = "Species (human or mouse)")

# 可选参数（带默认值）
parser$add_argument("--n_cells", type = "integer", default = 2000,
                   help = "Number of cells to sketch (default: 2000)")
parser$add_argument("--threads", type = "integer", default = 5,
                   help = "Number of threads to use (default: 5)")
parser$add_argument("--k", type = "integer", default = 100,
                   help = "Number of components to use (default: 100)")
parser$add_argument("--min_cells_per_gene", type = "integer", default = 5,
                   help = "Minimum number of cells to keep a gene (default: 5)")

# 解析参数
args <- parser$parse_args()

# 根据物种获取cisTarget文件路径
get_cistarget_file <- function(species) {
  if (species == "human") {
    return("cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
  } else if (species == "mouse") {
    return("cistarget/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
  } else {
    stop("Error: Unsupported species '", species, "'. Supported: human, mouse")
  }
}

# 参数验证函数
validate_parameters <- function(args) {
  # 检查输入文件是否存在
  if (!file.exists(args$input_file)) {
    stop("Error: Input file '", args$input_file, "' does not exist!")
  }
  
  # 检查输入文件扩展名
  if (!grepl("\\.(rds|qs)$", args$input_file, ignore.case = TRUE)) {
    stop("Error: Input file must have .rds or .qs extension!")
  }
  
  # 获取并检查cisTarget数据库文件
  cisdb_file <- get_cistarget_file(args$species)
  if (!file.exists(cisdb_file)) {
    stop("Error: cisTarget motif ranking file '", cisdb_file, "' does not exist!")
  }
  
  # 检查输出目录是否存在，如果不存在则创建
  if (!dir.exists(args$output_path)) {
    cat("Creating output directory:", args$output_path, "\n")
    dir.create(args$output_path, recursive = TRUE)
  }
  
  # 检查数值参数的有效性
  if (args$n_cells <= 0) {
    stop("Error: n_cells must be positive!")
  }
  
  if (args$threads <= 0) {
    stop("Error: threads must be positive!")
  }
  
  if (args$k <= 0) {
    stop("Error: k must be positive!")
  }
  
  if (args$min_cells_per_gene <= 0) {
    stop("Error: min_cells_per_gene must be positive!")
  }
  
  cat("Parameter validation passed!\n")
}

# 执行参数验证
validate_parameters(args)

# 获取cisTarget文件路径
cisdb_file <- get_cistarget_file(args$species)

# 打印参数信息
cat("Running SCENIC preprocessing with parameters:\n")
cat("Input file:", args$input_file, "\n")
cat("Output path:", args$output_path, "\n")
cat("Species:", args$species, "\n")
cat("cisTarget file:", cisdb_file, "\n")
cat("Number of subsampled cells:", args$n_cells, "\n")
cat("Threads:", args$threads, "\n")
cat("Components (k):", args$k, "\n")
cat("Min cells per gene:", args$min_cells_per_gene, "\n")
cat("---\n")

### 读入数据
cat("Reading Seurat object...\n")
if (grepl("\\.qs$", args$input_file, ignore.case = TRUE)) {
  seu <- qs::qread(args$input_file)
} else if (grepl("\\.rds$", args$input_file, ignore.case = TRUE)) {
  seu <- readRDS(args$input_file)
} else {
  stop("Error: Input file must have .rds or .qs extension!")
}

cat("Normalizing data...\n")
seu <- NormalizeData(seu)

## sketch cells
cat("Finding variable features...\n")
seu <- FindVariableFeatures(seu)

# 检测Seurat版本并选择合适的sketching方法
seurat_version <- packageVersion("Seurat")
cat("Seurat version:", as.character(seurat_version), "\n")

cat("Sampling cells...\n")
if (seurat_version >= "5.0.0") {
  # Seurat V5 使用LeverageScore方法
  cat("Using LeverageScore method for Seurat V5\n")
  seu <- SketchData(seu, ncells = args$n_cells, method = "LeverageScore")
  sketch.cells <- Cells(seu)
  cat("Extracting expression matrix...\n")
  expr.mat <- LayerData(seu, layer = "data", assay = "RNA")
} else {
  # Seurat V4 使用random sampling
  cat("Using random sampling for Seurat V4\n")
  set.seed(42)  # 设置随机种子以确保可重复性
  sketch.cells <- sample(Cells(seu), size = args$n_cells)
  expr.mat <- GetAssayData(seu, slot = "data", assay = "RNA")
}

expr.in.cells <- rowSums(expr.mat > 0)
expr.mat <- expr.mat[expr.in.cells >= args$min_cells_per_gene, ]

cat("Reading cisTarget database...\n")
cisdb <- arrow::read_feather(cisdb_file)
genes.use <- intersect(colnames(cisdb), rownames(expr.mat))
expr.mat <- expr.mat[genes.use, ]

cat("Running NMF with", args$threads, "threads...\n")
RcppML::setRcppMLthreads(args$threads)
# k << number of genes
# Run NMF multiple times with different seeds and average the results
model <- RcppML::nmf(expr.mat, k = args$k, verbose = T)

cat("Computing imputed matrix...\n")
imputed.mat <- model$w %*% diag(model$d) %*% model$h
colnames(imputed.mat) <- colnames(expr.mat)
rownames(imputed.mat) <- rownames(expr.mat)
imputed.mat = imputed.mat[, sketch.cells]
imputed.mat = round(t(imputed.mat), 2)

# 构建输出文件路径
output_file <- file.path(args$output_path, "imputed.mat.csv")

cat("Writing results to:", output_file, "\n")
write.csv(imputed.mat, output_file, row.names = T)

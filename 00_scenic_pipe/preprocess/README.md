# SCENIC Preprocessing Script


这个脚本用于SCENIC流程的预处理步骤，支持命令行参数传递。

## 安装依赖

确保安装了所需的R包：

```r
install.packages(c("argparse", "tidyverse"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Seurat", "Matrix"))
install.packages("RcppML")
```

## 使用方法

### 基本用法（只使用必须参数）

```bash
Rscript preprocess.R \
  --input_file "./scenic_pipe/input/ifnb_pbmc.seurat.rds" \
  --output_path "./scenic_pipe/output" \
  --cisdb_motif_ranking_file "cisTarget_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
```

### 完整用法（包含所有可选参数）

```bash
Rscript preprocess.R \
  --input_file "./scenic_pipe/input/ifnb_pbmc.seurat.rds" \
  --output_path "./scenic_pipe/output" \
  --cisdb_motif_ranking_file "cisTarget_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  --n_cells 3000 \
  --threads 8 \
  --k 150 \
  --min_cells_per_gene 10
```

## 参数说明

### 必须参数

- `--input_file`: 输入的Seurat对象文件路径（支持.rds或.qs格式）
- `--output_path`: 输出目录路径
- `--cisdb_motif_ranking_file`: cisTarget数据库motif排序文件路径

### 可选参数

- `--n_cells`: 用于sketch的细胞数量（默认：2000）
- `--threads`: 使用的线程数（默认：5）
- `--k`: NMF使用的组件数（默认：100）
- `--min_cells_per_gene`: 保留基因的最小细胞数（默认：5）

## 参数验证

脚本会自动验证以下内容：

1. 输入文件是否存在
2. 输入文件格式是否正确（.rds或.qs）
3. cisTarget数据库文件是否存在
4. 输出目录是否存在（如果不存在会自动创建）
5. 所有数值参数是否为正数

## 输出文件

脚本会在指定的输出目录中生成：
- `imputed.mat.csv`: 插补后的基因表达矩阵

## 错误处理

如果参数验证失败，脚本会显示详细的错误信息并停止执行。常见的错误包括：

- 文件不存在
- 文件格式不正确
- 参数值无效（如负数）

## 示例输出

```
Parameter validation passed!
Running SCENIC preprocessing with parameters:
Input file: ./scenic_pipe/input/ifnb_pbmc.seurat.rds
Output path: ./scenic_pipe/output
cisTarget file: cisTarget_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
Number of cells: 2000
Threads: 5
Components (k): 100
Min cells per gene: 5
---
Reading Seurat object...
Normalizing data...
Finding variable features...
Sketching cells...
Extracting expression matrix...
Reading cisTarget database...
Running NMF with 5 threads...
Computing imputed matrix...
Writing results to: ./scenic_pipe/output/imputed.mat.csv
Preprocessing completed successfully!
``` 
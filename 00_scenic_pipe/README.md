# SCENIC Pipeline — 基因调控网络 (GRN) 推断框架

本 pipeline 实现了完整的 SCENIC (Single-Cell Regulatory Network Inference and Clustering) 分析流程，用于从单细胞 RNA-seq 数据中推断基因调控网络 (Gene Regulatory Network, GRN)。

## 目录

- [流程概览](#流程概览)
- [文件结构](#文件结构)
- [环境配置](#环境配置)
- [配置文件说明](#配置文件说明)
- [使用方法](#使用方法)
- [分步详解](#分步详解)
- [输出文件](#输出文件)
- [常见问题](#常见问题)

---

## 流程概览

SCENIC 分析分为三个核心步骤：

`
Step 1: 数据预处理 (preprocess.R)
    Seurat 对象 → NormalizeData → SketchData → NMF 填补 → imputed.mat.csv

Step 2: pySCENIC GRN 推断 (run_pyscenic.py)
    imputed.mat.csv → GRNBoost2 共表达网络 → cisTarget motif 富集 → regulons

Step 3: AUCell 活性评分 (run_aucell.R)
    Seurat 对象 + regulons.gmt → AUCell/UCell 评分 → regulon_activity.qs
`

## 文件结构

`
scenic_pipe/
├── main.py              # 主控制脚本（argparse 调度三步流程）
├── run_scenic.R         # R 包装器（自动激活 conda 环境并调用 main.py）
├── configs.txt          # 默认配置模板
├── configs-seu.txt      # 甲状旁腺项目配置（使用示例）
├── environment.yml      # Conda 环境定义文件（pyscenic-env, Python 3.10）
├── README.md            # 本文档
├── preprocess/
│   ├── preprocess.R     # Step 1: R 预处理脚本
│   ├── run_pyscenic.py  # Step 2: Python pySCENIC 脚本
│   ├── run_aucell.R     # Step 3: AUCell 评分脚本
│   └── README.md        # 预处理模块说明
└── cistarget/           # cisTarget 数据库文件
    ├── hg38_*.feather   # Human (hg38) motif ranking 数据库
    ├── mm10_*.feather   # Mouse (mm10) motif ranking 数据库
    ├── *.tbl            # motif 注释表
    └── *.txt            # TF 列表文件
`

## 环境配置

### 1. R 环境

需要安装以下 R 包：

`
install.packages(c("Seurat", "Matrix", "tidyverse", "argparse", "qs", "arrow", "RcppML"))

# Bioconductor 包
BiocManager::install(c("AUCell", "clusterProfiler"))

# 可选：UCell（如需使用 UCell 方法）
remotes::install_github("carmonalab/UCell", ref = "v1.3")
`

### 2. Python 环境（pySCENIC）

#### 方法一：使用 Conda 环境文件（推荐）

`ash
# 使用提供的 environment.yml 创建环境
conda env create -f environment.yml

# 激活环境
conda activate pyscenic-env

# 验证安装
python -c "import pyscenic; print('pySCENIC version:', pyscenic.__version__)"
# 期望输出: pySCENIC version: 0.12.1
`

> **注意**：如果下载速度慢，可使用清华镜像源：
> `ash
> conda env create -f environment.yml -c conda-forge -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
> `

#### 方法二：手动安装

`ash
conda create -n pyscenic-env python=3.10
conda activate pyscenic-env
conda install -c conda-forge pyscenic
pip install pandas configparser
`

### 3. cisTarget 数据库

pipeline 已在 cistarget/ 目录中提供了必要的数据库文件，无需额外下载。支持的物种：

| 物种 | 基因组 | 数据库文件 |
|------|--------|-----------|
| Human | hg38 | hg38_*_10kbp_up_10kbp_down_*.feather |
| Mouse | mm10 | mm10_*_10kbp_up_10kbp_down_*.feather |

pipeline 会根据配置文件中的 species 参数自动选择对应的数据库文件。

### 4. 验证环境

`ash
cd scenic_pipe
conda activate pyscenic-env
python main.py --check-deps
`

---

## 配置文件说明

配置文件采用 INI 格式，包含四个部分。以 configs-seu.txt（甲状旁腺项目配置）为例：

`ini
[global]
species = human                # 物种：human 或 mouse
output_path = output/seu       # 输出目录路径

[preprocess]
rscript_path = Rscript         # Rscript 路径（或绝对路径如 /opt/R/4.4.2/bin/Rscript）
input_file = input/seu.qs      # 输入 Seurat 对象文件（.rds 或 .qs）
n_cells = 6000                 # SketchData 采样细胞数
threads = 10                   # 并行线程数
k = 100                        # NMF 因子数
min_cells_per_gene = 10        # 基因最少表达细胞数

[pyscenic]
num_workers = 15               # GRNBoost2 并行 worker 数
seed = 777                     # 随机种子
min_regulon_size = 10          # regulon 最小基因数

[aucell]
method = AUCell                # 评分方法：AUCell 或 UCell
min_size = 10                  # 基因集最小基因数
batch_size = 500               # AUCell 每批处理细胞数
cores = 10                     # 并行核数
`

### 参数调优建议

- **n_cells**: 建议 2000–6000，数据量大时取较大值以保留更多信息
- **k (NMF 因子数)**: 默认 100，细胞类型多时可适当增大
- **num_workers**: 根据 CPU 核心数设置，建议不超过物理核心数
- **batch_size**: 内存紧张时可减小此值

---

## 使用方法

### 方式一：通过 Python 主脚本运行（推荐）

`ash
# 激活 conda 环境
conda activate pyscenic-env

# 进入 pipeline 目录
cd scenic_pipe

# 使用默认配置运行完整流程
python main.py

# 使用自定义配置文件
python main.py --config configs-seu.txt

# 跳过已完成的步骤
python main.py --config configs-seu.txt --skip-preprocess    # 跳过预处理
python main.py --config configs-seu.txt --skip-pyscenic      # 跳过 pySCENIC
python main.py --config configs-seu.txt --skip-aucell        # 跳过 AUCell

# 遇到错误时继续执行后续步骤
python main.py --config configs-seu.txt --continue-on-error
`

### 方式二：通过 R 包装器运行

un_scenic.R 提供了 R 接口，自动检测并激活 conda 环境：

`
# 在 R 中 source 后调用
source("run_scenic.R")

# 使用默认配置运行
result <- run_scenic(config_file = "configs-seu.txt")

# 跳过预处理步骤
result <- run_scenic(config_file = "configs-seu.txt", skip_preprocess = TRUE)

# 指定 conda 路径和环境名称
result <- run_scenic(
  config_file = "configs-seu.txt",
  conda_path = "~/miniconda3",
  env_name = "pyscenic-env"
)
`

或通过命令行调用：

`ash
Rscript run_scenic.R configs-seu.txt
`

### 方式三：分步单独运行

也可以逐步手动执行各脚本：

`ash
conda activate pyscenic-env

# Step 1: 预处理
Rscript preprocess/preprocess.R \
  --input_file input/seu.qs \
  --output_path output/seu \
  --species human \
  --n_cells 6000 --threads 10 --k 100 --min_cells_per_gene 10

# Step 2: pySCENIC GRN 推断
python preprocess/run_pyscenic.py \
  --input_csv output/seu/imputed.mat.csv \
  --output_path output/seu \
  --species human \
  --num_workers 15 --seed 777

# Step 3: AUCell 活性评分
Rscript preprocess/run_aucell.R \
  --input_file input/seu.qs \
  --regulon_file output/seu/regulons.gmt \
  --output_path output/seu \
  --method AUCell --min_size 10 --batch_size 500 --cores 10
`

---

## 分步详解

### Step 1: 数据预处理 (preprocess/preprocess.R)

**功能**：从 Seurat 对象提取表达矩阵，经过标准化、特征选择、细胞抽样和 NMF 填补，生成适用于 GRNBoost2 的稠密表达矩阵。

**处理流程**：
1. 读取 Seurat 对象（.rds 或 .qs 格式）
2. NormalizeData() 标准化
3. FindVariableFeatures() 筛选高变基因
4. SketchData() 抽样 n_cells 个细胞（大数据集降采样）
5. 提取表达矩阵，过滤低表达基因（min_cells_per_gene）
6. 使用 RcppML::nmf() 进行非负矩阵分解（NMF）
7. 计算 NMF 填补矩阵 (W × H)，输出为 imputed.mat.csv

**输入**：Seurat 对象（需包含 RNA assay 和原始 counts）
**输出**：imputed.mat.csv（行=细胞，列=基因的稠密矩阵）

> **为什么使用 NMF 填补？**
> 单细胞数据存在大量 dropout（零值），直接使用稀疏矩阵会导致 GRNBoost2 推断的共表达网络质量低下。NMF 填补能恢复被 dropout 掩盖的真实表达信号，显著提升 GRN 推断准确性。

### Step 2: pySCENIC GRN 推断 (preprocess/run_pyscenic.py)

**功能**：基于填补后的表达矩阵，依次执行 GRNBoost2 共表达分析和 cisTarget motif 富集验证，输出经过验证的 regulons。

**处理流程**：
1. 读取 imputed.mat.csv
2. **GRNBoost2**（rboreto_with_multiprocessing.py）：基于梯度提升树算法推断 TF-target 共表达关系 → grn.tsv
3. **cisTarget**（pyscenic ctx）：使用 cisTarget 数据库验证 TF binding motif 富集，过滤假阳性 → ctx.tsv
4. 从 ctx 结果中提取 regulon 基因集 → egulons.gmt + egulons.txt

**输入**：imputed.mat.csv（来自 Step 1）
**输出**：
- grn.tsv：原始 TF-target 共表达网络
- ctx.tsv：cisTarget motif 验证后的 regulons
- egulons.gmt：GMT 格式的 regulon 基因集（用于 Step 3）
- egulons.txt：可读文本格式的 regulons

> **关于运行时间**：GRNBoost2 是计算密集型步骤，6000 个细胞的数据集大约需要 30–60 分钟（取决于 num_workers 和 CPU 性能）。

### Step 3: AUCell 活性评分 (preprocess/run_aucell.R)

**功能**：使用 AUCell 或 UCell 方法计算每个细胞在每个 regulon 上的活性得分，生成 regulon 活性矩阵。

**处理流程**：
1. 读取原始 Seurat 对象和 regulons.gmt 文件
2. 将 GMT 格式转换为基因集列表
3. 使用 AUCell 方法：
   - AUCell_buildRankings()：对每个细胞的基因按表达量排序
   - AUCell_calcAUC()：计算各 regulon 的 AUC 得分
4. 分批并行计算（batch_size 控制每批细胞数）
5. 保存 regulon 列表和活性矩阵

**输入**：
- Seurat 对象（原始 .rds/.qs 文件，同 Step 1 输入）
- egulons.gmt（来自 Step 2）

**输出**：
- egulons.rds：R list 格式的 regulon 基因列表
- egulon_activity.qs：**核心结果**——细胞×regulon 活性得分矩阵

---

## 输出文件

所有输出文件保存在配置文件中指定的 output_path 目录下：

| 步骤 | 文件 | 说明 |
|------|------|------|
| Step 1 | imputed.mat.csv | NMF 填补后的表达矩阵（行=细胞，列=基因） |
| Step 2 | grn.tsv | GRNBoost2 原始共表达网络（TF-target 关系及重要性得分） |
| Step 2 | ctx.tsv | cisTarget motif 验证后的 regulons |
| Step 2 | egulons.gmt | GMT 格式 regulon 基因集（用于 AUCell 输入） |
| Step 2 | egulons.txt | 可读文本格式 regulons（逗号分隔基因列表） |
| Step 3 | egulons.rds | R list 格式的 regulon 基因列表 |
| Step 3 | egulon_activity.qs | **核心结果**：细胞×regulon 活性得分矩阵 |

### 下游分析示例

`
# 读取 regulon 活性矩阵
ras.mat <- qs::qread("output/seu/regulon_activity.qs")
dim(ras.mat)  # 行=细胞, 列=regulon

# 查看所有 regulon 名称
colnames(ras.mat)

# 读取 regulon 基因列表
regulons <- readRDS("output/seu/regulons.rds")
names(regulons)                  # 查看所有 regulon
regulons[["STAT1(+)"]]          # 查看 STAT1 regulon 的靶基因

# 将活性矩阵添加到 Seurat 对象
seu[["SCENIC"]] <- CreateAssayObject(data = t(ras.mat))

# 可视化特定 regulon 的活性
FeaturePlot(seu, features = "STAT1(+)", reduction = "umap")
`

---

## 常见问题

### 1. Rscript 找不到

`ash
# 检查 Rscript 是否在 PATH 中
which Rscript

# 如果不在 PATH 中，在配置文件中指定绝对路径
rscript_path = /opt/R/4.4.2/bin/Rscript
`

### 2. pySCENIC 安装失败

`ash
# 删除旧环境重新创建
conda env remove -n pyscenic-env
conda env create -f environment.yml

# 验证安装
conda activate pyscenic-env
python -c "import pyscenic; print(pyscenic.__version__)"
`

### 3. 内存不足 (OOM)

- 减小 
_cells 参数（如从 6000 降至 3000）
- 减小 
um_workers 和 	hreads
- 增加 atch_size 不会减少内存，减小它才会

### 4. imputed.mat.csv 已存在时的行为

main.py 会自动检测 imputed.mat.csv 是否存在：如果已存在，将跳过 Step 1 直接进入 Step 2。也可手动通过 --skip-preprocess 跳过。

### 5. 如何选择 AUCell 还是 UCell

- **AUCell**（默认）：经典方法，基于 AUC 曲线下面积，适合大多数场景
- **UCell**：Mann-Whitney U 统计量方法，在小基因集上可能更稳健

在配置文件的 [aucell] 部分设置 method = AUCell 或 method = UCell。

### 6. 支持哪些输入格式

Seurat 对象支持 .rds（eadRDS）和 .qs（qs::qread）两种格式。推荐使用 .qs 格式，读写速度更快且文件体积更小。

---

## Seurat 输入要求

| 必需项 | 说明 |
|--------|------|
| RNA assay | 包含原始 counts 数据 |
| 基因命名 | Human: HGNC 符号 (如 TP53)；Mouse: MGI 符号 (如 Trp53) |

SCENIC 流程本身**不需要**预先的降维或聚类结果。但如果后续需要可视化 regulon 活性的空间分布，建议输入对象包含 UMAP 坐标和细胞类型注释。
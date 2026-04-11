#!/bin/bash
# ============================================================
# 03a_build_mmul_ref.sh
# 构建恒河猴 (Macaca mulatta) CellRanger 参考基因组
# 基于 Ensembl Mmul_10 (rheMac10) + Ensembl annotation
# 需在运行 03_GSE233962_cellranger.sh 前完成
# ============================================================

set -euo pipefail

WORK_DIR="/home/oyh/project/scRNA"
REF_DIR="${WORK_DIR}/resource/ref"
RMSK_DIR="${WORK_DIR}/resource/rmsk"
mkdir -p ${REF_DIR}/mmul10_build ${RMSK_DIR}
cd ${REF_DIR}/mmul10_build

# --- 1) 下载 Ensembl 基因组和注释 ---
ENSEMBL_VER="112"
GENOME_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VER}/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VER}/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.${ENSEMBL_VER}.gtf.gz"

echo "下载 Mmul_10 基因组..."
wget -c ${GENOME_URL}
echo "下载 GTF 注释..."
wget -c ${GTF_URL}

gunzip Macaca_mulatta.Mmul_10.dna.primary_assembly.fa.gz
gunzip Macaca_mulatta.Mmul_10.${ENSEMBL_VER}.gtf.gz

GENOME="Macaca_mulatta.Mmul_10.dna.primary_assembly.fa"
GTF="Macaca_mulatta.Mmul_10.${ENSEMBL_VER}.gtf"

# --- 2) 过滤GTF: 仅保留protein_coding等基因 (与10x官方mkref流程一致) ---
cellranger mkgtf ${GTF} ${GTF%.gtf}.filtered.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA \
    --attribute=gene_biotype:IG_C_gene \
    --attribute=gene_biotype:IG_D_gene \
    --attribute=gene_biotype:IG_J_gene \
    --attribute=gene_biotype:IG_V_gene \
    --attribute=gene_biotype:TR_C_gene \
    --attribute=gene_biotype:TR_D_gene \
    --attribute=gene_biotype:TR_J_gene \
    --attribute=gene_biotype:TR_V_gene

# --- 3) 构建 CellRanger 参考基因组 ---
echo "构建 CellRanger mkref..."
cellranger mkref \
    --genome=refdata-gex-Mmul10 \
    --fasta=${GENOME} \
    --genes=${GTF%.gtf}.filtered.gtf \
    --nthreads=10 \
    --memgb=60

# 移动到标准位置
mv refdata-gex-Mmul10 ${REF_DIR}/

# --- 4) 下载 RepeatMasker 注释 (用于velocyto) ---
echo "下载 rheMac10 RepeatMasker..."
wget -O ${RMSK_DIR}/mmul10_rmsk.gtf.gz \
    "https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.rmsk.gtf.gz"
gunzip ${RMSK_DIR}/mmul10_rmsk.gtf.gz

echo "恒河猴参考基因组构建完成: ${REF_DIR}/refdata-gex-Mmul10"

#!/bin/bash
# ============================================================
# 00_build_star_index.sh
# 构建STAR基因组索引（GRCh38 + GENCODE v44）
# 只需运行一次，后续所有bulk样本共用
# ============================================================

# --- 路径配置 ---
WORK_DIR="/home/oyh/project/bulk"
REF_DIR="${WORK_DIR}/reference"
THREADS=16

mkdir -p ${REF_DIR}

# --- 下载参考基因组和注释 ---
cd ${REF_DIR}

# GRCh38 primary assembly
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# GENCODE v44 GTF注释
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

# --- 构建STAR索引 ---
STAR --runMode genomeGenerate \
     --runThreadN ${THREADS} \
     --genomeDir ${REF_DIR}/star_index \
     --genomeFastaFiles ${REF_DIR}/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile ${REF_DIR}/gencode.v44.annotation.gtf \
     --sjdbOverhang 149

echo "STAR index built: ${REF_DIR}/star_index"

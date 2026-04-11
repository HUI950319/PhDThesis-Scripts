#!/bin/bash
# ============================================================
# 02_GSE190773_cellranger.sh
# 公开数据集 GSE190773 甲状旁腺scRNA-seq
# SRA已下载 → fasterq-dump → fastq重命名 → CellRanger count → velocyto
# ============================================================

set -euo pipefail

# --- 路径配置 ---
WORK_DIR="/home/oyh/project/scRNA"
cd ${WORK_DIR}

export PATH=${WORK_DIR}/cellranger_soft/cellranger-8.0.1:$PATH

REF="${WORK_DIR}/resource/ref/refdata-gex-GRCh38-2020-A"
GTF="${REF}/genes/genes.gtf"
RMSK="${WORK_DIR}/resource/rmsk/hg38_rmsk.gtf"
SRA_DIR="${WORK_DIR}/sra/GSE190773"
FASTQ_DIR="${WORK_DIR}/fastq_public/GSE190773"
# sample_map.tsv: SRR编号\t样本名
MAP_FILE="${SRA_DIR}/sample_map.tsv"
THREADS=10
MEM=120

mkdir -p ${FASTQ_DIR}

# --- 1) fasterq-dump + fastq重命名 ---
while IFS=$'\t' read -r SRR SNAME; do
    [[ "${SRR}" =~ ^# ]] && continue
    echo "========== fasterq-dump: ${SRR} → ${SNAME} =========="

    mkdir -p ${FASTQ_DIR}/${SNAME}
    fasterq-dump ${SRA_DIR}/${SRR}/${SRR}.sra \
        -O ${FASTQ_DIR}/${SNAME} \
        -e ${THREADS} --split-3

    cd ${FASTQ_DIR}/${SNAME}
    mv ${SRR}_1.fastq ${SNAME}_S1_L001_R1_001.fastq
    mv ${SRR}_2.fastq ${SNAME}_S1_L001_R2_001.fastq
    pigz -p ${THREADS} *.fastq
    cd ${WORK_DIR}
done < ${MAP_FILE}

# --- 2) CellRanger count + velocyto ---
while IFS=$'\t' read -r SRR SNAME; do
    [[ "${SRR}" =~ ^# ]] && continue
    echo "========== CellRanger: ${SNAME} =========="

    cellranger count \
        --id=${SNAME} \
        --transcriptome="${REF}" \
        --fastqs="${FASTQ_DIR}/${SNAME}/" \
        --create-bam=true \
        --nosecondary \
        --localcores=${THREADS} \
        --localmem=${MEM}

    echo "========== Velocyto: ${SNAME} =========="
    velocyto run10x -@ ${THREADS} -m ${RMSK} ${SNAME} ${GTF}
done < ${MAP_FILE}

echo "GSE190773 pipeline complete."

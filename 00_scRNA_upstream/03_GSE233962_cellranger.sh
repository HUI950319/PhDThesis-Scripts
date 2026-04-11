#!/bin/bash
# ============================================================
# 03_GSE233962_cellranger.sh
# 公开数据集 GSE233962 (Venkat et al. 2024, Genome Research)
# 4个人类样本(Y7,Y9,Y11,Y13) + 4个恒河猴样本(Mmul1-4)
# 人类 → GRCh38 | 恒河猴 → Macaca mulatta (Mmul10/Ensembl)
# SRA已下载 → fasterq-dump → fastq重命名 → CellRanger → velocyto
# ============================================================

set -euo pipefail

# --- 路径配置 ---
WORK_DIR="/home/oyh/project/scRNA"
cd ${WORK_DIR}

export PATH=${WORK_DIR}/cellranger_soft/cellranger-8.0.1:$PATH

# 人类参考基因组
REF_HUMAN="${WORK_DIR}/resource/ref/refdata-gex-GRCh38-2020-A"
GTF_HUMAN="${REF_HUMAN}/genes/genes.gtf"
RMSK_HUMAN="${WORK_DIR}/resource/rmsk/hg38_rmsk.gtf"

# 恒河猴参考基因组 (需先运行 03a_build_mmul_ref.sh 构建)
REF_MMUL="${WORK_DIR}/resource/ref/refdata-gex-Mmul10"
GTF_MMUL="${REF_MMUL}/genes/genes.gtf"
RMSK_MMUL="${WORK_DIR}/resource/rmsk/mmul10_rmsk.gtf"

SRA_DIR="${WORK_DIR}/sra/GSE233962"
FASTQ_DIR="${WORK_DIR}/fastq_public/GSE233962"
# sample_map.tsv: SRR编号\t样本名\t物种(human/mmul)
MAP_FILE="${SRA_DIR}/sample_map.tsv"
THREADS=10
MEM=120

mkdir -p ${FASTQ_DIR}

# --- 1) fasterq-dump + fastq重命名 ---
while IFS=$'\t' read -r SRR SNAME SPECIES; do
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

# --- 2) CellRanger count + velocyto (按物种选择参考基因组) ---
while IFS=$'\t' read -r SRR SNAME SPECIES; do
    [[ "${SRR}" =~ ^# ]] && continue

    if [[ "${SPECIES}" == "human" ]]; then
        REF=${REF_HUMAN}; GTF=${GTF_HUMAN}; RMSK=${RMSK_HUMAN}
    elif [[ "${SPECIES}" == "mmul" ]]; then
        REF=${REF_MMUL}; GTF=${GTF_MMUL}; RMSK=${RMSK_MMUL}
    else
        echo "未知物种: ${SPECIES}，跳过 ${SNAME}"; continue
    fi

    echo "========== CellRanger [${SPECIES}]: ${SNAME} =========="
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

echo "GSE233962 pipeline complete (Human + NHP)."

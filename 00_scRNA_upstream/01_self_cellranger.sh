#!/bin/bash
# ============================================================
# 01_self_cellranger.sh
# 自测15例甲状旁腺scRNA-seq上游处理
# CellRanger count + velocyto loom生成
# PH1-6 (PHPT), PT1-3 (Normal), SH1-6 (SHPT)
# ============================================================

set -euo pipefail

# --- 路径配置 ---
WORK_DIR="/home/oyh/project/scRNA"
cd ${WORK_DIR}

export PATH=${WORK_DIR}/cellranger_soft/cellranger-8.0.1:$PATH

REF="${WORK_DIR}/resource/ref/refdata-gex-GRCh38-2020-A"
GTF="${REF}/genes/genes.gtf"
RMSK="${WORK_DIR}/resource/rmsk/hg38_rmsk.gtf"
FASTQ_DIR="${WORK_DIR}/fastq_10x"
THREADS=10
MEM=120

# --- 样本列表 ---
SAMPLES=(PH1 PH2 PH3 PH4 PH5 PH6 PT1 PT2 PT3 SH1 SH2 SH3 SH4 SH5 SH6)

# --- CellRanger count + velocyto ---
for S in "${SAMPLES[@]}"; do
    echo "========== CellRanger: ${S} =========="
    cellranger count \
        --id=${S} \
        --transcriptome="${REF}" \
        --fastqs="${FASTQ_DIR}/${S}/" \
        --create-bam=true \
        --nosecondary \
        --localcores=${THREADS} \
        --localmem=${MEM}

    echo "========== Velocyto: ${S} =========="
    velocyto run10x -@ ${THREADS} -m ${RMSK} ${S} ${GTF}

    echo "========== Done: ${S} =========="
done

echo "All 15 samples complete."
echo "CellRanger outputs: ${WORK_DIR}/{PH,PT,SH}*/outs/"
echo "Loom files: ${WORK_DIR}/{PH,PT,SH}*/velocyto/*.loom"

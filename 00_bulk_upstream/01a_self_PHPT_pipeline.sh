#!/bin/bash
# ============================================================
# 01a_self_PHPT_pipeline.sh
# 自测原发bulk RNA-seq（PHPT×12）
# fastq → fastp QC → STAR比对 → featureCounts定量
# ============================================================

set -euo pipefail

# --- 路径配置 ---
WORK_DIR="/home/oyh/project/bulk"
RAW_DIR="${WORK_DIR}/raw_data/self_PHPT"
CLEAN_DIR="${WORK_DIR}/clean_data/self_PHPT"
ALIGN_DIR="${WORK_DIR}/aligned/self_PHPT"
COUNTS_DIR="${WORK_DIR}/counts/self_PHPT"
QC_DIR="${WORK_DIR}/qc/self_PHPT"
STAR_INDEX="${WORK_DIR}/reference/star_index"
GTF="${WORK_DIR}/reference/gencode.v44.annotation.gtf"
THREADS=16

mkdir -p ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR} ${QC_DIR}

# --- 样本列表（根据实际文件名修改） ---
SAMPLES=(
    PHPT_01 PHPT_02 PHPT_03 PHPT_04 PHPT_05 PHPT_06
    PHPT_07 PHPT_08 PHPT_09 PHPT_10 PHPT_11 PHPT_12
)

# --- 循环处理 ---
for SAMPLE in "${SAMPLES[@]}"; do
    echo "========== Processing: ${SAMPLE} =========="
    R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"

    # fastp质控
    fastp \
        -i ${R1} -I ${R2} \
        -o ${CLEAN_DIR}/${SAMPLE}_R1_clean.fastq.gz \
        -O ${CLEAN_DIR}/${SAMPLE}_R2_clean.fastq.gz \
        -h ${QC_DIR}/${SAMPLE}_fastp.html \
        -j ${QC_DIR}/${SAMPLE}_fastp.json \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread ${THREADS}

    # STAR比对
    STAR --runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn \
             ${CLEAN_DIR}/${SAMPLE}_R1_clean.fastq.gz \
             ${CLEAN_DIR}/${SAMPLE}_R2_clean.fastq.gz \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${ALIGN_DIR}/${SAMPLE}_ \
         --quantMode GeneCounts \
         --outSAMattributes NH HI AS NM MD

    samtools index ${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam
    echo "========== Done: ${SAMPLE} =========="
done

# --- featureCounts定量 ---
BAM_FILES=$(ls ${ALIGN_DIR}/*_Aligned.sortedByCoord.out.bam | tr '\n' ' ')

featureCounts \
    -a ${GTF} \
    -o ${COUNTS_DIR}/self_PHPT_counts_raw.txt \
    -T ${THREADS} \
    -p --countReadPairs \
    -s 2 \
    -g gene_id \
    -t exon \
    --extraAttributes gene_name \
    ${BAM_FILES}

cut -f1,8- ${COUNTS_DIR}/self_PHPT_counts_raw.txt | tail -n +2 > ${COUNTS_DIR}/self_PHPT_counts.txt

# --- MultiQC ---
multiqc ${QC_DIR} ${ALIGN_DIR} ${COUNTS_DIR} -o ${QC_DIR}/multiqc_self_PHPT

echo "Self PHPT pipeline complete."
echo "Counts: ${COUNTS_DIR}/self_PHPT_counts.txt"

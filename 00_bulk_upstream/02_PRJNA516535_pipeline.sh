#!/bin/bash
# ============================================================
# 02_PRJNA516535_pipeline.sh
# 公开数据集 PRJNA516535（10个PHPT样本）
# SRA下载 → fastq-dump → fastp → STAR → featureCounts
# ============================================================

set -euo pipefail

# --- 路径配置 ---
WORK_DIR="/home/oyh/project/bulk"
SRA_DIR="${WORK_DIR}/raw_data/PRJNA516535/sra"
RAW_DIR="${WORK_DIR}/raw_data/PRJNA516535/fastq"
CLEAN_DIR="${WORK_DIR}/clean_data/PRJNA516535"
ALIGN_DIR="${WORK_DIR}/aligned/PRJNA516535"
COUNTS_DIR="${WORK_DIR}/counts/PRJNA516535"
QC_DIR="${WORK_DIR}/qc/PRJNA516535"
STAR_INDEX="${WORK_DIR}/reference/star_index"
GTF="${WORK_DIR}/reference/gencode.v44.annotation.gtf"
THREADS=16

mkdir -p ${SRA_DIR} ${RAW_DIR} ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR} ${QC_DIR}

# --- SRR列表（PRJNA516535 / GSE125433） ---
# 文献：Comparative Gene Expression Profiles in Parathyroid Adenoma
#        and Normal Parathyroid Tissue (PMID: 30832348)
# 共15样本：10 parathyroid adenoma + 5 normal parathyroid
# 测序平台：Illumina HiSeq2500, TruSeq RNA Access Library Prep
# SRP编号：SRP175300
#
# 【获取SRR编号】在服务器上运行以下任一命令：
#   pysradb srp-to-srr SRP175300 --saveto PRJNA516535_runs.tsv
#   或访问: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA516535
#   下载 Accession List，筛选10个adenoma样本的SRR编号填入下方
#
SRR_LIST=(
    SRR8475556 SRR8475557 SRR8475558 SRR8475559 SRR8475560
    SRR8475561 SRR8475562 SRR8475563 SRR8475564 SRR8475566
)

# --- 1) 批量下载SRA并转fastq ---
for SRR in "${SRR_LIST[@]}"; do
    echo "========== Downloading: ${SRR} =========="

    # prefetch下载.sra文件
    prefetch ${SRR} -O ${SRA_DIR}

    # 转换为fastq（双端）
    fasterq-dump ${SRA_DIR}/${SRR}/${SRR}.sra \
        -O ${RAW_DIR} \
        -e ${THREADS} \
        --split-3

    # 压缩fastq
    pigz -p ${THREADS} ${RAW_DIR}/${SRR}_1.fastq
    pigz -p ${THREADS} ${RAW_DIR}/${SRR}_2.fastq
done

# --- 2) 质控 + 比对 ---
for SRR in "${SRR_LIST[@]}"; do
    echo "========== Processing: ${SRR} =========="
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"

    # fastp质控
    fastp \
        -i ${R1} -I ${R2} \
        -o ${CLEAN_DIR}/${SRR}_R1_clean.fastq.gz \
        -O ${CLEAN_DIR}/${SRR}_R2_clean.fastq.gz \
        -h ${QC_DIR}/${SRR}_fastp.html \
        -j ${QC_DIR}/${SRR}_fastp.json \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread ${THREADS}

    # STAR比对
    STAR --runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn \
             ${CLEAN_DIR}/${SRR}_R1_clean.fastq.gz \
             ${CLEAN_DIR}/${SRR}_R2_clean.fastq.gz \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${ALIGN_DIR}/${SRR}_ \
         --quantMode GeneCounts \
         --outSAMattributes NH HI AS NM MD

    samtools index ${ALIGN_DIR}/${SRR}_Aligned.sortedByCoord.out.bam

    echo "========== Done: ${SRR} =========="
done

# --- 3) featureCounts批量定量 ---
BAM_FILES=$(ls ${ALIGN_DIR}/*_Aligned.sortedByCoord.out.bam | tr '\n' ' ')

featureCounts \
    -a ${GTF} \
    -o ${COUNTS_DIR}/PRJNA516535_counts_raw.txt \
    -T ${THREADS} \
    -p --countReadPairs \
    -s 2 \
    -g gene_id \
    -t exon \
    --extraAttributes gene_name \
    ${BAM_FILES}

cut -f1,8- ${COUNTS_DIR}/PRJNA516535_counts_raw.txt | tail -n +2 > ${COUNTS_DIR}/PRJNA516535_counts.txt

# --- MultiQC汇总 ---
multiqc ${QC_DIR} ${ALIGN_DIR} ${COUNTS_DIR} -o ${QC_DIR}/multiqc_PRJNA516535

echo "PRJNA516535 pipeline complete."
echo "Counts matrix: ${COUNTS_DIR}/PRJNA516535_counts.txt"

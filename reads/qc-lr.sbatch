#!/bin/bash -l
#SBATCH -J qc-lr
#SBATCH --mem 16gb
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%A_%a.out
#SBATCH --error  logs/%x-%A_%a.err
#
# Long-read QC with fastplong: 1 kb minimum length plus its default N and
# quality filters. Every tool reads the same QC'd reads.
#   holdout  one array task per haplotype (10)
#   wol2     100 genomes per task (160 tasks)
# Output: ${READS_DIR}/long-read/<set>/qc/<name>.fastq.gz
#
# Usage (from the repo root):
#   sbatch --array=1-10  reads/qc-lr.sbatch holdout
#   sbatch --array=1-160 reads/qc-lr.sbatch wol2

source ./common.sh
conda_activate hf-qc

SET="${1:?usage: qc-lr.sbatch <holdout|wol2>}"
THREADS=${SLURM_CPUS_PER_TASK:-16}
OUT_DIR="${READS_DIR}/long-read/${SET}/qc"
mkdir -p "${OUT_DIR}"

case "${SET}" in
    holdout)
        base=$(awk -F'\t' -v n="$(( SLURM_ARRAY_TASK_ID + 1 ))" 'NR == n { print $1 }' "${REPO}/references/holdout10.tsv")
        IN="${READS_DIR}/long-read/holdout/raw/${base}.fastq.gz"; OUT="${OUT_DIR}/${base}.fastq.gz"
        require_files "${IN}"
        [[ -s ${OUT} ]] && { logger "${OUT} exists, skipping"; exit 0; }
        fastplong -i "${IN}" -o "${OUT}" --length_required 1000 --thread "${THREADS}" \
            -j "${OUT_DIR}/${base}.json" -h "${OUT_DIR}/${base}.html" ;;
    wol2)
        start=$(( (SLURM_ARRAY_TASK_ID - 1) * 100 + 1 ))
        while read -r acc; do
            IN="${READS_DIR}/long-read/wol2/sim/${acc}.fastq.gz"; OUT="${OUT_DIR}/${acc}.fastq.gz"
            [[ -s ${OUT} ]] && continue
            [[ -e ${IN} ]] || { logger "no simulated reads for ${acc}"; continue; }
            fastplong -i "${IN}" -o "${OUT}" --length_required 1000 --thread "${THREADS}" \
                --json /dev/null --html /dev/null \
                || { logger "fastplong failed: ${acc}"; rm -f "${OUT}"; }
        done < <(sed -n "${start},$(( start + 99 ))p" "${REPO}/reads/wol2-genomes.list") ;;
    *)  logger "unknown set: ${SET}"; exit 1 ;;
esac
logger "done"

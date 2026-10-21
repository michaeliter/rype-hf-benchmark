#!/bin/bash -l
#SBATCH -J qc-sr
#SBATCH --mem 16gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%A_%a.out
#SBATCH --error  logs/%x-%A_%a.err
#
# Short-read QC with fastp: adapter trimming against known_adapters.fna and a
# 100 bp minimum length, mates kept or dropped together. Every tool reads the
# same QC'd pairs.
#   holdout  one array task per haplotype (10)
#   wol2     WoL2 per-genome pairs from WOL2_SR_READS_DIR, 100 genomes per
#            task (reads/wol2-genomes.list: 15,953 genomes, 160 tasks)
# Output: ${READS_DIR}/short-read/<set>/qc/<name>_R{1,2}.fastq.gz
#
# Usage (from the repo root):
#   sbatch --array=1-10  reads/qc-sr.sbatch holdout
#   sbatch --array=1-160 reads/qc-sr.sbatch wol2

source ./common.sh
conda_activate human-filtration

SET="${1:?usage: qc-sr.sbatch <holdout|wol2>}"
THREADS=${SLURM_CPUS_PER_TASK:-8}
OUT_DIR="${READS_DIR}/short-read/${SET}/qc"
mkdir -p "${OUT_DIR}"
require_files "${ADAPTERS}"

qc_pair () {  # <name> <r1> <r2> <json>
    local o1="${OUT_DIR}/$1_R1.fastq.gz" o2="${OUT_DIR}/$1_R2.fastq.gz"
    [[ -s ${o1} && -s ${o2} ]] && { logger "$1 exists, skipping"; return 0; }
    [[ -e $2 && -e $3 ]] || { logger "missing input pair for $1"; return 1; }
    fastp -i "$2" -I "$3" -o "${o1}" -O "${o2}" -l 100 -w "${THREADS}" \
        --adapter_fasta "${ADAPTERS}" --json "$4" --html /dev/null \
        || { rm -f "${o1}" "${o2}"; return 1; }
}

case "${SET}" in
    holdout)
        base=$(awk -F'\t' -v n="$(( SLURM_ARRAY_TASK_ID + 1 ))" 'NR == n { print $1 }' "${REPO}/references/holdout10.tsv")
        RAW="${READS_DIR}/short-read/holdout/raw"
        qc_pair "${base}" "${RAW}/${base}_R1.fastq.gz" "${RAW}/${base}_R2.fastq.gz" "${OUT_DIR}/${base}.fastp.json" ;;
    wol2)
        # A genome that fails fails the task, so a resubmit refills the gap.
        start=$(( (SLURM_ARRAY_TASK_ID - 1) * 100 + 1 ))
        fail=0
        while read -r acc; do
            qc_pair "${acc}" "${WOL2_SR_READS_DIR}/${acc}_R1.fastq.gz" "${WOL2_SR_READS_DIR}/${acc}_R2.fastq.gz" /dev/null \
                || { logger "fastp failed: ${acc}"; fail=$((fail + 1)); }
        done < <(sed -n "${start},$(( start + 99 ))p" "${REPO}/reads/wol2-genomes.list")
        (( fail == 0 )) || { logger "${fail} genome(s) failed; resubmit this task"; exit 1; } ;;
    *)  logger "unknown set: ${SET}"; exit 1 ;;
esac
logger "done"

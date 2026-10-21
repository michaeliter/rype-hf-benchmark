#!/bin/bash -l
#SBATCH -J build-concat-inputs
#SBATCH --mem 32gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 16:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# Builds the benchmark inputs: one concatenated FASTQ (pair) per modality and
# set from the QC'd reads, every header rewritten from "@<id> ..." to
# "@<tag>::<id> ...":
#   holdout  tag = haplotype (references/holdout10.tsv)
#   wol2     tag = chunk_NN, from reads/wol2-chunks/<modality>/chunk_NN.list
# Every tool then filters each whole set in one call (the index is loaded
# once), and FN/FP are counted per tag from the survivors. "::" contains no
# whitespace, so the tag stays part of the read ID through every tool.
#
#   ${CONCAT_DIR}/short-read/{holdout,wol2}_R{1,2}.fastq.gz
#   ${CONCAT_DIR}/long-read/{holdout,wol2}.fastq.gz
#
# Put CONCAT_DIR on fast parallel storage: every benchmark arm streams from it.
#
# Usage (from the repo root): sbatch reads/build-concat-inputs.sbatch

source ./common.sh
mkdir -p "${CONCAT_SR}" "${CONCAT_LR}"

relabel () {  # <tag> <src.fastq.gz>  -> stdout
    zcat -- "$2" | awk -v tag="$1" 'NR % 4 == 1 { sub(/^@/, ""); print "@" tag "::" $0; next } { print }'
}
holdout_names () { awk -F'\t' 'NR > 1 { print $1 }' "${REPO}/references/holdout10.tsv"; }

sr_holdout () {
    local out1="${CONCAT_SR}/holdout_R1.fastq.gz" out2="${CONCAT_SR}/holdout_R2.fastq.gz" qc="${READS_DIR}/short-read/holdout/qc"
    [[ -s ${out1} && -s ${out2} ]] && return 0
    : > "${out1}.tmp"; : > "${out2}.tmp"
    for base in $(holdout_names); do
        relabel "${base}" "${qc}/${base}_R1.fastq.gz" | gzip -c >> "${out1}.tmp"
        relabel "${base}" "${qc}/${base}_R2.fastq.gz" | gzip -c >> "${out2}.tmp"
    done
    mv "${out1}.tmp" "${out1}"; mv "${out2}.tmp" "${out2}"
}
sr_wol2 () {
    local out1="${CONCAT_SR}/wol2_R1.fastq.gz" out2="${CONCAT_SR}/wol2_R2.fastq.gz" qc="${READS_DIR}/short-read/wol2/qc" list acc
    [[ -s ${out1} && -s ${out2} ]] && return 0
    : > "${out1}.tmp"; : > "${out2}.tmp"
    for list in "${REPO}"/reads/wol2-chunks/short-read/chunk_*.list; do
        while read -r acc; do
            [[ -s ${qc}/${acc}_R1.fastq.gz && -s ${qc}/${acc}_R2.fastq.gz ]] || { logger "missing QC'd reads for ${acc}"; return 1; }
            relabel "$(basename "${list}" .list)" "${qc}/${acc}_R1.fastq.gz" | gzip -c >> "${out1}.tmp"
            relabel "$(basename "${list}" .list)" "${qc}/${acc}_R2.fastq.gz" | gzip -c >> "${out2}.tmp"
        done < "${list}"
    done
    mv "${out1}.tmp" "${out1}"; mv "${out2}.tmp" "${out2}"
}
lr_holdout () {
    local out="${CONCAT_LR}/holdout.fastq.gz" qc="${READS_DIR}/long-read/holdout/qc"
    [[ -s ${out} ]] && return 0
    : > "${out}.tmp"
    for base in $(holdout_names); do
        relabel "${base}" "${qc}/${base}.fastq.gz" | gzip -c >> "${out}.tmp"
    done
    mv "${out}.tmp" "${out}"
}
lr_wol2 () {
    local out="${CONCAT_LR}/wol2.fastq.gz" qc="${READS_DIR}/long-read/wol2/qc" list acc
    [[ -s ${out} ]] && return 0
    : > "${out}.tmp"
    for list in "${REPO}"/reads/wol2-chunks/long-read/chunk_*.list; do
        while read -r acc; do
            [[ -s ${qc}/${acc}.fastq.gz ]] || { logger "missing QC'd reads for ${acc}"; return 1; }
            relabel "$(basename "${list}" .list)" "${qc}/${acc}.fastq.gz" | gzip -c >> "${out}.tmp"
        done < "${list}"
    done
    mv "${out}.tmp" "${out}"
}

# The four targets share no files: build them in parallel.
declare -A PID=()
for step in sr_holdout sr_wol2 lr_holdout lr_wol2; do
    ( set -euo pipefail; "${step}" ) & PID[${step}]=$!
done
failed=0
for step in sr_holdout sr_wol2 lr_holdout lr_wol2; do
    if wait "${PID[${step}]}"; then logger "${step}: OK"; else logger "${step}: FAILED"; failed=1; fi
done
(( failed == 0 )) || exit 1
logger "inputs built under ${CONCAT_DIR}"

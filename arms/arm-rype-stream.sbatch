#!/bin/bash -l
#SBATCH -J arm-rype-stream
#SBATCH --mem 64gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# rype and rype -> minimap2 through rype-hf, the streaming pipeline
# (rype classify filter | minimap2 | samtools fastq | bgzip -l1), 8 CPU / 64 GB.
#
# Timed region = the whole pipeline, gzipped FASTQ in -> BGZF survivor FASTQ
# out. Short-read survivors are written as separate R1/R2 files (rype-hf -O),
# long-read survivors as one file. Recorded per run:
#   tree_rss.tsv   wall seconds + peak summed RSS of the whole process tree
#   time-v.stderr  GNU time -v (max RSS of the single largest process)
#   stats.tsv      read counts per stage, from rype-hf
#   host.txt       node + rype commit
# Only replicate 1 keeps its survivors (for FN/FP).
#
# Window sweep: an optional 5th argument runs another window's index instead
# of w=20, into <set>/<mode>/w<W>/. Those runs count survivors per source tag
# (tag_counts.tsv, from R1 for short reads) outside the timed region and then
# delete the survivors. w=20 is the operating-point rep1 run.
#
# Usage (from the repo root):
#   sbatch arms/arm-rype-stream.sbatch <sr|lr> <holdout|wol2> <rype|rype-mm2> <rep> [window]

source ./common.sh

MOD="${1:?usage: arm-rype-stream.sbatch <sr|lr> <holdout|wol2> <rype|rype-mm2> <rep> [window]}"
SET="${2:?missing set}"
MODE="${3:?missing mode}"
REP="${4:-1}"
W="${5:-${OP_WINDOW}}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
INDEX="$(rype_index "${W}")"

case "${MOD}" in
    sr) INPUT_ARGS=(-1 "${CONCAT_SR}/${SET}_R1.fastq.gz" -2 "${CONCAT_SR}/${SET}_R2.fastq.gz")
        SURV=(survivors_R1.fastq.gz survivors_R2.fastq.gz)
        THRESHOLD=${SR_THRESHOLD}; MMI=${SR_RESCUE_MMI}; PRESET=sr; MOD_DIR=short-read ;;
    lr) INPUT_ARGS=(-1 "${CONCAT_LR}/${SET}.fastq.gz")
        SURV=(survivors.fastq.gz)
        THRESHOLD=${LR_THRESHOLD}; MMI=${LR_RESCUE_MMI}; PRESET=map-hifi; MOD_DIR=long-read ;;
    *)  logger "unknown modality: ${MOD}"; exit 1 ;;
esac
case "${SET}" in holdout|wol2) ;; *) logger "unknown set: ${SET}"; exit 1 ;; esac
RESCUE_ARGS=()
case "${MODE}" in
    rype)     ;;
    rype-mm2) RESCUE_ARGS=(-m "${MMI}" -x "${PRESET}") ;;
    *)        logger "unknown mode: ${MODE}"; exit 1 ;;
esac
require_files "${RYPE_BIN}" "${RYPE_HF}" "${INDEX}/manifest.toml" "${INPUT_ARGS[1]}" "${INPUT_ARGS[@]: -1}"
[[ ${MODE} == rype ]] || require_files "${MMI}"

if [[ ${W} == "${OP_WINDOW}" ]]; then
    OUT_DIR="${RESULTS_DIR}/${MOD_DIR}/${SET}/${MODE}/rep${REP}"
else
    OUT_DIR="${RESULTS_DIR}/${MOD_DIR}/${SET}/${MODE}/w${W}"
fi
mkdir -p "${OUT_DIR}"
if [[ -s ${OUT_DIR}/tree_rss.tsv ]] && awk 'NR==2 && $3==0 {ok=1} END {exit !ok}' "${OUT_DIR}/tree_rss.tsv" \
        && [[ ${W} == "${OP_WINDOW}" || -s ${OUT_DIR}/tag_counts.tsv ]]; then
    logger "${MOD} ${SET} ${MODE} rep${REP} w${W}: already complete, skipping"
    exit 0
fi

logger "${MOD} ${SET} ${MODE} rep${REP} on $(hostname): w=${W} t=${THRESHOLD} threads=${THREADS}"
RYPE_COMMIT=$(git -C "$(dirname "${RYPE_BIN}")/../.." rev-parse --short HEAD 2>/dev/null || echo '?')
printf "host\t%s\nslurm_job\t%s\nrype_commit\t%s\nindex\t%s\n" "$(hostname)" "${SLURM_JOB_ID:-NA}" "${RYPE_COMMIT}" "${INDEX}" > "${OUT_DIR}/host.txt"

# Per-run scratch for rype-hf's small stage logs (no read data goes here).
job_scratch arm-rype-stream
OUTPUT_ARGS=(-o "${OUT_DIR}/${SURV[0]}")
[[ ${MOD} == sr ]] && OUTPUT_ARGS+=(-O "${OUT_DIR}/${SURV[1]}")

"${PY}" "${TREE_RSS}" "${OUT_DIR}/tree_rss.tsv" -- \
    ${BINTIME} -o "${OUT_DIR}/time-v.stderr" \
    "${RYPE_HF}" -i "${INDEX}" "${INPUT_ARGS[@]}" -t "${THRESHOLD}" \
        "${RESCUE_ARGS[@]}" -p "${THREADS}" -M "${RYPE_MAX_MEMORY}" \
        -s "${OUT_DIR}/stats.tsv" "${OUTPUT_ARGS[@]}"

cat "${OUT_DIR}/tree_rss.tsv" >&2
cat "${OUT_DIR}/stats.tsv" >&2
if [[ ${W} != "${OP_WINDOW}" ]]; then
    # untimed: per-tag survivor counts (R1 for short reads: one record per
    # pair), then drop the survivors
    if [[ -s ${OUT_DIR}/${SURV[0]} ]]; then
        "${BGZIP_BIN}" -dc -@ "${THREADS}" "${OUT_DIR}/${SURV[0]}"
    fi \
        | awk 'NR % 4 == 1 { split($0, a, "::"); c[substr(a[1], 2)]++ }
               END { print "tag\treads"; for (k in c) print k "\t" c[k] }' \
        > "${OUT_DIR}/tag_counts.tsv.tmp"
    mv "${OUT_DIR}/tag_counts.tsv.tmp" "${OUT_DIR}/tag_counts.tsv"
    rm -f "${OUT_DIR}"/survivors*.fastq.gz
elif [[ ${REP} -ne 1 ]]; then
    rm -f "${OUT_DIR}"/survivors*.fastq.gz
fi
logger "${MOD} ${SET} ${MODE} rep${REP}: done"

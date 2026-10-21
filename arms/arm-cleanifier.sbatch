#!/bin/bash -l
#SBATCH -J arm-cleanifier
#SBATCH --mem 64gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# cleanifier, run the way its paper and docs recommend, with either index
# (references/build-cleanifier-index.sbatch):
#   exact    (3,4) bucketed Cuckoo hash table      -> <set>/cleanifier/
#   filter   windowed Cuckoo filter, FPR 2^-14     -> <set>/cleanifier-filter/
# The filter call is common.sh's cleanifier_args: threshold 0.5, sampling mode
# (both defaults), 8 threads, gz output at cleanifier's defaults, --buffersize
# 22 for long reads only.
#
# Measured like the rype arms: tree_rss.tsv (wall + peak summed RSS of the
# process tree, which includes the bgzip processes cleanifier spawns) and
# time-v.stderr. Kept reads stay in the result directory for the FN/FP counts.
#
# Usage (from the repo root): sbatch arms/arm-cleanifier.sbatch <exact|filter> <sr|lr> <holdout|wol2>

source ./common.sh

KIND="${1:?usage: arm-cleanifier.sbatch <exact|filter> <sr|lr> <holdout|wol2>}"
MOD="${2:?missing modality}"
SET="${3:?missing set}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

case "${KIND}" in
    exact)  INDEX=${CLEANIFIER_EXACT_IDX};  APPROACH=cleanifier;        require_files "${INDEX}.hash" "${INDEX}.info" ;;
    filter) INDEX=${CLEANIFIER_FILTER_IDX}; APPROACH=cleanifier-filter; require_files "${INDEX}.filter" "${INDEX}.info" ;;
    *)      logger "unknown index kind: ${KIND}"; exit 1 ;;
esac
case "${MOD}" in
    sr) R1="${CONCAT_SR}/${SET}_R1.fastq.gz"; R2="${CONCAT_SR}/${SET}_R2.fastq.gz"; MOD_DIR=short-read ;;
    lr) R1="${CONCAT_LR}/${SET}.fastq.gz";    R2=;                                 MOD_DIR=long-read ;;
    *)  logger "unknown modality: ${MOD}"; exit 1 ;;
esac
case "${SET}" in holdout|wol2) ;; *) logger "unknown set: ${SET}"; exit 1 ;; esac
require_files "${R1}" ${R2:+"${R2}"}

OUT_DIR="${RESULTS_DIR}/${MOD_DIR}/${SET}/${APPROACH}"
mkdir -p "${OUT_DIR}"
if [[ -s ${OUT_DIR}/tree_rss.tsv ]] && awk 'NR==2 && $3==0 {ok=1} END {exit !ok}' "${OUT_DIR}/tree_rss.tsv"; then
    logger "${APPROACH} ${MOD} ${SET}: already complete, skipping"
    exit 0
fi
rm -f "${OUT_DIR}"/clean_*

printf "host\t%s\nslurm_job\t%s\nindex\t%s\n" "$(hostname)" "${SLURM_JOB_ID:-NA}" "${INDEX}" > "${OUT_DIR}/host.txt"
cleanifier_args "${MOD}" "${INDEX}" "${R1}" "${R2}" "${OUT_DIR}/clean" "${THREADS}"
logger "${APPROACH} ${MOD} ${SET} on $(hostname): cleanifier ${CLEANIFIER_ARGS[*]}"

conda_activate cleanifier
"${PY}" "${TREE_RSS}" "${OUT_DIR}/tree_rss.tsv" -- \
    ${BINTIME} -o "${OUT_DIR}/time-v.stderr" cleanifier "${CLEANIFIER_ARGS[@]}" > /dev/null

ls "${OUT_DIR}"/clean_*keep* > /dev/null || { logger "cleanifier wrote no keep file"; exit 1; }
cat "${OUT_DIR}/tree_rss.tsv" >&2
logger "${APPROACH} ${MOD} ${SET}: done"

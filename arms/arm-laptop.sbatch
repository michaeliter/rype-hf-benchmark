#!/bin/bash -l
#SBATCH -J arm-laptop
#SBATCH --mem 16gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# Laptop-class arm: 8 CPU and a small memory allocation, same concatenated
# inputs and indices as the 64 GB runs. The tier is the job's --mem, set by
# run-laptop.sh (16 GB by default here; 8 GB for the long-read-only tier), and
# results go to laptop-<N>gb/. Only tools whose footprint fits are run; Movi2
# and mm2->mm2->Movi (~37 GB) do not fit and are reported as such, and neither
# does cleanifier exact at 8 GB (its table is ~15 GB).
#
#   rype       rype-hf, rype only
#   rype-mm2   rype-hf, rype -> minimap2, streaming as in the 64 GB arm. The
#              full minimap2 index is resident for the whole run, so rype-hf's
#              auto budget gives rype 90% of what is left minus 512 MiB
#              (~3.6 GiB for SR).
#              Short read at 16 GB only: minimap2's long-read (T2T map-hifi)
#              index alone peaks at 13.8 GB, the SR index at 11.4 GB.
#   cleanifier          exact index      } same call as arm-cleanifier.sbatch
#   cleanifier-filter   FPR 2^-14 filter } (common.sh cleanifier_args)
#
# Timed region and metrics match arm-rype-stream: tree_rss.tsv (wall + peak
# summed RSS of the process tree) and time-v.stderr. Accuracy comes from the
# 64 GB runs (the outputs do not depend on the memory budget), so survivors
# go to per-job scratch and are deleted; rype's stats.tsv is kept so its
# counts can be checked against the 64 GB run.
#
# Usage (from the repo root):
#   sbatch [--mem 8gb] arms/arm-laptop.sbatch <sr|lr> <holdout|wol2> <rype|rype-mm2|cleanifier|cleanifier-filter>

source ./common.sh

MOD="${1:?usage: arm-laptop.sbatch <sr|lr> <holdout|wol2> <tool>}"
SET="${2:?missing set}"
TOOL="${3:?missing tool}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
RYPE_INDEX_W20="$(rype_index "${OP_WINDOW}")"

case "${MOD}" in
    sr) R1="${CONCAT_SR}/${SET}_R1.fastq.gz"; R2="${CONCAT_SR}/${SET}_R2.fastq.gz"
        THRESHOLD=${SR_THRESHOLD}; MMI=${SR_RESCUE_MMI}; PRESET=sr; MOD_DIR=short-read ;;
    lr) R1="${CONCAT_LR}/${SET}.fastq.gz"; R2=
        THRESHOLD=${LR_THRESHOLD}; MMI=${LR_RESCUE_MMI}; PRESET=map-hifi; MOD_DIR=long-read ;;
    *)  logger "unknown modality: ${MOD}"; exit 1 ;;
esac
case "${SET}" in holdout|wol2) ;; *) logger "unknown set: ${SET}"; exit 1 ;; esac
require_files "${R1}" ${R2:+"${R2}"}
MEM_GB=$(( ${SLURM_MEM_PER_NODE:?run under sbatch, which sets the memory tier} / 1024 ))
if [[ ${TOOL} == rype-mm2 ]] && [[ ${MOD} == lr || ${MEM_GB} -lt 16 ]]; then
    logger "${MOD} rype-mm2 does not fit in ${MEM_GB} GB (minimap2 alone peaks at 11.4 GB for SR, 13.8 GB for LR)"
    exit 1
fi

OUT_DIR="${RESULTS_DIR}/laptop-${MEM_GB}gb/${MOD_DIR}/${SET}/${TOOL}"
mkdir -p "${OUT_DIR}"
if [[ -s ${OUT_DIR}/tree_rss.tsv ]] && awk 'NR==2 && $3==0 {ok=1} END {exit !ok}' "${OUT_DIR}/tree_rss.tsv"; then
    logger "laptop-${MEM_GB}gb ${MOD} ${SET} ${TOOL}: already complete, skipping"
    exit 0
fi

job_scratch arm-laptop

RYPE_COMMIT=$(git -C "$(dirname "${RYPE_BIN}")/../.." rev-parse --short HEAD 2>/dev/null || echo '?')
printf "host\t%s\nslurm_job\t%s\nmem\t%s\nrype_commit\t%s\n" "$(hostname)" "${SLURM_JOB_ID:-NA}" \
    "${SLURM_MEM_PER_NODE:-NA}" "${RYPE_COMMIT}" > "${OUT_DIR}/host.txt"
logger "laptop-${MEM_GB}gb ${MOD} ${SET} ${TOOL} on $(hostname), threads=${THREADS}, mem=${SLURM_MEM_PER_NODE}MB"

timed () { "${PY}" "${TREE_RSS}" "${OUT_DIR}/tree_rss.tsv" -- ${BINTIME} -o "${OUT_DIR}/time-v.stderr" "$@"; }

case "${TOOL}" in
    rype|rype-mm2)
        require_files "${RYPE_BIN}" "${RYPE_HF}" "${RYPE_INDEX_W20}"
        if [[ ${TOOL} == rype ]]; then
            EXTRA=(-M "${RYPE_MAX_MEMORY}")
        else
            require_files "${MMI}"
            EXTRA=(-M "${RYPE_MAX_MEMORY}" -m "${MMI}" -x "${PRESET}")
        fi
        # short reads: separate R1/R2 survivors, as in the 64 GB arm
        if [[ -n ${R2} ]]; then
            OUTS=(-o "${TMPDIR}/survivors_R1.fastq.gz" -O "${TMPDIR}/survivors_R2.fastq.gz")
        else
            OUTS=(-o "${TMPDIR}/survivors.fastq.gz")
        fi
        timed "${RYPE_HF}" -i "${RYPE_INDEX_W20}" -1 "${R1}" ${R2:+-2 "${R2}"} -t "${THRESHOLD}" \
            "${EXTRA[@]}" -p "${THREADS}" -s "${OUT_DIR}/stats.tsv" "${OUTS[@]}"
        cat "${OUT_DIR}/stats.tsv" >&2 ;;
    cleanifier|cleanifier-filter)
        if [[ ${TOOL} == cleanifier ]]; then
            INDEX=${CLEANIFIER_EXACT_IDX}; require_files "${INDEX}.hash" "${INDEX}.info"
        else
            INDEX=${CLEANIFIER_FILTER_IDX}; require_files "${INDEX}.filter" "${INDEX}.info"
        fi
        cleanifier_args "${MOD}" "${INDEX}" "${R1}" "${R2}" "${TMPDIR}/clean" "${THREADS}"
        conda_activate cleanifier
        timed cleanifier "${CLEANIFIER_ARGS[@]}" > /dev/null ;;
    *)  logger "unknown tool: ${TOOL}"; exit 1 ;;
esac

cat "${OUT_DIR}/tree_rss.tsv" >&2
logger "laptop-${MEM_GB}gb ${MOD} ${SET} ${TOOL}: done"

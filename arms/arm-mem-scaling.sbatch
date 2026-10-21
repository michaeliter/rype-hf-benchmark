#!/bin/bash -l
#SBATCH -J arm-mem-scaling
#SBATCH --mem 64gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# Memory scaling of long-read rype (the long-read operating point: rype alone,
# w=20, t=0.05) at 8 CPU and a range of job memory limits, for the
# supplementary memory-scaling figure. The tier is the job's --mem, set by
# run-mem-scaling.sh (64 down to 1 GB). rype-hf sizes rype's --max-memory from
# the job limit with the same rule as every other run: 90% minus 512 MiB.
# Long read only: the short-read operating point streams through minimap2,
# whose index alone needs ~11 GB.
#
# Recorded per run, in mem-scaling/long-read/<set>/<N>gb/rep<R>/:
#   tree_rss.tsv     wall seconds, peak summed RSS of the process tree, exit code
#   time-v.stderr    /usr/bin/time -v
#   stats.tsv        read counts (checked against the 64 GB operating-point run)
#   rype-hf.stderr   rype-hf + rype stderr: --max-memory, pass budget, passes
#   host.txt         node, job, memory limit, rype commit
# Survivors go to per-job scratch and are deleted. A run that fails or is
# OOM-killed keeps its tree_rss.tsv (non-zero exit) and is reported as such.
#
# Usage (from the repo root): sbatch --mem <N>gb arms/arm-mem-scaling.sbatch <holdout|wol2> <rep>

source ./common.sh

SET="${1:?usage: arm-mem-scaling.sbatch <holdout|wol2> <rep>}"
REP="${2:?missing rep}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
RYPE_INDEX_W20="$(rype_index "${OP_WINDOW}")"
case "${SET}" in holdout|wol2) ;; *) logger "unknown set: ${SET}"; exit 1 ;; esac
R1="${CONCAT_LR}/${SET}.fastq.gz"
require_files "${R1}" "${RYPE_BIN}" "${RYPE_HF}" "${RYPE_INDEX_W20}/manifest.toml"
MEM_GB=$(( ${SLURM_MEM_PER_NODE:?run under sbatch, which sets the memory tier} / 1024 ))

OUT_DIR="${RESULTS_DIR}/mem-scaling/long-read/${SET}/${MEM_GB}gb/rep${REP}"
mkdir -p "${OUT_DIR}"
if [[ -s ${OUT_DIR}/tree_rss.tsv ]] && awk 'NR==2 && $3==0 {ok=1} END {exit !ok}' "${OUT_DIR}/tree_rss.tsv"; then
    logger "mem-scaling lr ${SET} ${MEM_GB}gb rep${REP}: already complete, skipping"
    exit 0
fi
rm -f "${OUT_DIR}"/{tree_rss.tsv,time-v.stderr,stats.tsv,rype-hf.stderr}

job_scratch arm-mem-scaling

RYPE_COMMIT=$(git -C "$(dirname "${RYPE_BIN}")/../.." rev-parse --short HEAD 2>/dev/null || echo '?')
printf "host\t%s\nslurm_job\t%s\nmem_mb\t%s\nrype_commit\t%s\n" "$(hostname)" "${SLURM_JOB_ID:-NA}" \
    "${SLURM_MEM_PER_NODE}" "${RYPE_COMMIT}" > "${OUT_DIR}/host.txt"
logger "mem-scaling lr ${SET} ${MEM_GB}gb rep${REP} on $(hostname), threads=${THREADS}"

status=0
"${PY}" "${TREE_RSS}" "${OUT_DIR}/tree_rss.tsv" -- \
    ${BINTIME} -o "${OUT_DIR}/time-v.stderr" \
    "${RYPE_HF}" -i "${RYPE_INDEX_W20}" -1 "${R1}" -t "${LR_THRESHOLD}" -M "${RYPE_MAX_MEMORY}" \
        -p "${THREADS}" -s "${OUT_DIR}/stats.tsv" -o "${TMPDIR}/survivors.fastq.gz" \
    2> "${OUT_DIR}/rype-hf.stderr" || status=$?

cat "${OUT_DIR}/rype-hf.stderr" >&2
[[ -s ${OUT_DIR}/tree_rss.tsv ]] && cat "${OUT_DIR}/tree_rss.tsv" >&2
logger "mem-scaling lr ${SET} ${MEM_GB}gb rep${REP}: exit ${status}"
exit "${status}"

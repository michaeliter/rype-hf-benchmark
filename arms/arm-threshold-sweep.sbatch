#!/bin/bash -l
#SBATCH -J arm-threshold-sweep
#SBATCH --mem 64gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# rype-alone threshold sweep at w=20 for the supplementary threshold figure.
# Scores once, rescores per threshold: `rype classify run` at the lowest grid
# threshold (0.0001) reports the score of every read or pair that is host at
# any threshold of the grid (classify run and classify filter share rype's
# classification core; pairs are scored as one query in both). DuckDB then
# counts, per threshold and source tag, the reads scoring >= threshold. Scores
# are written as Parquet (full f64): the TSV output rounds to 4 decimals,
# which would move reads across thresholds such as 0.05.
#
# FN/FP only: nothing is timed. eval/eval_op.py checks that the 0.0001 (short read)
# and 0.05 (long read) rows match the operating-point rype-hf runs.
#
# Output: <mod>/<set>/threshold-sweep/
#   scores.parquet   read_id, bucket_name, score (reads scoring >= 0.0001)
#   host_counts.tsv  threshold \t tag \t host (reads or pairs scoring >= threshold)
#   host.txt         node + rype commit
#   time-v.stderr    GNU time -v of classify run (peak memory)
#
# Usage (from the repo root): sbatch arms/arm-threshold-sweep.sbatch <sr|lr> <holdout|wol2>

source ./common.sh

THR_GRID=(0.0001 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5 0.7)
MOD="${1:?usage: arm-threshold-sweep.sbatch <sr|lr> <holdout|wol2>}"
SET="${2:?missing set}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
INDEX="$(rype_index "${OP_WINDOW}")"

case "${MOD}" in
    sr) INPUT_ARGS=(-1 "${CONCAT_SR}/${SET}_R1.fastq.gz" -2 "${CONCAT_SR}/${SET}_R2.fastq.gz")
        MOD_DIR=short-read ;;
    lr) INPUT_ARGS=(-1 "${CONCAT_LR}/${SET}.fastq.gz")
        MOD_DIR=long-read ;;
    *)  logger "unknown modality: ${MOD}"; exit 1 ;;
esac
case "${SET}" in holdout|wol2) ;; *) logger "unknown set: ${SET}"; exit 1 ;; esac
require_files "${RYPE_BIN}" "${DUCKDB_BIN}" "${INDEX}/manifest.toml" "${INPUT_ARGS[1]}" "${INPUT_ARGS[@]: -1}"

OUT_DIR="${RESULTS_DIR}/${MOD_DIR}/${SET}/threshold-sweep"
mkdir -p "${OUT_DIR}"
if [[ -s ${OUT_DIR}/host_counts.tsv ]]; then
    logger "${MOD} ${SET}: already complete, skipping"
    exit 0
fi

# rype's own --max-memory auto misses the Slurm limit (it sits on the job
# cgroup, not the task's), so the budget is explicit: half the job limit.
# classify run sizes its input batches from the budget without the cap that
# classify filter has, and runs past it: at 90% of 64 GB all four sets were
# OOM-killed. The budget only sets the number of passes, not which reads
# score as host.
BUDGET_MIB=$(( ${SLURM_MEM_PER_NODE:-65536} / 2 ))
logger "${MOD} ${SET} on $(hostname): w=${OP_WINDOW} t=${THR_GRID[0]} threads=${THREADS} budget=${BUDGET_MIB} MiB"
RYPE_COMMIT=$(git -C "$(dirname "${RYPE_BIN}")/../.." rev-parse --short HEAD 2>/dev/null || echo '?')
printf "host\t%s\nslurm_job\t%s\nrype_commit\t%s\nindex\t%s\n" "$(hostname)" "${SLURM_JOB_ID:-NA}" \
    "${RYPE_COMMIT}" "${INDEX}" > "${OUT_DIR}/host.txt"

RAYON_NUM_THREADS=${THREADS} ${BINTIME} -o "${OUT_DIR}/time-v.stderr" "${RYPE_BIN}" classify run -i "${INDEX}" "${INPUT_ARGS[@]}" \
    -t "${THR_GRID[0]}" --max-memory "${BUDGET_MIB}M" -o "${OUT_DIR}/scores.tmp.parquet"
mv "${OUT_DIR}/scores.tmp.parquet" "${OUT_DIR}/scores.parquet"

# DuckDB sizes itself from node RAM, not the job limit: cap it and let it spill
# (short-read holdout has a score row for nearly all of its 100M pairs).
job_scratch arm-threshold-sweep
THR_VALUES=$(printf "('%s')," "${THR_GRID[@]}")
"${DUCKDB_BIN}" -c "
SET threads = ${THREADS};
SET memory_limit = '40GB';
SET temp_directory = '${TMPDIR}';
COPY (
    WITH best AS (SELECT read_id, max(score) AS score
                  FROM read_parquet('${OUT_DIR}/scores.parquet') GROUP BY read_id),
         grid(threshold) AS (VALUES ${THR_VALUES%,})
    SELECT g.threshold, split_part(b.read_id, '::', 1) AS tag, count(*) AS host
    FROM best b JOIN grid g ON b.score >= CAST(g.threshold AS DOUBLE)
    GROUP BY g.threshold, tag ORDER BY CAST(g.threshold AS DOUBLE), tag
) TO '${OUT_DIR}/host_counts.tsv.tmp' (FORMAT CSV, DELIMITER '\t', HEADER);"
mv "${OUT_DIR}/host_counts.tsv.tmp" "${OUT_DIR}/host_counts.tsv"
awk -F'\t' 'NR > 1 { n[$1] += $3 } END { for (t in n) print "host at t=" t ": " n[t] }' \
    "${OUT_DIR}/host_counts.tsv" | sort -t= -k2 -g >&2
logger "${MOD} ${SET}: done"

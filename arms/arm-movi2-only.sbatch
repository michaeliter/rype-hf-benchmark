#!/bin/bash -l
#SBATCH -J arm-movi2-only
#SBATCH --mem 64gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 48:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# Movi 2 alone: `movi query --filter -v` against the hprc-yr1 index, no
# minimap2, one call over the whole concatenated set. Movi reads plain
# FASTQ, so the input is first written out uncompressed (short read:
# interleaved with seqtk mergepe; long read: zcat), as its own timed stage.
# gzip of Movi's output is inside the timed call, since every other tool
# compresses its own output. Movi --filter writes FASTA, despite the file
# name; split_survivors.py and eval_common.py detect it.
#
# Output: <mod>/<set>/movi2-only/
#   resource/ALL.tsv   one row per timed stage (wall time, peak RSS)
#   resource/*.time-v.stderr
#   survivors/<tag>.survivors.fastq.gz   split back per source tag
#
# Usage (from the repo root): sbatch arms/arm-movi2-only.sbatch <sr|lr> <holdout|wol2>

source ./common.sh

MOD="${1:?usage: arm-movi2-only.sbatch <sr|lr> <holdout|wol2>}"
SET="${2:?missing set}"
THREADS=${SLURM_CPUS_PER_TASK:-8}
case "${MOD}" in
    sr) R1="${CONCAT_SR}/${SET}_R1.fastq.gz"; R2="${CONCAT_SR}/${SET}_R2.fastq.gz"; MOD_DIR=short-read ;;
    lr) R1="${CONCAT_LR}/${SET}.fastq.gz";    R2=;                                 MOD_DIR=long-read ;;
    *)  logger "unknown modality: ${MOD}"; exit 1 ;;
esac
case "${SET}" in holdout|wol2) ;; *) logger "unknown set: ${SET}"; exit 1 ;; esac
require_files "${R1}" ${R2:+"${R2}"} "${MOVI_BIN}" "${MOVI_INDEX}"

OUT_DIR="${RESULTS_DIR}/${MOD_DIR}/${SET}/movi2-only"
RES_DIR="${OUT_DIR}/resource"; SURV_DIR="${OUT_DIR}/survivors"
mkdir -p "${RES_DIR}" "${SURV_DIR}"
job_scratch "arm-movi2-only-${MOD}-${SET}"
PLAIN="${TMPDIR}/${SET}.fastq"
CONCAT_SURV="${TMPDIR}/${SET}.survivors.fastq.gz"

log_stage () { printf "%s\t%s\t%s\t%s\n" ALL "$1" \
    "$(grep -m1 'Elapsed (wall clock)' "$2" | sed 's/.*): //')" \
    "$(grep -m1 'Maximum resident set size' "$2" | sed 's/.*): //')" >> "${RES_DIR}/ALL.tsv"; }
printf "base\tstage\twall_time\tpeak_rss_kb\n" > "${RES_DIR}/ALL.tsv"

if [[ ${MOD} == sr ]]; then
    ts="${RES_DIR}/${SET}.interleave.time-v.stderr"; logger "${SET}: interleave"
    ${BINTIME} bash -c "source '${CONDA_SH}' && conda activate human-filtration && seqtk mergepe '${R1}' '${R2}' > '${PLAIN}'" \
        2> >(tee "${ts}" >&2)
    log_stage interleave "${ts}"
else
    ts="${RES_DIR}/${SET}.decompress.time-v.stderr"; logger "${SET}: decompress"
    ${BINTIME} bash -c "zcat '${R1}' > '${PLAIN}'" 2> >(tee "${ts}" >&2)
    log_stage decompress "${ts}"
fi

ts="${RES_DIR}/${SET}.movi2_filter.time-v.stderr"; logger "${SET}: movi query --filter -v"
LD_LIBRARY_PATH="${MOVI_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}" \
${BINTIME} bash -c 'out="$1"; shift; "$@" | gzip -c > "${out}"' _ "${CONCAT_SURV}" \
    "${MOVI_BIN}" query --filter -v --index "${MOVI_INDEX}" --read "${PLAIN}" -t "${THREADS}" \
    2> >(tee "${ts}" >&2)
log_stage movi2_filter "${ts}"

logger "${SET}: splitting survivors per source tag"
"${PY}" "${SPLIT_SURVIVORS}" "${CONCAT_SURV}" "${SURV_DIR}" ".survivors.fastq.gz"
logger "${SET}: done"

#!/bin/bash -l
#SBATCH -J arm-sota
#SBATCH --mem 64gb
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 48:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# State of the art: minimap2 -> minimap2 -> Movi 2, one call per stage over
# the whole concatenated set, each stage timed:
#   short read  mm2 -ax sr vs human-GRC-db.mmi (input interleaved with seqtk
#               mergepe), keep pairs with both mates unmapped (-f 12 -F 256)
#               -> mm2 -ax sr vs human-GCA-phix-db.mmi, same filter
#               (the two Qiita human-pangenome minimap2 databases)
#   long read   mm2 -ax map-hifi vs GRCh38_hifi.mmi, keep unmapped
#               (-f 4 -F 0x900) -> mm2 -ax map-hifi vs T2T_hifi.mmi, same
#   then        movi query --filter -v against hprc-yr1, output gzipped
#               inside the timed call (FASTA, despite the file name)
#
# Output: <mod>/<set>/sota/{resource/ALL.tsv, resource/*.time-v.stderr,
# survivors/<tag>.survivors.fastq.gz}
#
# Usage (from the repo root): sbatch arms/arm-sota.sbatch <sr|lr> <holdout|wol2>

source ./common.sh

MOD="${1:?usage: arm-sota.sbatch <sr|lr> <holdout|wol2>}"
SET="${2:?missing set}"
THREADS=${SLURM_CPUS_PER_TASK:-8}
case "${MOD}" in
    sr) R1="${CONCAT_SR}/${SET}_R1.fastq.gz"; R2="${CONCAT_SR}/${SET}_R2.fastq.gz"; MOD_DIR=short-read
        IDX1=${QIITA_HG38_MMI}; IDX2=${QIITA_T2T_MMI}; PRESET=sr; KEEP="-f 12 -F 256"
        STAGE1=align_hg38 ;;
    lr) R1="${CONCAT_LR}/${SET}.fastq.gz"; R2=; MOD_DIR=long-read
        IDX1=${LR_SOTA_HG38_MMI}; IDX2=${LR_RESCUE_MMI}; PRESET=map-hifi; KEEP="-f 4 -F 0x900"
        STAGE1=align_grch38 ;;
    *)  logger "unknown modality: ${MOD}"; exit 1 ;;
esac
case "${SET}" in holdout|wol2) ;; *) logger "unknown set: ${SET}"; exit 1 ;; esac
require_files "${R1}" ${R2:+"${R2}"} "${IDX1}" "${IDX2}" "${MOVI_BIN}" "${MOVI_INDEX}"

OUT_DIR="${RESULTS_DIR}/${MOD_DIR}/${SET}/sota"
RES_DIR="${OUT_DIR}/resource"; SURV_DIR="${OUT_DIR}/survivors"
mkdir -p "${RES_DIR}" "${SURV_DIR}"
job_scratch "arm-sota-${MOD}-${SET}"
STEP1_OUT="${TMPDIR}/${SET}.step1.fastq"; STEP2_OUT="${TMPDIR}/${SET}.step2.fastq"
CONCAT_SURV="${TMPDIR}/${SET}.survivors.fastq.gz"
ACTIVATE="source '${CONDA_SH}' && conda activate human-filtration"

log_stage () { printf "%s\t%s\t%s\t%s\n" ALL "$1" \
    "$(grep -m1 'Elapsed (wall clock)' "$2" | sed 's/.*): //')" \
    "$(grep -m1 'Maximum resident set size' "$2" | sed 's/.*): //')" >> "${RES_DIR}/ALL.tsv"; }
printf "base\tstage\twall_time\tpeak_rss_kb\n" > "${RES_DIR}/ALL.tsv"

ts="${RES_DIR}/${SET}.${STAGE1}.time-v.stderr"; logger "${SET}: minimap2 step 1 (${PRESET})"
if [[ ${MOD} == sr ]]; then
    ILV="${TMPDIR}/${SET}.interleaved.fastq"
    ${BINTIME} bash -c "${ACTIVATE} && seqtk mergepe '${R1}' '${R2}' > '${ILV}' && minimap2 -2 -ax sr -t ${THREADS} '${IDX1}' '${ILV}' | samtools fastq -@ ${THREADS} ${KEEP} > '${STEP1_OUT}'" \
        2> >(tee "${ts}" >&2)
else
    ${BINTIME} bash -c "${ACTIVATE} && minimap2 -2 -ax map-hifi -t ${THREADS} '${IDX1}' '${R1}' | samtools fastq -@ ${THREADS} ${KEEP} > '${STEP1_OUT}'" \
        2> >(tee "${ts}" >&2)
fi
log_stage "${STAGE1}" "${ts}"

ts="${RES_DIR}/${SET}.align_t2t.time-v.stderr"; logger "${SET}: minimap2 step 2 (${PRESET})"
${BINTIME} bash -c "${ACTIVATE} && minimap2 -2 -ax ${PRESET} -t ${THREADS} '${IDX2}' '${STEP1_OUT}' | samtools fastq -@ ${THREADS} ${KEEP} > '${STEP2_OUT}'" \
    2> >(tee "${ts}" >&2)
log_stage align_t2t "${ts}"

ts="${RES_DIR}/${SET}.movi2_filter.time-v.stderr"; logger "${SET}: movi query --filter -v"
LD_LIBRARY_PATH="${MOVI_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}" \
${BINTIME} bash -c 'out="$1"; shift; "$@" | gzip -c > "${out}"' _ "${CONCAT_SURV}" \
    "${MOVI_BIN}" query --filter -v --index "${MOVI_INDEX}" --read "${STEP2_OUT}" -t "${THREADS}" \
    2> >(tee "${ts}" >&2)
log_stage movi2_filter "${ts}"

logger "${SET}: splitting survivors per source tag"
"${PY}" "${SPLIT_SURVIVORS}" "${CONCAT_SURV}" "${SURV_DIR}" ".survivors.fastq.gz"
logger "${SET}: done"

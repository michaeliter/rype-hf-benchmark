#!/bin/bash -l
#SBATCH -J simulate-lr-holdout
#SBATCH --mem 64gb
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 48:00:00
#SBATCH --output logs/%x-%A_%a.out
#SBATCH --error  logs/%x-%A_%a.err
#
# Long-read holdout (FN) reads: PacBio HiFi via pbsim3 multi-pass subreads
# (ERRHMM Sequel model, 10 passes, 1x depth, length 15,000 +- 5,000) and ccs
# consensus, for each of the 10 holdout haplotypes in references/holdout10.tsv.
# pbsim3 writes one subread BAM per contig; every contig's ccs must succeed
# before the output is written. ccs temp files go to SCRATCH_DIR, not /tmp.
# The pbsim BAMs are cached so a rerun after a ccs failure skips pbsim.
#
# Output: ${READS_DIR}/long-read/holdout/raw/<haplotype>.fastq.gz
#
# Usage (from the repo root): sbatch --array=1-10 reads/simulate-lr-holdout.sbatch

source ./common.sh
conda_activate hifi-sim

MODEL="${CONDA_PREFIX}/data/ERRHMM-SEQUEL.model"
THREADS=${SLURM_CPUS_PER_TASK:-16}
FILE=$(awk -F'\t' -v n="$(( SLURM_ARRAY_TASK_ID + 1 ))" 'NR == n { print $2 }' "${REPO}/references/holdout10.tsv")
REF="${HOLDOUT_REF_DIR}/${FILE}"
require_files "${REF}" "${MODEL}"
base="${FILE%.fa.gz}"
OUT_DIR="${READS_DIR}/long-read/holdout/raw"
OUT="${OUT_DIR}/${base}.fastq.gz"
mkdir -p "${OUT_DIR}"
[[ -s ${OUT} ]] && { logger "${OUT} exists, skipping"; exit 0; }

WORK="${SCRATCH_DIR}/simulate-lr-holdout/${base}"     # ccs workspace for this run
SAVED="${SCRATCH_DIR}/simulate-lr-holdout/saved/${base}"   # pbsim BAM cache
rm -rf "${WORK}"; mkdir -p "${WORK}/cctmp"
export TMPDIR="${WORK}/cctmp"

shopt -s nullglob
if compgen -G "${SAVED}/sd_*.bam" > /dev/null; then
    logger "${base}: reusing cached pbsim BAMs"
else
    gunzip -c "${REF}" > "${WORK}/genome.fa"
    logger "${base}: pbsim3 (errhmm SEQUEL, 10 passes, 1x)"
    ( cd "${WORK}" && pbsim --strategy wgs --method errhmm --errhmm "${MODEL}" \
        --pass-num 10 --depth 1 --length-mean 15000 --length-sd 5000 \
        --genome genome.fa --prefix sd --id-prefix "${base}" )
    bams=( "${WORK}"/sd_*.bam )
    [[ ${#bams[@]} -gt 0 ]] || { logger "pbsim produced no BAM"; exit 1; }
    mkdir -p "${SAVED}" && mv "${WORK}"/sd_*.bam "${SAVED}/"
fi

bams=( "${SAVED}"/sd_*.bam )
logger "${base}: ccs over ${#bams[@]} per-contig subread BAMs"
: > "${WORK}/all.fastq"
ok=0
for bam in "${bams[@]}"; do
    out="${WORK}/$(basename "${bam%.bam}").ccs.fastq.gz"
    if ccs --num-threads "${THREADS}" "${bam}" "${out}" 2>> "${WORK}/ccs.log"; then
        ok=$((ok + 1))
        [[ -s ${out} ]] && zcat "${out}" >> "${WORK}/all.fastq"   # a contig with no consensus is fine
    else
        logger "ccs failed on $(basename "${bam}")"
    fi
done
[[ ${ok} -eq ${#bams[@]} ]] || { logger "${base}: ccs finished ${ok}/${#bams[@]} contigs; resubmit to retry ccs"; exit 1; }

gzip -c "${WORK}/all.fastq" > "${OUT}"
logger "${base}: done -> ${OUT} ($(( $(zcat "${OUT}" | wc -l) / 4 )) HiFi reads)"
rm -rf "${WORK}" "${SAVED}"

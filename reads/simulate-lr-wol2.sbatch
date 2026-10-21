#!/bin/bash -l
#SBATCH -J simulate-lr-wol2
#SBATCH --mem 32gb
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%A_%a.out
#SBATCH --error  logs/%x-%A_%a.err
#
# Long-read WoL2 (FP) reads: pbsim3 + ccs per WoL2 genome with the same
# settings as the holdout (simulate-lr-holdout.sbatch), 100 genomes per array
# task (reads/wol2-genomes.list: 15,953 genomes, 160 tasks). A genome is
# written only if ccs succeeded on every contig; others are skipped with a
# warning and retried on resubmit. 5 genomes produced no reads in the
# manuscript set (15,948 of 15,953; reads/wol2-chunks/long-read/).
#
# Output: ${READS_DIR}/long-read/wol2/sim/<accession>.fastq.gz
#
# Usage (from the repo root): sbatch --array=1-160 reads/simulate-lr-wol2.sbatch

source ./common.sh
set +e    # per-genome failures are reported and skipped
conda_activate hifi-sim

MODEL="${CONDA_PREFIX}/data/ERRHMM-SEQUEL.model"
THREADS=${SLURM_CPUS_PER_TASK:-1}
OUT_DIR="${READS_DIR}/long-read/wol2/sim"
mkdir -p "${OUT_DIR}"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * 100 + 1 ))
job_scratch simulate-lr-wol2

while read -r acc; do
    OUT="${OUT_DIR}/${acc}.fastq.gz"
    [[ -s ${OUT} ]] && continue
    GZ="${WOL2_GENOME_DIR}/${acc}.fasta.gz"
    [[ -e ${GZ} ]] || { logger "missing genome ${GZ}"; continue; }
    W="${TMPDIR}/${acc}"; mkdir -p "${W}/cctmp"
    gunzip -c "${GZ}" > "${W}/genome.fa"
    ( cd "${W}" && pbsim --strategy wgs --method errhmm --errhmm "${MODEL}" \
        --pass-num 10 --depth 1 --length-mean 15000 --length-sd 5000 \
        --genome genome.fa --prefix sd --id-prefix "${acc}" ) > /dev/null 2>> "${W}/pbsim.log" \
        || { logger "pbsim failed: ${acc}"; rm -rf "${W}"; continue; }
    : > "${W}/all.fastq"
    shopt -s nullglob
    bams=( "${W}"/sd_*.bam ); ok=0
    for bam in "${bams[@]}"; do
        if TMPDIR="${W}/cctmp" ccs --num-threads "${THREADS}" "${bam}" "${bam%.bam}.ccs.fastq.gz" 2>> "${W}/ccs.log"; then
            ok=$((ok + 1))
            [[ -s ${bam%.bam}.ccs.fastq.gz ]] && zcat "${bam%.bam}.ccs.fastq.gz" >> "${W}/all.fastq"
        fi
    done
    shopt -u nullglob
    if [[ ${#bams[@]} -gt 0 && ${ok} -eq ${#bams[@]} && -s ${W}/all.fastq ]]; then
        gzip -c "${W}/all.fastq" > "${OUT}"
    else
        logger "${acc}: ccs ${ok}/${#bams[@]} contigs (or no reads), skipped"
    fi
    rm -rf "${W}"
done < <(sed -n "${start},$(( start + 99 ))p" "${REPO}/reads/wol2-genomes.list")
logger "task ${SLURM_ARRAY_TASK_ID}: done"

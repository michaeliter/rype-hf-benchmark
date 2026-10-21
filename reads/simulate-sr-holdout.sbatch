#!/bin/bash -l
#SBATCH -J simulate-sr-holdout
#SBATCH --time 24:00:00
#SBATCH --mem 16gb
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --output logs/%x-%A_%a.out
#SBATCH --error  logs/%x-%A_%a.err
#
# Short-read holdout (FN) reads: 1x-coverage 150 bp paired-end Illumina reads
# with ART (HiSeq 2500 profile, 450 +- 50 bp fragments) for each of the 10
# holdout haplotypes in references/holdout10.tsv. ART seeds from the clock;
# the seed is printed to the job log.
#
# Usage (from the repo root): sbatch --array=1-10 reads/simulate-sr-holdout.sbatch

source ./common.sh
conda_activate art

FILE=$(awk -F'\t' -v n="$(( SLURM_ARRAY_TASK_ID + 1 ))" 'NR == n { print $2 }' "${REPO}/references/holdout10.tsv")
REF="${HOLDOUT_REF_DIR}/${FILE}"
require_files "${REF}"
base="${FILE%.fa.gz}"
OUT_DIR="${READS_DIR}/short-read/holdout/raw"
R1="${OUT_DIR}/${base}_R1.fastq.gz"
R2="${OUT_DIR}/${base}_R2.fastq.gz"
mkdir -p "${OUT_DIR}"
[[ -s ${R1} && -s ${R2} ]] && { logger "${base} already simulated, skipping"; exit 0; }

job_scratch simulate-sr-holdout
gunzip -c "${REF}" > "${TMPDIR}/${base}.fa"
logger "ART HS25 1x PE150 on ${base}"
art_illumina -ss HS25 -i "${TMPDIR}/${base}.fa" -p -l 150 -f 1 -m 450 -s 50 \
    -o "${TMPDIR}/${base}_" -na
gzip -c "${TMPDIR}/${base}_1.fq" > "${R1}"
gzip -c "${TMPDIR}/${base}_2.fq" > "${R2}"
logger "done -> ${R1} ${R2}"

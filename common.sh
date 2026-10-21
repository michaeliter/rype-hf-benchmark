#!/bin/bash
# Sourced by every script: site paths (config.sh), the derived layout of
# indices / reads / results, the operating points and small helpers.
#
# Run everything from the repository root: the run-*.sh submitters cd
# there, and a batch job starts in the directory it was submitted from.

export REPO=$(pwd)
[[ -f ${REPO}/config.sh ]] || { echo "FATAL: run from the repository root (no ${REPO}/config.sh)" >&2; exit 1; }
source "${REPO}/config.sh"
export TREE_RSS=${REPO}/eval/tree_rss.py
export SPLIT_SURVIVORS=${REPO}/eval/split_survivors.py

# --- indices -------------------------------------------------------------------
# rype: one bucket over the 97 references, k=64, one index per window w.
rype_index () { echo "${INDEX_DIR}/rype/human-all97-w${1}.ryxdi"; }
export SR_RESCUE_MMI=${INDEX_DIR}/minimap2/T2T_sr.mmi       # T2T-CHM13v2.0, sr preset
export LR_RESCUE_MMI=${INDEX_DIR}/minimap2/T2T_hifi.mmi     # T2T-CHM13v2.0, map-hifi
export LR_SOTA_HG38_MMI=${INDEX_DIR}/minimap2/GRCh38_hifi.mmi
export CLEANIFIER_EXACT_IDX=${INDEX_DIR}/cleanifier/human-all97-exact-k29w33-n3.2e9    # .hash + .info
export CLEANIFIER_FILTER_IDX=${INDEX_DIR}/cleanifier/human-all97-filter-k29w33-fpr14   # .filter + .info

# --- benchmark inputs (reads/build-concat-inputs.sbatch) ---------------------------
# One ID-prefixed file (pair) per modality and set: every header is
# "@<tag>::<id>", tag = holdout haplotype or WoL2 chunk.
export CONCAT_SR=${CONCAT_DIR}/short-read    # {holdout,wol2}_R{1,2}.fastq.gz
export CONCAT_LR=${CONCAT_DIR}/long-read     # {holdout,wol2}.fastq.gz

# --- operating points ----------------------------------------------------------
# w=20 for both modalities; short read thresholded at 0.0001 with the minimap2
# rescue pass, long read at 0.05 with rype alone.
export OP_WINDOW=20
export SR_THRESHOLD=0.0001
export LR_THRESHOLD=0.05
# rype-hf sizes rype's memory budget from the job limit (rype-hf README).
export RYPE_MAX_MEMORY=auto

# cleanifier filter call, per its docs: default threshold 0.5 and sampling mode,
# gz output at cleanifier's defaults (bgzip -l1, 2 threads per file). Only long
# reads raise --buffersize (to 22): the default 2^16-byte buffer garbles
# ~15 kb HiFi reads, and the docs recommend 20-24 for long reads.
# Sets CLEANIFIER_ARGS.
# Usage: cleanifier_args <sr|lr> <index> <r1> <r2|""> <out_prefix> <threads>
function cleanifier_args () {
    CLEANIFIER_ARGS=(filter --index "$2" --fastq "$3")
    [[ -n $4 ]] && CLEANIFIER_ARGS+=(--pairs "$4")
    CLEANIFIER_ARGS+=(--out "$5" --threshold 0.5 --threads "$6" --compression gz)
    [[ $1 == lr ]] && CLEANIFIER_ARGS+=(--buffersize 22)
    return 0
}

# Activate one of envs/*.yml by name (cleanifier must run inside its env:
# called by path it exits 0 and writes nothing).
function conda_activate () {
    set +u
    source "${CONDA_SH}"
    conda activate "$1"
    set -u
}

function logger () {
    echo "$(date) :: ${*}" 1>&2
}

function require_files () {
    local f
    for f in "$@"; do
        [[ -e ${f} ]] || { logger "FATAL: missing ${f}"; exit 1; }
    done
}

# Per-job scratch directory under SCRATCH_DIR, removed on exit.
function job_scratch () {
    export TMPDIR="${SCRATCH_DIR}/${1}/${SLURM_JOB_ID:-$$}"
    mkdir -p "${TMPDIR}"
    trap 'rm -rf "${TMPDIR}"' EXIT
}

set -euo pipefail

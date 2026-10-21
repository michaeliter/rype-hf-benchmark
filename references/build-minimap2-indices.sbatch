#!/bin/bash -l
#SBATCH -J build-minimap2-indices
#SBATCH --mem 64gb
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 04:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# minimap2 indices built here (the two Qiita short-read SOTA indices are
# downloaded, see README):
#   T2T_sr.mmi       T2T-CHM13v2.0, sr preset        rype -> mm2, short read
#   T2T_hifi.mmi     T2T-CHM13v2.0, map-hifi preset  rype -> mm2, long read;
#                                                    SOTA long read, step 2
#   GRCh38_hifi.mmi  GRCh38.p14, map-hifi preset     SOTA long read, step 1
#
# The indices used in the manuscript were written by miint's
# save_minimap2_index (DuckDB) with the same presets; this is the equivalent
# minimap2 CLI call.
#
# Usage (from the repo root): sbatch references/build-minimap2-indices.sbatch

source ./common.sh

THREADS=${SLURM_CPUS_PER_TASK:-16}
T2T="${REF_DIR}/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna"
HG38="${REF_DIR}/GCF_000001405.40_GRCh38.p14_genomic.fna"
require_files "${MINIMAP2_BIN}" "${T2T}" "${HG38}"
mkdir -p "$(dirname "${SR_RESCUE_MMI}")"

build_one () {  # <preset> <reference> <out.mmi>
    [[ -s $3 ]] && { logger "$3 exists, skipping"; return 0; }
    logger "minimap2 -x $1 -d $3 $2"
    "${MINIMAP2_BIN}" -x "$1" -t "${THREADS}" -d "$3.tmp" "$2"
    mv "$3.tmp" "$3"
}
build_one sr       "${T2T}"  "${SR_RESCUE_MMI}"
build_one map-hifi "${T2T}"  "${LR_RESCUE_MMI}"
build_one map-hifi "${HG38}" "${LR_SOTA_HG38_MMI}"
logger "done"

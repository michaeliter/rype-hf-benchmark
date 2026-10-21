#!/bin/bash -l
#SBATCH -J build-rype-index
#SBATCH --mem 128gb
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 24:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# rype host index for one window w: one bucket over the 97 references of
# human-all97.tsv, k=64, orient on. Writes the TOML config from the reference
# list, then `rype index from-config` into ${INDEX_DIR}/rype/.
#
# The benchmark uses w in {10, 20, 25, 35, 45} (short read) and
# {20, 50, 100, 150, 200} (long read); w=20 is the operating point.
# (The w=10 index was built with rype's default --row-group-size and no
# --max-memory; both change how the build is batched, not the index content.)
#
# Usage (from the repo root): sbatch references/build-rype-index.sbatch <w>

source ./common.sh

W="${1:?usage: build-rype-index.sbatch <window>}"
NAME="human-all97-w${W}"
OUT_DIR="${INDEX_DIR}/rype"
[[ -d ${OUT_DIR}/${NAME}.ryxdi ]] && { logger "${OUT_DIR}/${NAME}.ryxdi exists, skipping"; exit 0; }
mkdir -p "${OUT_DIR}"

{
    printf '[index]\nk = 64\nwindow = %s\nsalt = 0x5555555555555555\noutput = "%s"\norient_sequences = true\n\n' "${W}" "${NAME}"
    printf '[buckets.%s]\nfiles = [\n' "${NAME}"
    awk -F'\t' -v dir="${REF_DIR}" 'NR > 1 { printf "    \"%s/%s\",\n", dir, $1 }' "${REPO}/references/human-all97.tsv"
    printf ']\n'
} > "${OUT_DIR}/${NAME}.toml"
awk -F'\t' -v dir="${REF_DIR}" 'NR > 1 { print dir "/" $1 }' "${REPO}/references/human-all97.tsv" \
    | while read -r f; do require_files "${f}"; done

# rype resolves [index].output against the working directory.
cd "${OUT_DIR}"
logger "rype index from-config: ${NAME}"
"${RYPE_BIN}" index from-config \
    --config "${NAME}.toml" \
    --orient \
    --timing \
    --max-memory 100GB \
    --row-group-size 1000000 \
    --verbose
logger "done: ${OUT_DIR}/${NAME}.ryxdi"

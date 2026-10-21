#!/bin/bash -l
#SBATCH -J build-cleanifier-index
#SBATCH --mem 40gb
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# cleanifier indices over the same 97 references as rype, built the way the
# cleanifier paper and docs recommend (Zentgraf et al., Bioinformatics 2026,
# btaf632):
#   exact    (3,4) bucketed Cuckoo hash table, default fill 0.85     (~15.5 GB)
#   filter   (2,2) windowed Cuckoo filter, FPR 2^-14, 95% load       (~6.7 GB)
# Both: mask k29 w33 and 15 subtables, as in the authors' prebuilt index;
# capacity -n 3.2e9 (the references hold 3,148,289,270 distinct gapped
# k-mers); --maxwalk 10000, as in the authors' index, and no --maxfailures
# -1, so a failed insertion fails the build instead of dropping k-mers.
#
# Usage (from the repo root): sbatch references/build-cleanifier-index.sbatch <exact|filter>

source ./common.sh

KIND="${1:?usage: build-cleanifier-index.sbatch <exact|filter>}"
case "${KIND}" in
    exact)  OUT=${CLEANIFIER_EXACT_IDX};  EXT=hash
            TABLE=(--exact --bucketsize 4 --fill 0.85) ;;
    filter) OUT=${CLEANIFIER_FILTER_IDX}; EXT=filter
            TABLE=(--fpr 14 --windowsize 2 --fill 0.95) ;;
    *)      logger "unknown index kind: ${KIND}"; exit 1 ;;
esac
[[ -s ${OUT}.${EXT} ]] && { logger "index exists, refusing to overwrite: ${OUT}.${EXT}"; exit 0; }
mkdir -p "$(dirname "${OUT}")"
mapfile -t REFS < <(awk -F'\t' -v dir="${REF_DIR}" 'NR > 1 { print dir "/" $1 }' "${REPO}/references/human-all97.tsv")
[[ ${#REFS[@]} -eq 97 ]] || { logger "expected 97 reference files, got ${#REFS[@]}"; exit 1; }
require_files "${REFS[@]}"

conda_activate cleanifier
logger "cleanifier $(cleanifier --version 2>&1 | head -1) on $(hostname): ${KIND}"
${BINTIME} -o "${OUT}.build-time-v.stderr" cleanifier index \
    --index "${OUT}" \
    --files "${REFS[@]}" \
    --mask "######_#######_###_#######_######" \
    -n 3200000000 \
    "${TABLE[@]}" \
    --subtables 15 \
    --maxwalk 10000 \
    --statistics summary
ls -la "${OUT}".* >&2
logger "done: ${OUT}"

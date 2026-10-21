#!/bin/bash
# Threshold sweep (threshold supplementary figure): rype alone at w=20 over
# 11 thresholds. One 8 CPU / 64 GB job per modality and set scores every read
# once (`rype classify run` at t=0.0001) and counts host reads per threshold.
# -> 4 jobs, then the eval and figures (afterok).
#
# Usage: ./run-threshold-sweep.sh
cd "$(dirname "$0")" && source common.sh
mkdir -p logs

JOBS=()
for MOD in sr lr; do
    for SET in holdout wol2; do
        jid=$(sbatch --parsable -J "threshold-sweep-${MOD}-${SET}" arms/arm-threshold-sweep.sbatch "${MOD}" "${SET}")
        JOBS+=("${jid}"); echo "${jid}  ${MOD} ${SET}"
    done
done
jid=$(sbatch --parsable --dependency="afterok:$(IFS=:; echo "${JOBS[*]}")" eval/eval.sbatch)
echo "eval + figures: ${jid} (after ${#JOBS[@]} jobs)"

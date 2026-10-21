#!/bin/bash
# Main benchmark (Fig. 2 and the window-independent rows of every figure):
# every tool x {short, long read} x {holdout, wol2}, each task its own
# 8 CPU / 64 GB job, then the eval and figures (afterok):
#   rype, rype-mm2                arms/arm-rype-stream.sbatch (REPS runs, default 1)
#   cleanifier exact, filter      arms/arm-cleanifier.sbatch
#   Movi 2 alone                  arms/arm-movi2-only.sbatch
#   mm2 -> mm2 -> Movi 2 (SOTA)   arms/arm-sota.sbatch
#
# Usage: [REPS=1] ./run-main.sh
cd "$(dirname "$0")" && source common.sh
mkdir -p logs
REPS=${REPS:-1}

JOBS=()
submit () { local jid; jid=$(sbatch --parsable "$@"); JOBS+=("${jid}"); echo "${jid}  ${*: -4}"; }
for MOD in sr lr; do
    for SET in holdout wol2; do
        for MODE in rype rype-mm2; do
            for REP in $(seq 1 "${REPS}"); do
                submit -J "${MODE}-${MOD}-${SET}-r${REP}" arms/arm-rype-stream.sbatch "${MOD}" "${SET}" "${MODE}" "${REP}"
            done
        done
        for KIND in exact filter; do
            submit -J "cleanifier-${KIND}-${MOD}-${SET}" arms/arm-cleanifier.sbatch "${KIND}" "${MOD}" "${SET}"
        done
        submit -J "movi2-only-${MOD}-${SET}" arms/arm-movi2-only.sbatch "${MOD}" "${SET}"
        submit -J "sota-${MOD}-${SET}" arms/arm-sota.sbatch "${MOD}" "${SET}"
    done
done
jid=$(sbatch --parsable --dependency="afterok:$(IFS=:; echo "${JOBS[*]}")" eval/eval.sbatch)
echo "eval + figures: ${jid} (after ${#JOBS[@]} jobs)"

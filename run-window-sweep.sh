#!/bin/bash
# Window sweep (window supplementary figure): rype and rype-mm2 at every
# other window, one 8 CPU / 64 GB job each; w=20 is the run-main.sh run.
#   short read  w in {10, 25, 35, 45}     x {rype, rype-mm2} x {holdout, wol2}
#   long read   w in {50, 100, 150, 200}  x {rype, rype-mm2} x {holdout, wol2}
# -> 32 jobs, then the eval and figures (afterok).
#
# Usage: ./run-window-sweep.sh
cd "$(dirname "$0")" && source common.sh
mkdir -p logs

declare -A WINDOWS=([sr]="10 25 35 45" [lr]="50 100 150 200")
JOBS=()
for MOD in sr lr; do
    for W in ${WINDOWS[${MOD}]}; do
        for SET in holdout wol2; do
            for MODE in rype rype-mm2; do
                jid=$(sbatch --parsable -J "sweep-${MODE}-${MOD}-${SET}-w${W}" \
                    arms/arm-rype-stream.sbatch "${MOD}" "${SET}" "${MODE}" 1 "${W}")
                JOBS+=("${jid}"); echo "${jid}  ${MODE} ${MOD} ${SET} w${W}"
            done
        done
    done
done
jid=$(sbatch --parsable --dependency="afterok:$(IFS=:; echo "${JOBS[*]}")" eval/eval.sbatch)
echo "eval + figures: ${jid} (after ${#JOBS[@]} jobs)"

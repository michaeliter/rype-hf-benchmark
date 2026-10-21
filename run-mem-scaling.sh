#!/bin/bash
# Memory scaling (memory-scaling supplementary figure): long-read rype at
# 8 CPU and job memory limits 64, 32, 16, 8, 4, 2, 1 GB x {holdout, wol2}
# x 3 replicates = 42 jobs, then the eval and figures (afterany, so a run
# that OOMs still shows up, as a failed point). Run after run-main.sh: every
# run's counts are checked against the 64 GB operating-point run.
#
# Usage: [TIERS="64 32 16 8 4 2 1"] [REPS=3] ./run-mem-scaling.sh
cd "$(dirname "$0")" && source common.sh
mkdir -p logs
TIERS=${TIERS:-"64 32 16 8 4 2 1"}
REPS=${REPS:-3}

JOBS=()
for MEM_GB in ${TIERS}; do
    for SET in holdout wol2; do
        for REP in $(seq 1 "${REPS}"); do
            jid=$(sbatch --parsable --mem "${MEM_GB}gb" -J "mem-scaling-${SET}-${MEM_GB}gb-r${REP}" \
                arms/arm-mem-scaling.sbatch "${SET}" "${REP}")
            JOBS+=("${jid}"); echo "${jid}  ${SET} ${MEM_GB} GB rep${REP}"
        done
    done
done
jid=$(sbatch --parsable --dependency="afterany:$(IFS=:; echo "${JOBS[*]}")" eval/eval.sbatch)
echo "eval + figures: ${jid} (after ${#JOBS[@]} jobs)"

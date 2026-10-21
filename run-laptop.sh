#!/bin/bash
# Laptop tiers (laptop supplementary figures): 8 CPU and 16 or 8 GB, each
# task its own job, then the eval and figures (afterany, so an OOM still
# shows up as a failed row). Run after run-main.sh: the tables compare with
# the 64 GB runs.
#   16 GB  {rype, rype-mm2, cleanifier, cleanifier-filter} x {sr, lr} x
#          {holdout, wol2}, except long-read rype-mm2 (minimap2's map-hifi
#          index alone peaks at 13.8 GB) -> 14 jobs
#   8 GB   long read only: {rype, cleanifier-filter} x {holdout, wol2} -> 4 jobs
#
# Usage: ./run-laptop.sh [16|8]      (default 16)
cd "$(dirname "$0")" && source common.sh
mkdir -p logs
MEM_GB=${1:-16}
case "${MEM_GB}" in
    16) MODS=(sr lr); TOOLS=(rype rype-mm2 cleanifier cleanifier-filter) ;;
    8)  MODS=(lr);    TOOLS=(rype cleanifier-filter) ;;
    *)  echo "usage: run-laptop.sh [16|8]" >&2; exit 2 ;;
esac

JOBS=()
for MOD in "${MODS[@]}"; do
    for SET in holdout wol2; do
        for TOOL in "${TOOLS[@]}"; do
            [[ ${MOD} == lr && ${TOOL} == rype-mm2 ]] && continue
            jid=$(sbatch --parsable --mem "${MEM_GB}gb" -J "laptop${MEM_GB}-${TOOL}-${MOD}-${SET}" \
                arms/arm-laptop.sbatch "${MOD}" "${SET}" "${TOOL}")
            JOBS+=("${jid}"); echo "${jid}  laptop-${MEM_GB}gb ${TOOL} ${MOD} ${SET}"
        done
    done
done
jid=$(sbatch --parsable --dependency="afterany:$(IFS=:; echo "${JOBS[*]}")" eval/eval.sbatch)
echo "eval + figures: ${jid} (after ${#JOBS[@]} jobs)"

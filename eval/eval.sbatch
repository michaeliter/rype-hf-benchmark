#!/bin/bash -l
#SBATCH -J eval
#SBATCH --mem 16gb
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH --output logs/%x-%j.out
#SBATCH --error  logs/%x-%j.err
#
# Builds every figure table in data/ from the finished runs, then draws the
# figures. Counting only: no classification. Record counts run in a thread
# pool sized to the job's CPUs.
#
#   eval_external.py   read totals, Movi 2 alone and SOTA   -> ${RESULTS_DIR}/<mod>/tables/
#   eval_op.py         + rype, rype-mm2, cleanifier         -> data/<mod>/
#   eval_laptop.py     16 GB and 8 GB tiers                 -> data/<mod>/laptop-<N>gb_runtime.tsv
#   eval_mem_scaling.py                                     -> data/long-read/mem_scaling.tsv
#   figures/make-figures.sh                                 -> figures/output/
# A stage whose runs are incomplete warns and leaves those rows out.
#
# Usage (from the repo root): sbatch eval/eval.sbatch

source ./common.sh
cd "${REPO}"
"${PY}" eval/eval_external.py sr lr
"${PY}" eval/eval_op.py sr lr
"${PY}" eval/eval_laptop.py --mem 16
"${PY}" eval/eval_laptop.py --mem 8
"${PY}" eval/eval_mem_scaling.py
PY="${PY}" figures/make-figures.sh
logger "done"

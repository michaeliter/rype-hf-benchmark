#!/bin/bash
# Draws every manuscript figure from the tables in data/ into figures/output/.
# Needs only the hf-repro environment (envs/hf-repro.yml), not the cluster:
#   PY=/path/to/envs/hf-repro/bin/python3 figures/make-figures.sh
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
"${PY}" fig2_benchmarking.py          # main benchmarking figure
"${PY}" ed_fig1_window.py             # window sweep
"${PY}" ed_fig2_threshold.py          # threshold sweep
"${PY}" supp_laptop.py --mem 16       # laptop tier, 16 GB
"${PY}" supp_laptop.py --mem 8        # laptop tier, 8 GB (long read)
"${PY}" supp_mem_scaling.py           # memory scaling (long read)

#!/bin/bash
# Site configuration: every path and tool location the pipeline uses. Edit the
# values below for your system before running anything; every script sources
# this file (through common.sh). Paths may be absolute or relative to the repo.

# --- working space (large: ~100 GB of reads, ~550 GB of results) -------------
export WORK_DIR=${WORK_DIR:-${REPO}/work}
export READS_DIR=${WORK_DIR}/reads           # simulated + QC'd reads
export CONCAT_DIR=${WORK_DIR}/concat         # concatenated benchmark inputs
export INDEX_DIR=${WORK_DIR}/indices         # every host index
export RESULTS_DIR=${WORK_DIR}/results       # per-run outputs
export SCRATCH_DIR=${WORK_DIR}/tmp           # per-job scratch (fast disk)

# --- inputs you provide (see README, "Inputs") --------------------------------
export REF_DIR=/path/to/human-refs               # the 97 files of references/human-all97.tsv
export HOLDOUT_REF_DIR=/path/to/hprc-r2          # the 10 assemblies of references/holdout10.tsv (*.fa.gz)
export WOL2_GENOME_DIR=/path/to/wol2-genomes     # <accession>.fasta.gz, long-read simulation
export WOL2_SR_READS_DIR=/path/to/wol2-illumina  # <accession>_R{1,2}.fastq.gz, short-read FP set
export ADAPTERS=/path/to/known_adapters.fna      # fastp adapter list (README, "Inputs")
export MOVI_INDEX=/path/to/movi/hprc-yr1         # Movi 2 index (README, "Indices")
export QIITA_HG38_MMI=/path/to/human-GRC-db.mmi        # SOTA short read, step 1
export QIITA_T2T_MMI=/path/to/human-GCA-phix-db.mmi    # SOTA short read, step 2

# --- tools not installed through conda ---------------------------------------
export RYPE_BIN=/path/to/RYpe/target/release/rype   # michaeliter/RYpe @ f2ed8b8
export RYPE_HF=/path/to/rype-hf-cli/bin/rype-hf     # michaeliter/rype-hf-cli @ 4bbdae7
export MOVI_BIN=/path/to/Movi/build/movi            # mohsenzakeri/Movi @ 9bc0acc
export BINTIME="/usr/bin/time -v"                   # GNU time

# --- conda environments (envs/*.yml) ------------------------------------------
export CONDA_SH=/path/to/miniforge3/etc/profile.d/conda.sh
export CONDA_ENVS=/path/to/miniforge3/envs
export PY=${CONDA_ENVS}/hf-repro/bin/python3
export DUCKDB_BIN=${CONDA_ENVS}/hf-repro/bin/duckdb   # README, "Environments"
export MINIMAP2_BIN=${CONDA_ENVS}/minimap2/bin/minimap2
export SAMTOOLS_BIN=${CONDA_ENVS}/minimap2/bin/samtools
export BGZIP_BIN=${CONDA_ENVS}/minimap2/bin/bgzip
export PIGZ_BIN=${CONDA_ENVS}/pigz/bin/pigz
# Movi links against libgomp; the hf-repro env provides it.
export MOVI_LD_LIBRARY_PATH=${CONDA_ENVS}/hf-repro/lib

# --- Slurm -------------------------------------------------------------------
# Every job asks for 8 CPUs and its own memory (64 GB for the benchmark arms).
# sbatch reads SBATCH_PARTITION itself; uncomment to pick a partition.
# export SBATCH_PARTITION=short

# rype-hf-benchmark

Scripts, tables and figures for the human host-filtration benchmark of
[rype](https://github.com/the-miint/rype) through
[rype-hf](https://github.com/michaeliter/rype-hf-cli), its streaming
host-filtration pipeline. The benchmark compares it with cleanifier (exact and
probabilistic), Movi 2 and minimap2 → minimap2 → Movi 2 (the state of the
art, SOTA) on simulated short (Illumina) and long (PacBio HiFi) reads:

- **holdout** reads come from 10 human haplotypes that are not in any index.
  Every read is host, so every survivor is a false negative (FN).
- **wol2** reads come from 15,953 microbial genomes of the Web of Life 2
  (WoL2). No read is host, so every removed read is a false positive (FP).

Every tool runs as one 8 CPU / 64 GB Slurm job per dataset, on the same
inputs.

## Repository layout

| Path | Contents |
|---|---|
| `config.sh` | every site path; edit before running |
| `common.sh` | derived paths, operating points, shared helpers |
| `references/` | reference lists and index builds |
| `reads/` | read simulation, QC, WoL2 chunk lists, benchmark inputs |
| `arms/` | one batch script per tool |
| `run-*.sh` | submit the benchmark runs |
| `eval/` | turn run outputs into the tables in `data/` |
| `figures/` | draw the figures from `data/` into `figures/output/` |
| `data/` | the tables behind every figure |
| `envs/` | conda environments and `versions.tsv` |

## Figures

| Figure | Script | Tables (`data/<modality>/`) | Runs |
|---|---|---|---|
| Benchmark at the operating points (main figure) | `fig2_benchmarking.py` | `op_runtime`, `op_fn`, `op_fp`, `storage` | `run-main.sh` |
| FN, FP and runtime vs rype window size | `ed_fig1_window.py` | `fn`, `fp`, `runtime`, `memory` | `run-main.sh`, `run-window-sweep.sh` |
| FN and FP vs rype threshold | `ed_fig2_threshold.py` | `threshold_sweep_w20`, `movi_ref` | `run-threshold-sweep.sh` |
| Laptop tiers, 16 GB and 8 GB | `supp_laptop.py --mem 16`, `--mem 8` | `laptop-16gb_runtime`, `laptop-8gb_runtime` | `run-laptop.sh 16`, `run-laptop.sh 8` |
| Long-read rype memory scaling, 64 GB to 1 GB | `supp_mem_scaling.py` | `long-read/mem_scaling` | `run-mem-scaling.sh` |

To redraw every figure from the committed tables, without a cluster:

```bash
conda env create -f envs/hf-repro.yml
PY=$(conda run -n hf-repro which python3) figures/make-figures.sh
```

## Operating points

| | rype window | rype threshold | minimap2 rescue pass |
|---|---|---|---|
| Short read | w = 20 | 0.0001 | yes: T2T-CHM13v2.0, `sr` preset |
| Long read | w = 20 | 0.05 | no (rype alone) |

The main figure's rype bar is rype → mm2 for short reads and rype alone for
long reads. The window and threshold sweeps show both rype modes.

## Software

Exact versions are in [`envs/versions.tsv`](envs/versions.tsv).
- **Conda environments.** Each environment in `envs/` was exported twice: as
  `<name>.yml` (`conda env export --no-builds`) and as `<name>.explicit.txt`
  (`conda list --explicit`). Create them under their own names:
  ```bash
  for f in envs/*.yml; do conda env create -f "$f"; done
  ```
- **hf-repro.** Its pip packages (matplotlib, numpy, pandas) are listed only
  in the `.yml`. DuckDB 1.5.4 is the official release binary, copied into the
  environment:
  ```bash
  conda activate hf-repro
  curl -sL https://github.com/duckdb/duckdb/releases/download/v1.5.4/duckdb_cli-linux-amd64.gz \
      | gunzip > "${CONDA_PREFIX}/bin/duckdb"
  chmod +x "${CONDA_PREFIX}/bin/duckdb" && duckdb --version    # v1.5.4 (Variegata)
  ```
- **rype** 1.0.0-rc.2, commit `f2ed8b8`, branch `michaeliter/classify-filter`
  of [michaeliter/RYpe](https://github.com/michaeliter/RYpe). This branch adds
  `rype classify filter`. It was built with rust 1.93.1 (rustup), using the
  `rype-build` environment's gcc 14.3 as C compiler and linker:
  ```bash
  git clone -b michaeliter/classify-filter https://github.com/michaeliter/RYpe.git
  cd RYpe && git checkout f2ed8b8
  conda activate rype-build
  CC=x86_64-conda-linux-gnu-cc CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=x86_64-conda-linux-gnu-cc \
      cargo build --release                              # target/release/rype
  ```
- **rype-hf** 0.1.0, commit `4bbdae7` of
  [michaeliter/rype-hf-cli](https://github.com/michaeliter/rype-hf-cli). It is
  a bash script and needs no build.
- **Movi 2**, commit `9bc0acc` of
  [mohsenzakeri/Movi](https://github.com/mohsenzakeri/Movi). Build it with
  cmake as in its README. It needs `libgomp` at run time; `config.sh` takes it
  from the hf-repro environment.
- **GNU time** (`/usr/bin/time -v`) and Slurm.

## Reproducing the benchmark

Edit `config.sh` first. Run every command from the repository root. Each
batch script's header documents its arguments.

### 1. Inputs

| Input | Source | `config.sh` |
|---|---|---|
| 97 host references ([`references/human-all97.tsv`](references/human-all97.tsv)) | T2T-CHM13v2.0 ([GCA_009914755.4](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/)); GRCh38.p14 ([GCF_000001405.40](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/)); [GENCODE 47](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/) transcripts; 47 HPRC year-1 samples, both haplotypes, from the [year-1 assembly index](https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index) (`f1_assembly_v2_genbank`). All files decompressed. | `REF_DIR` |
| 10 holdout haplotypes ([`references/holdout10.tsv`](references/holdout10.tsv)) | HPRC release 2 (`assemblies_release2_v1.0.index.csv`), gzipped FASTA | `HOLDOUT_REF_DIR` |
| WoL2 genomes ([`reads/wol2-genomes.list`](reads/wol2-genomes.list)) | [Web of Life 2](https://biocore.github.io/wol/), one `<accession>.fasta.gz` per genome | `WOL2_GENOME_DIR` |
| WoL2 short reads | one simulated Illumina pair per genome, `<accession>_R{1,2}.fastq.gz` | `WOL2_SR_READS_DIR` |
| fastp adapters | `ref/known_adapters.fna` of [cguccione/human_host_filtration](https://github.com/cguccione/human_host_filtration) @ `9807633` | `ADAPTERS` |

### 2. Indices

```bash
mkdir -p logs
for w in 10 20 25 35 45 50 100 150 200; do sbatch references/build-rype-index.sbatch $w; done
sbatch references/build-minimap2-indices.sbatch          # T2T_sr, T2T_hifi, GRCh38_hifi
sbatch references/build-cleanifier-index.sbatch exact
sbatch references/build-cleanifier-index.sbatch filter
```

The scripts above build these indices:
- **rype:** one bucket over all 97 references, k = 64, one index per window.
- **cleanifier:** over the same 97 references, built the way the cleanifier
  paper recommends (see the script header).
- **minimap2 (3 indices):** the two T2T-CHM13v2.0 indices for the rype → mm2
  rescue pass (`sr` and `map-hifi` presets). The long-read SOTA uses the
  T2T `map-hifi` index plus a GRCh38 `map-hifi` index.

Two more index types come from elsewhere:
- **SOTA short-read indices.** These are Qiita's human-pangenome minimap2
  databases, built by `scripts/download_references.sh` and
  `scripts/create_minimap2_indexes.sh` of
  [cguccione/human_host_filtration](https://github.com/cguccione/human_host_filtration):
  - `human-GRC-db.mmi`: GRCh38, `sr` preset;
  - `human-GCA-phix-db.mmi`: T2T-CHM13v2.0 plus phiX174, `sr` preset.

  Set `QIITA_HG38_MMI` and `QIITA_T2T_MMI` to them.
- **Movi 2 index.** This is a prebuilt `hprc-yr1` index: Movi 2.0.0,
  regular-thresholds mode, over CHM13 and both haplotypes of the 47 HPRC
  year-1 samples. Set `MOVI_INDEX` to it.

### 3. Reads

```bash
sbatch --array=1-10  reads/simulate-sr-holdout.sbatch    # ART, 1x PE150
sbatch --array=1-10  reads/qc-sr.sbatch holdout          # fastp
sbatch --array=1-160 reads/qc-sr.sbatch wol2
sbatch --array=1-10  reads/simulate-lr-holdout.sbatch    # pbsim3 + ccs, 1x HiFi
sbatch --array=1-160 reads/simulate-lr-wol2.sbatch
sbatch --array=1-10  reads/qc-lr.sbatch holdout          # fastplong
sbatch --array=1-160 reads/qc-lr.sbatch wol2
sbatch reads/build-concat-inputs.sbatch
```

Each QC step waits for its simulation, and the last step waits for all QC.
- **Concatenation.** `build-concat-inputs.sbatch` writes one input file (or
  pair) per read type and dataset. Every read ID is prefixed with its source
  tag, `<haplotype>::` or `chunk_NN::`, so each tool filters a whole dataset
  in one call and FN/FP are counted per tag.
- **WoL2 chunk lists.** The genome-to-chunk assignment used in the manuscript
  is in `reads/wol2-chunks/`. FP is pooled over chunks.
- **Simulation is not seeded.** ART seeds from the clock, and pbsim3 was run
  without a seed. Regenerated reads therefore differ read for read from the
  manuscript set, though not in distribution.

### 4. Benchmark runs

```bash
./run-main.sh               # every tool at the operating points        (24 jobs)
./run-window-sweep.sh       # rype and rype-mm2 at the other windows    (32 jobs)
./run-threshold-sweep.sh    # rype at w = 20, 11 thresholds             (4 jobs)
./run-laptop.sh 16          # 8 CPU / 16 GB                             (14 jobs)
./run-laptop.sh 8           # 8 CPU / 8 GB, long read                   (4 jobs)
./run-mem-scaling.sh        # long-read rype, 64 GB to 1 GB, 3 reps     (42 jobs)
```

- **Order.** Run `run-main.sh` first: the laptop and memory-scaling tables
  compare against its 64 GB runs.
- **Eval.** Each submitter queues `eval/eval.sbatch` after its jobs. The eval
  rebuilds every table in `data/` and redraws the figures.
- **Results.** Runs write to `${RESULTS_DIR}/<modality>/<dataset>/<tool>/`.

What each run measures:
- **Wall time** covers the whole call: gzipped FASTQ in, compressed
  survivors out.
- **Peak memory** is the summed RSS of the whole process tree
  (`eval/tree_rss.py`), so a pipeline's concurrent processes are all counted.
  For the multi-stage Movi 2 and SOTA arms, wall time is the sum of the
  stages and memory is the largest stage's.
- **FN and FP** are counted from each run's survivors: per haplotype for FN,
  pooled over the WoL2 chunks for FP. Short reads count pairs.

Compute at the operating points, from `data/*/op_runtime.tsv`, for one run on
8 CPUs:
- rype and rype → mm2: 5 to 25 minutes per dataset;
- cleanifier: 6 to 12 minutes;
- Movi 2 alone and SOTA: 0.4 to 3.7 hours.

### 5. Tables and figures

`eval/eval.sbatch` runs these in order:
1. `eval/eval_external.py`: read totals, and the Movi 2 and SOTA tables.
2. `eval/eval_op.py`: adds rype, rype → mm2 and cleanifier; writes `data/`.
3. `eval/eval_laptop.py` and `eval/eval_mem_scaling.py`.
4. `figures/make-figures.sh`.

A stage whose runs are not all complete prints a warning and leaves those
rows out. At the operating threshold, the threshold sweep must match the
operating-point run exactly; the eval warns if it does not.

## License

BSD 3-Clause. See [LICENSE](LICENSE).

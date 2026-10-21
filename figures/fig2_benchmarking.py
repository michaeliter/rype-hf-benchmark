#!/usr/bin/env python3
"""Main benchmarking figure: rype vs external host-filtration tools, short read
+ long read side by side, at each modality's operating point.

Four panels:
  (a) wall-clock time per 100k reads/pairs, holdout (human) + wol2 (microbial)
  (b) peak RSS, human holdout
  (c) on-disk reference index footprint
  (d) FN/FP error rates

rype numbers come from the streaming rype-hf pipeline (`rype classify
filter`, piped into minimap2 for the short-read rescue); one run per
configuration (mean if more replicates exist), memory is the peak summed RSS
of the whole pipeline's process tree. Reads the operating-point tables
(data/<mod>/op_*.tsv, eval/eval_op.py). Which tools appear, and in what
order, comes from figures/series.py.

Short read: rype bar = rype + mm2 rescue (w=20, threshold 0.0001).
Long read:  rype bar = rype alone (w=20, threshold 0.05, the FN=0 knee).
All tools: 8 CPU / 64 GB Slurm allocation, same concatenated inputs.
"""
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.patches as mpatches
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_style import (set_nature_style, create_figure, save_nature_figure,  # noqa: E402
                          add_panel_labels, lighten_color)
from figdata import DATA_ROOT, FIGDIR  # noqa: E402
import series as S  # noqa: E402

DATA = {"Short read": DATA_ROOT / "short-read", "Long read": DATA_ROOT / "long-read"}
GROUPS = ["Short read", "Long read"]
RYPE_APPROACH = {"Short read": "rype-mm2", "Long read": "rype"}
OP_WINDOW = 20
PER_READS = 100_000
GB, GIB_PER_KB = 1e9, 1.0 / 1048576.0


def rows(path):
    with open(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def approach_of(k, g):
    return RYPE_APPROACH[g] if k == "rype" else k


def color_of(k, g):
    return S.COLOR[approach_of(k, g)]


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
runtime = defaultdict(list)     # (group, dataset, approach) -> [s per 100k reads/pairs]
memory = defaultdict(list)      # (group, approach) -> [GiB], holdout only
for g in GROUPS:
    for r in rows(DATA[g] / "op_runtime.tsv"):
        runtime[(g, r["dataset"], r["approach"])].append(
            float(r["wall_seconds"]) / int(r["reads"]) * PER_READS)
        if r["dataset"] == "holdout" and r["peak_rss_kb"] not in ("", "NA"):
            memory[(g, r["approach"])].append(int(r["peak_rss_kb"]) * GIB_PER_KB)

fnr, fpr = {}, {}
for g in GROUPS:
    per_base = defaultdict(list)
    for r in rows(DATA[g] / "op_fn.tsv"):
        per_base[r["approach"]].append(float(r["fnr"]))
    for a, v in per_base.items():
        fnr[(g, a)] = statistics.mean(v)          # mean over holdout haplotypes
    for r in rows(DATA[g] / "op_fp.tsv"):
        fpr[(g, r["approach"])] = float(r["fpr"])

storage = {g: {} for g in GROUPS}
for g in GROUPS:
    for r in rows(DATA[g] / "storage.tsv"):
        key = f"rype_ryxdi_w{r['window']}" if r["component"] == "rype_ryxdi" else r["component"]
        storage[g][key] = int(r["bytes"])


def index_components(k, g):
    st = storage[g]
    a = approach_of(k, g)
    if k == "rype":
        comps = [("rype index", st[f"rype_ryxdi_w{OP_WINDOW}"])]
        if a == "rype-mm2":   # T2T-only .mmi (SR: not the SOTA arms' T2T + phiX file)
            comps.append(("T2T (minimap2)", st.get("rype_mm2_mmi", st["sota_t2t_mmi"])))
        return comps
    if a == "mm2-mm2-movi":
        hg38 = st.get("sota_hg38_mmi", st.get("sota_grch38_mmi"))
        return [("hg38 (minimap2)", hg38), ("T2T (minimap2)", st["sota_t2t_mmi"]),
                ("hprc-yr1 (Movi 2)", st["movi2_index"])]
    return [(label, st[comp]) for label, comp in S.INDEX_COMPONENTS[a]]


def complete(k):
    for g in GROUPS:
        a = approach_of(k, g)
        if not (runtime[(g, "holdout", a)] and runtime[(g, "wol2", a)]
                and (g, a) in fnr and (g, a) in fpr and memory[(g, a)]):
            return False
    return True


SERIES = ["rype"] + S.EXTERNAL
missing = [k for k in SERIES if not complete(k)]
for k in missing:
    print(f"WARNING: {k} has incomplete operating-point data -- left out of the figure", file=sys.stderr)
SERIES = [k for k in SERIES if k not in missing]
assert "rype" in SERIES, "rype operating-point runs missing; run arms/ + eval/eval_op.py first"

N = len(SERIES)
BAR_W = 0.8 / N
OFF = {k: (i - (N - 1) / 2) * BAR_W for i, k in enumerate(SERIES)}


def mean_sd(v):
    return statistics.mean(v), (statistics.pstdev(v) if len(v) > 1 else 0.0)


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
set_nature_style()
fig, axes = create_figure(width="double", nrows=2, ncols=2, aspect_ratio=0.9)
ax_rt, ax_mem, ax_store, ax_acc = axes.flatten()
# one run per configuration, so no error bars
bar_kw = dict(width=BAR_W, edgecolor="white", linewidth=0.5)
xs = np.arange(len(GROUPS))

# (a) runtime: modality x dataset
RT_GROUPS = [(g, ds) for g in GROUPS for ds in ("holdout", "wol2")]
RT_LABEL = {"holdout": "human", "wol2": "microbial"}
xs_rt = np.arange(len(RT_GROUPS))
rt_stats = {(g, ds, k): mean_sd(runtime[(g, ds, approach_of(k, g))]) for g, ds in RT_GROUPS for k in SERIES}
for k in SERIES:
    ms = [rt_stats[(g, ds, k)][0] for g, ds in RT_GROUPS]
    ax_rt.bar(xs_rt + OFF[k], ms, color=[color_of(k, g) for g, _ in RT_GROUPS], **bar_kw)
ax_rt.set_xticks(xs_rt)
ax_rt.set_xticklabels([f"{g}\n{RT_LABEL[ds]}" for g, ds in RT_GROUPS], fontsize=5.5)
ax_rt.set_ylabel(f"Wall-clock time per {PER_READS:,} reads/pairs (s)")
ax_rt.set_title("Runtime")
ax_rt.set_yscale("log")
all_rt = [v[0] for v in rt_stats.values()]
ax_rt.set_ylim(bottom=min(all_rt) / 3.0, top=max(all_rt) * 6.0)


def fmt_seconds(v):
    return f"{v:.0f}" if v >= 100 else f"{v:.1f}" if v >= 1 else f"{v:.2f}"


# value above every bar, rype's in bold
for k in SERIES:
    for x, (g, ds) in zip(xs_rt + OFF[k], RT_GROUPS):
        m = rt_stats[(g, ds, k)][0]
        ax_rt.annotate(fmt_seconds(m), (x, m), textcoords="offset points", xytext=(0, 3),
                       ha="center", va="bottom", fontsize=4.5, rotation=90,
                       fontweight="bold" if k == "rype" else "normal")
legend_handles = [mpatches.Patch(color=S.COLOR[RYPE_APPROACH[g]], label=f"{S.LABEL[RYPE_APPROACH[g]]} ({g.split()[0].lower()})")
                  for g in GROUPS]
legend_handles += [mpatches.Patch(color=S.COLOR[k], label=S.LABEL[k]) for k in SERIES if k != "rype"]
ax_rt.legend(handles=legend_handles, loc="upper left", fontsize=5)

# (b) peak RSS, holdout
mem_stats = {(g, k): mean_sd(memory[(g, approach_of(k, g))]) for g in GROUPS for k in SERIES}
for k in SERIES:
    ms = [mem_stats[(g, k)][0] for g in GROUPS]
    ax_mem.bar(xs + OFF[k], ms, color=[color_of(k, g) for g in GROUPS], **bar_kw)
    for x, m in zip(xs + OFF[k], ms):
        ax_mem.annotate(f"{m:.0f}", (x, m), textcoords="offset points", xytext=(0, 3),
                        ha="center", fontsize=4.5, rotation=90,
                        fontweight="bold" if k == "rype" else "normal")
ax_mem.set_xticks(xs); ax_mem.set_xticklabels(GROUPS)
ax_mem.set_ylabel("Peak resident memory (GiB)")
ax_mem.set_title("Peak RSS, human holdout")
ax_mem.set_ylim(top=max(m for m, _ in mem_stats.values()) * 1.35)

# (c) index footprint, stacked components
for k in SERIES:
    for gi, g in enumerate(GROUPS):
        x, bottom = gi + OFF[k], 0.0
        comps = index_components(k, g)
        total_gb = sum(b for _, b in comps) / GB
        for j, (clabel, b) in enumerate(comps):
            gb = b / GB
            shade = lighten_color(color_of(k, g), 0.55 * j / max(len(comps) - 1, 1))
            ax_store.bar(x, gb, bottom=bottom, width=BAR_W, color=shade, edgecolor="white", linewidth=0.5)
            # single-component bars are identified by the legend; only label stacks
            if len(comps) > 1 and gb / total_gb > 0.06:
                ax_store.annotate(f"{clabel}\n{gb:.2g} GB", (x, bottom + gb / 2), ha="center",
                                  va="center", fontsize=3.6, color="black")
            bottom += gb
        ax_store.annotate(f"{bottom:.1f}", (x, bottom), textcoords="offset points", xytext=(0, 3),
                          ha="center", fontsize=4.5, rotation=90,
                          fontweight="bold" if k == "rype" else "normal")
all_comp_gb = [b / GB for k in SERIES for g in GROUPS for _, b in index_components(k, g)]
all_tot_gb = [sum(b for _, b in index_components(k, g)) / GB for k in SERIES for g in GROUPS]
ax_store.set_xticks(xs); ax_store.set_xticklabels(GROUPS)
ax_store.set_ylabel("Reference index footprint (GB)")
ax_store.set_title("On-disk index size")
ax_store.set_yscale("log")
ax_store.set_ylim(bottom=min(all_comp_gb) / 3.0, top=max(all_tot_gb) * 3.0)

# (d) error rates on a log axis (rates near 0 span orders of magnitude, which
# recall/specificity near 1 would hide)
metrics = [("FN rate\n(1 − recall)", fnr), ("FP rate\n(1 − specificity)", fpr)]
dgroups = [(g, label, table) for g in GROUPS for label, table in metrics]
xs4 = np.arange(len(dgroups))
vals = {k: [table[(g, approach_of(k, g))] for g, _, table in dgroups] for k in SERIES}
floor_d = min(v for vs in vals.values() for v in vs if v > 0) / 3.0
for k in SERIES:
    plot_vals = [v if v > 0 else floor_d for v in vals[k]]
    ax_acc.bar(xs4 + OFF[k], plot_vals, width=BAR_W, color=[color_of(k, g) for g, _, _ in dgroups],
               edgecolor="white", linewidth=0.5)
    for x, v, pv in zip(xs4 + OFF[k], vals[k], plot_vals):
        ax_acc.annotate("0" if v == 0 else f"{v:.2e}", (x, pv), textcoords="offset points",
                        xytext=(0, 3), ha="center", fontsize=4, rotation=90, va="bottom",
                        fontweight="bold" if k == "rype" else "normal")
ax_acc.set_xticks(xs4); ax_acc.set_xticklabels([f"{g}\n{label}" for g, label, _ in dgroups], fontsize=5.5)
ax_acc.set_ylabel("Error rate")
ax_acc.set_yscale("log")
ax_acc.set_ylim(bottom=floor_d / 1.5, top=max(v for vs in vals.values() for v in vs) * 8)
ax_acc.set_title(f"FN/FP rate, w={OP_WINDOW}")

add_panel_labels(axes)
save_nature_figure(fig, str(FIGDIR / "fig_summary_sr_lr_rype_vs_sota"))

for g, ds in RT_GROUPS:
    print(f"runtime {g} {RT_LABEL[ds]} (s/100k): " +
          " ".join(f"{k}={rt_stats[(g, ds, k)][0]:.1f}" for k in SERIES))
for g in GROUPS:
    print(f"=== {g.lower()} ===")
    print("peak RSS (GiB): " + " ".join(f"{k}={mem_stats[(g, k)][0]:.1f}" for k in SERIES))
    print("storage (GB): " + " ".join(f"{k}={sum(b for _, b in index_components(k, g)) / GB:.1f}" for k in SERIES))
    print("FN rate: " + " ".join(f"{k}={fnr[(g, approach_of(k, g))]:.3e}" for k in SERIES))
    print("FP rate: " + " ".join(f"{k}={fpr[(g, approach_of(k, g))]:.3e}" for k in SERIES))

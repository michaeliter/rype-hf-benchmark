#!/usr/bin/env python3
"""Supplementary laptop-tier figure: the operating points at 8 CPU and 16 or 8 GB.

Two panels, styled like the main benchmarking figure (fig2_benchmarking.py):
  (a) wall-clock time per 100k reads/pairs, holdout (human) + wol2 (microbial),
      for the tools that ran in this tier. A black tick on each bar marks the
      same tool's 8 CPU / 64 GB run.
  (b) peak RSS on the human holdout (summed over the process tree), with the
      tier's job limit as a dashed line. Tools that do not fit are drawn as
      unfilled, hatched bars at their 64 GB peak.

Series and colors come from figures/series.py, as in the main figure; the
short-read rype bar is rype + mm2, the long-read one rype alone. FN/FP and
index size are not repeated: survivors are identical to the 64 GB runs
(checked through counts_match in the table).

Reads data/<mod>/laptop-<N>gb_runtime.tsv (eval/eval_laptop.py); writes
figures/output/fig_supp_laptop_<N>gb.{pdf,png}.

Usage: supp_laptop.py [--mem 16|8]      (default 16)
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_style import set_nature_style, create_figure, save_nature_figure, add_panel_labels  # noqa: E402
import series as S  # noqa: E402
from figdata import DATA_ROOT, FIGDIR  # noqa: E402

MODS = {"Short read": "short-read", "Long read": "long-read"}
TIER_GROUPS = {16: ["Short read", "Long read"], 8: ["Long read"]}
RYPE_APPROACH = {"Short read": "rype-mm2", "Long read": "rype"}
RT_LABEL = {"holdout": "human", "wol2": "microbial"}
PER_READS = 100_000
GIB_PER_KB = 1.0 / 1048576.0

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("--mem", type=int, choices=sorted(TIER_GROUPS), default=16)
MEM = ap.parse_args().mem
GROUPS = TIER_GROUPS[MEM]


def approach_of(k, g):
    return RYPE_APPROACH[g] if k == "rype" else k


def num(v):
    return None if v in ("", "NA", None) else float(v)


table = {}   # (group, dataset, approach) -> row
for g in GROUPS:
    path = DATA_ROOT / MODS[g] / f"laptop-{MEM}gb_runtime.tsv"
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            table[(g, r["dataset"], r["approach"])] = r

SERIES = ["rype"] + S.EXTERNAL
for g in GROUPS:
    for ds in RT_LABEL:
        for k in SERIES:
            assert (g, ds, approach_of(k, g)) in table, f"{g} {ds} {k}: no row; run eval/eval_laptop.py --mem {MEM}"
N = len(SERIES)
BAR_W = 0.8 / N
OFF = {k: (i - (N - 1) / 2) * BAR_W for i, k in enumerate(SERIES)}


def fmt_seconds(v):
    return f"{v:.0f}" if v >= 100 else f"{v:.1f}" if v >= 1 else f"{v:.2f}"


def weight(k):
    return "bold" if k == "rype" else "normal"


def per_reads(r, col):
    s, n = num(r[col]), num(r["reads"])
    return None if s is None or n is None else s / n * PER_READS


set_nature_style()
fig, (ax_rt, ax_mem) = create_figure(width="double", nrows=1, ncols=2, aspect_ratio=0.5)
bar_kw = dict(width=BAR_W, edgecolor="white", linewidth=0.5)
note_kw = dict(textcoords="offset points", xytext=(0, 2), ha="center", va="bottom",
               fontsize=5.5, rotation=90, color="0.35")
NOTE = {"missing": "not run"}

# (a) runtime: modality x dataset, laptop bars with 64 GB ticks. Tools that do
# not fit in this tier are left out here; panel (b) shows why.
RT_GROUPS = [(g, ds) for g in GROUPS for ds in RT_LABEL]
RT_SERIES = [k for k in SERIES if any(table[(g, ds, approach_of(k, g))]["status"] != "does_not_fit"
                                      for g, ds in RT_GROUPS)]
RT_OFF = {k: (i - (len(RT_SERIES) - 1) / 2) * BAR_W for i, k in enumerate(RT_SERIES)}
xs_rt = np.arange(len(RT_GROUPS))
rt_vals = []
for k in RT_SERIES:
    for x, (g, ds) in zip(xs_rt + RT_OFF[k], RT_GROUPS):
        r = table[(g, ds, approach_of(k, g))]
        color = S.COLOR[approach_of(k, g)]
        if r["status"] == "ok":
            v = per_reads(r, "wall_seconds")
            ax_rt.bar(x, v, color=color, **bar_kw)
            rt_vals.append(v)
            srv = per_reads(r, "server_wall_seconds")
            if srv is not None:
                ax_rt.plot([x - BAR_W * 0.4, x + BAR_W * 0.4], [srv, srv], color="black", linewidth=0.8)
                rt_vals.append(srv)
            # value above the bar or its 64 GB tick, whichever is higher; rype's in bold
            ax_rt.annotate(fmt_seconds(v), (x, max(v, srv or 0)), textcoords="offset points",
                           xytext=(0, 3), ha="center", va="bottom", fontsize=6, rotation=90,
                           fontweight=weight(k))
for k in RT_SERIES:   # notes after the bars, once the y-range is known
    for x, (g, ds) in zip(xs_rt + RT_OFF[k], RT_GROUPS):
        r = table[(g, ds, approach_of(k, g))]
        if r["status"] != "ok":
            ax_rt.annotate(NOTE.get(r["status"], r["status"]), (x, min(rt_vals) / 2.5), **note_kw)
ax_rt.set_xticks(xs_rt)
ax_rt.set_xticklabels([f"{g}\n{RT_LABEL[ds]}" for g, ds in RT_GROUPS], fontsize=7)
ax_rt.set_ylabel(f"Wall-clock time per\n{PER_READS:,} reads/pairs (s)")
ax_rt.set_title(f"Runtime, 8 CPU / {MEM} GB")
ax_rt.set_yscale("log")
ax_rt.set_ylim(bottom=min(rt_vals) / 3.0, top=max(rt_vals) * 6.0)
handles = [mpatches.Patch(color=S.COLOR[RYPE_APPROACH[g]], label=f"{S.LABEL[RYPE_APPROACH[g]]} ({g.split()[0].lower()})")
           for g in GROUPS]
handles += [mpatches.Patch(color=S.COLOR[k], label=S.LABEL[k]) for k in SERIES if k != "rype"]
handles.append(mpatches.Patch(fill=False, hatch="////", edgecolor="0.4", linewidth=0.6,
                              label=f"does not fit in {MEM} GB (64 GB peak)"))
handles.append(mlines.Line2D([], [], color="black", linewidth=0.8, label="same tool, 8 CPU / 64 GB"))
ax_rt.legend(handles=handles, loc="upper left", fontsize=6.5)

# (b) peak RSS, holdout, against the job limit. Labels sit on a white patch
# so the dashed limit line passes behind them, not through them.
label_bg = dict(boxstyle="square,pad=0.05", facecolor="white", edgecolor="none")
xs = np.arange(len(GROUPS))
mem_vals = [float(MEM)]
for k in SERIES:
    for x, g in zip(xs + OFF[k], GROUPS):
        r = table[(g, "holdout", approach_of(k, g))]
        color = S.COLOR[approach_of(k, g)]
        if r["status"] == "ok":
            v = num(r["peak_rss_kb"]) * GIB_PER_KB
            ax_mem.bar(x, v, color=color, **bar_kw)
            ax_mem.annotate(f"{v:.1f}", (x, v), textcoords="offset points", xytext=(0, 3),
                            ha="center", fontsize=6, rotation=90, fontweight=weight(k),
                            bbox=label_bg, zorder=4)
            mem_vals.append(v)
        elif r["status"] == "does_not_fit" and num(r["server_peak_rss_kb"]) is not None:
            v = num(r["server_peak_rss_kb"]) * GIB_PER_KB
            ax_mem.bar(x, v, width=BAR_W, fill=False, hatch="////", edgecolor=color, linewidth=0.6)
            ax_mem.annotate(f"{v:.0f}", (x, v), textcoords="offset points", xytext=(0, 3),
                            ha="center", fontsize=6, rotation=90, color="0.35", bbox=label_bg, zorder=4)
            mem_vals.append(v)
        else:
            ax_mem.annotate(NOTE.get(r["status"], r["status"]), (x, 0), **note_kw)
ax_mem.axhline(MEM, color="0.3", linestyle="--", linewidth=0.6, zorder=3)
ax_mem.annotate(f"{MEM} GiB\njob limit", (1.0, MEM), xycoords=("axes fraction", "data"),
                textcoords="offset points", xytext=(3, 0), ha="left", va="center", fontsize=6.5, color="0.3",
                annotation_clip=False)
ax_mem.set_xticks(xs); ax_mem.set_xticklabels(GROUPS)
ax_mem.set_ylabel("Peak resident memory (GiB)")
ax_mem.set_title("Peak RSS, human holdout")
ax_mem.set_ylim(top=max(mem_vals) * 1.3)

add_panel_labels(np.array([ax_rt, ax_mem]))
save_nature_figure(fig, str(FIGDIR / f"fig_supp_laptop_{MEM}gb"))

for g, ds in RT_GROUPS:
    print(f"{MEM} GB runtime {g} {RT_LABEL[ds]} (s/100k): " + " ".join(
        f"{k}={per_reads(table[(g, ds, approach_of(k, g))], 'wall_seconds') or table[(g, ds, approach_of(k, g))]['status']}"
        for k in SERIES))

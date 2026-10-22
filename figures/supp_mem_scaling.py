#!/usr/bin/env python3
"""Supplementary memory-scaling figure: long-read rype at 8 CPU and job memory
limits from 64 down to 1 GB, three replicates per limit.

Three panels, x = job memory limit (log2):
  (a) wall-clock time per 100k reads, human holdout and microbial wol2
      (line: mean of the replicates that finished; dots: each replicate)
  (b) peak RSS of the process tree, with the job limit (dashed); rype-hf
      gives rype 90% of the limit minus 512 MiB (stated in the methods)
  (c) rype classification passes: a smaller budget means more passes over
      the index, which is what costs time
Runs that did not finish (OOM-killed or failed) are marked with an x at the
bottom of panel (a), and at their recorded peak in panel (b) when there is one.
Survivors are identical to the 64 GB operating-point run at every limit that
finished (counts_match in the table), so accuracy is not repeated.

Reads data/long-read/mem_scaling.tsv (eval/eval_mem_scaling.py); writes
figures/output/fig_supp_mem_scaling.{pdf,png}.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.lines as mlines
from matplotlib.ticker import NullLocator
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_style import set_nature_style, create_figure, save_nature_figure, add_panel_labels  # noqa: E402
import series as S  # noqa: E402
from figdata import DATA_ROOT, FIGDIR  # noqa: E402

PER_READS = 100_000
GIB_PER_KB = 1.0 / 1048576.0
DATASETS = {"holdout": dict(label="human (holdout)", short="human", marker="o", ls="-"),
            "wol2": dict(label="microbial (wol2)", short="microbial", marker="s", ls="--")}
COLOR = S.COLOR["rype"]

FAIL_COLOR = "0.35"

rows = list(csv.DictReader(open(DATA_ROOT / "long-read" / "mem_scaling.tsv"), delimiter="\t"))
tiers = sorted({int(r["tier_gb"]) for r in rows})
by = defaultdict(list)          # (dataset, tier) -> rows
for r in rows:
    by[(r["dataset"], int(r["tier_gb"]))].append(r)


def ok(rs):
    return [r for r in rs if r["status"] == "ok"]


def failed(rs):
    """Runs that did not finish: OOM-killed, failed, or killed with no record.
    'incomplete' (rype-hf cut short after the run) and 'missing' are left out."""
    return [r for r in rs if r["status"] not in ("ok", "missing", "incomplete")]


def per_reads(r):
    return float(r["wall_seconds"]) / int(r["reads"]) * PER_READS


def num(v):
    return None if v in ("", "NA", None) else float(v)


set_nature_style()
fig, (ax_rt, ax_mem, ax_pass) = create_figure(width="double", nrows=1, ncols=3, aspect_ratio=0.42)
line_kw = dict(color=COLOR, lw=0.9, ms=3)
dot_kw = dict(color=COLOR, ms=2.2, alpha=0.55, lw=0, zorder=4)

# (a) runtime
rt_all = []
for ds, st in DATASETS.items():
    xs, ys = [], []
    for t in tiers:
        vals = [per_reads(r) for r in ok(by[(ds, t)])]
        if vals:
            xs.append(t); ys.append(np.mean(vals)); rt_all += vals
            ax_rt.plot([t] * len(vals), vals, marker=st["marker"], mfc="none", **dot_kw)
    ax_rt.plot(xs, ys, marker=st["marker"], ls=st["ls"], **line_kw)
lo, hi = min(rt_all), max(rt_all)
log_rt = hi / lo > 4
if log_rt:
    ax_rt.set_yscale("log")
    ax_rt.set_ylim(lo / 2.5, hi * 2)
    fail_y = {"holdout": lo / 1.9, "wol2": lo / 1.9 * 1.45}
else:
    ax_rt.set_ylim(0, hi * 1.25)
    fail_y = {"holdout": hi * 0.04, "wol2": hi * 0.10}
for ds, st in DATASETS.items():
    for t in tiers:
        n = len(failed(by[(ds, t)]))
        if n:
            ax_rt.plot(t, fail_y[ds], marker="x", color=FAIL_COLOR, ms=3.5, mew=0.8, lw=0)
            ax_rt.annotate(f"{st['short']} {n}/{len(by[(ds, t)])}", (t, fail_y[ds]),
                           textcoords="offset points", xytext=(4, 0), ha="left", va="center",
                           fontsize=6, color=FAIL_COLOR)
ax_rt.set_ylabel(f"Time (s) per\n{PER_READS // 1000}k reads")
ax_rt.set_title("Runtime")
handles = [mlines.Line2D([], [], marker=st["marker"], ls=st["ls"], label=st["label"], **line_kw)
           for st in DATASETS.values()]
any_failed = any(failed(rs) for rs in by.values())
if any_failed:
    handles.append(mlines.Line2D([], [], marker="x", color=FAIL_COLOR, lw=0, ms=3.5, mew=0.8,
                                 label="did not finish (n/3)"))
ax_rt.legend(handles=handles, loc="upper right", fontsize=6.5)

# (b) peak RSS against the job limit
mem_all = []
for ds, st in DATASETS.items():
    xs, ys = [], []
    for t in tiers:
        vals = [int(r["peak_rss_kb"]) * GIB_PER_KB for r in ok(by[(ds, t)])]
        if vals:
            xs.append(t); ys.append(np.mean(vals)); mem_all += vals
            ax_mem.plot([t] * len(vals), vals, marker=st["marker"], mfc="none", **dot_kw)
        for r in failed(by[(ds, t)]):
            if num(r["peak_rss_kb"]):
                v = num(r["peak_rss_kb"]) * GIB_PER_KB
                ax_mem.plot(t, v, marker="x", color=FAIL_COLOR, ms=3.5, mew=0.8, lw=0, zorder=5)
                mem_all.append(v)
    ax_mem.plot(xs, ys, marker=st["marker"], ls=st["ls"], **line_kw)
tx = np.array(tiers, float)
ax_mem.plot(tx, tx, color="0.3", ls=(0, (4, 2)), lw=0.8, zorder=1)
mem_handles = [mlines.Line2D([], [], color="0.3", ls=(0, (4, 2)), lw=0.8, label="job memory limit")]
if any(num(r["peak_rss_kb"]) for rs in by.values() for r in failed(rs)):
    mem_handles.append(mlines.Line2D([], [], marker="x", color=FAIL_COLOR, lw=0, ms=3.5, mew=0.8,
                                     label="peak before OOM kill"))
ax_mem.legend(handles=mem_handles, loc="upper left", fontsize=6.5)
ax_mem.set_yscale("log", base=2)
mem_floor = min(mem_all)
ax_mem.set_ylim(mem_floor / 1.6, tx[-1] * 1.6)
ax_mem.set_ylabel("Peak RSS (GiB)")
ax_mem.set_title("Peak RSS")

# (c) classification passes
for ds, st in DATASETS.items():
    xs, ys = [], []
    for t in tiers:
        vals = [int(r["passes"]) for r in ok(by[(ds, t)]) if r["passes"] != "NA"]
        if vals:
            xs.append(t); ys.append(np.mean(vals))
            ax_pass.plot([t] * len(vals), vals, marker=st["marker"], mfc="none", **dot_kw)
    ax_pass.plot(xs, ys, marker=st["marker"], ls=st["ls"], **line_kw)
ax_pass.set_yscale("log")
ax_pass.set_ylabel("Passes")
ax_pass.set_title("Passes over the index")

for ax in (ax_rt, ax_mem, ax_pass):
    ax.set_xscale("log", base=2)
    ax.set_xticks(tiers)
    ax.set_xticklabels([str(t) for t in tiers])
    ax.set_xlim(tiers[0] / 1.5, tiers[-1] * 1.5)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Job memory limit (GB)")
mem_ticks = [t for t in (0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
             if ax_mem.get_ylim()[0] <= t <= ax_mem.get_ylim()[1]]
ax_mem.set_yticks(mem_ticks)
ax_mem.set_yticklabels([f"{t:g}" for t in mem_ticks])
ax_mem.yaxis.set_minor_locator(NullLocator())

fig.suptitle("Long-read RYpe, 8 CPU", fontsize=9)
add_panel_labels([ax_rt, ax_mem, ax_pass], x=-0.3, y=1.22)
save_nature_figure(fig, str(FIGDIR / "fig_supp_mem_scaling"))
for ds in DATASETS:
    summ = []
    for t in tiers:
        rs = by[(ds, t)]
        v = [per_reads(r) for r in ok(rs)]
        summ.append(f"{t}GB={np.mean(v):.1f}" if v else f"{t}GB=({len(failed(rs))} failed)")
    print(f"runtime {ds} (s/100k): " + " ".join(summ))

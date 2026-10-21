#!/usr/bin/env python3
"""Supplementary figure: FN, FP and runtime vs rype window size, short read and
long read side by side.

3x2 panels:
  (a, b) FN rate vs window (10 holdout haplotypes, mean +- s.d.)
  (c, d) FP rate vs window (wol2, pooled)
  (e, f) wall time vs window, holdout and wol2 (one run each)
rype and rype + mm2 per window; window-independent tools (series.py) as flat
lines. Short-read windows {10, 20, 25, 35, 45}; long-read {20, 50, 100, 150, 200}.

Reads data/<mod>/{fn,fp,runtime}.tsv (eval/eval_op.py).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_style import (set_nature_style, create_figure, save_nature_figure,  # noqa: E402
                          add_panel_labels)
from figdata import FIGDIR, LR, SR  # noqa: E402

set_nature_style()


def present(mod):
    """Approaches from series.py that actually have rows in this modality's data
    (a newly registered tool shows up once eval_op.py has written its rows)."""
    have = {a for a, _ in mod.load_fn()}
    return [a for a in mod.APPROACHES if a in have]


def clip_lower(means, sds, floor):
    means = np.asarray(means, float); sds = np.asarray(sds, float)
    lower = means - np.maximum(means - sds, floor)
    lower = np.clip(lower, 0.0, None)
    return np.vstack([lower, sds])


def plot_fn(ax, mod, title, unit="pair"):
    fn = mod.load_fn()
    allm = [mod.mean_sd(fn[(a, w)]["fnr"])[0] for a in mod.APPROACHES for w in mod.WINDOWS
            if fn.get((a, w))]
    floor = min(m for m in allm if m > 0) / 3.0
    any_zero = False
    for a in present(mod):
        means = np.array([mod.mean_sd(fn.get((a, w), {"fnr": []})["fnr"])[0] for w in mod.WINDOWS])
        sds = np.array([mod.mean_sd(fn.get((a, w), {"fnr": []})["fnr"])[1] for w in mod.WINDOWS])
        if a in mod.WINDOW_INDEPENDENT:
            m, s = float(np.nanmean(means)), float(np.nanmean(sds))
            ax.axhline(m, color=mod.COLOR[a], ls=":", lw=1.0, label=mod.LABEL[a], zorder=2)
            ax.axhspan(max(m - s, floor), m + s, color=mod.COLOR[a], alpha=0.08, zorder=0)
        else:
            # exact-0 means can't be shown on a log axis; clip to the shared floor
            # and ring them with an open marker so they read as "measured zero",
            # not as missing data (matplotlib otherwise silently drops the point).
            is_zero = means == 0
            any_zero = any_zero or is_zero.any()
            means_plot = np.where(is_zero, floor, means)
            ax.errorbar(mod.WINDOWS, means_plot, yerr=clip_lower(means_plot, sds, floor),
                        color=mod.COLOR[a], marker=mod.MARKER[a], markersize=3,
                        linewidth=0.9, capsize=2, label=mod.LABEL[a], zorder=3)
            for w, y, z in zip(mod.WINDOWS, means_plot, is_zero):
                if z:
                    ax.plot(w, y, marker="o", ms=6, mfc="none", mec=mod.COLOR[a], mew=1.0, zorder=4)
    ax.set_yscale("log"); ax.set_ylim(bottom=floor)
    if any_zero:
        ax.annotate("open marker = measured FN = 0", (0.02, 0.02), xycoords="axes fraction",
                    fontsize=4.5, color="0.4", ha="left", va="bottom")
    ax.set_xlabel("rype window size, w"); ax.set_ylabel(f"Per-{unit} false-negative rate")
    ax.set_xticks(mod.WINDOWS); ax.set_title(title)


def plot_fp(ax, mod, title, unit="pair"):
    fp = mod.load_fp()
    allm = [fp[(a, w)]["fpr"] for a in mod.APPROACHES for w in mod.WINDOWS if (a, w) in fp]
    floor = min(m for m in allm if m > 0) / 3.0
    for a in present(mod):
        ys = np.array([fp.get((a, w), {"fpr": np.nan})["fpr"] for w in mod.WINDOWS])
        if a in mod.WINDOW_INDEPENDENT:
            m = float(np.nanmean(ys))
            ax.axhline(m, color=mod.COLOR[a], ls=":", lw=1.0, label=mod.LABEL[a], zorder=2)
        else:
            ax.plot(mod.WINDOWS, ys, color=mod.COLOR[a], marker=mod.MARKER[a], markersize=3,
                    linewidth=0.9, label=mod.LABEL[a], zorder=3)
    ax.set_yscale("log"); ax.set_ylim(bottom=floor)
    ax.set_xlabel("rype window size, w"); ax.set_ylabel(f"Per-{unit} false-positive rate")
    ax.set_xticks(mod.WINDOWS); ax.set_title(title)


def plot_runtime(ax, mod, title, unit="pair"):
    # One run per tool, window and dataset, so no error bars.
    rt = mod.load_runtime()      # {(dataset, approach, window): [sec per 1M reads/pairs, ...]}
    for a in present(mod):
        for ds, ls in (("holdout", "-"), ("wol2", "--")):
            ys = np.array([mod.mean_sd(rt.get((ds, a, w), []))[0] for w in mod.WINDOWS])
            if a in mod.WINDOW_INDEPENDENT:
                m = float(np.nanmean([y for y in ys if y == y]))
                if m == m:
                    ax.axhline(m, color=mod.COLOR[a], ls=ls, lw=1.0)
            else:
                ax.plot(mod.WINDOWS, ys, color=mod.COLOR[a], ls=ls,
                        marker=mod.MARKER[a], ms=3, lw=0.9)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("rype window size, w"); ax.set_ylabel(f"Wall time (s per 1M {unit}s)")
    ax.set_xticks(mod.WINDOWS); ax.set_title(title)
    from matplotlib.lines import Line2D
    h = [Line2D([0], [0], color=mod.COLOR[a], marker=mod.MARKER[a], lw=0.9, label=mod.LABEL[a])
         for a in present(mod)]
    h += [Line2D([0], [0], color="0.3", ls="-", label="holdout (FN set)"),
          Line2D([0], [0], color="0.3", ls="--", label="wol2 (FP set)")]
    ax.legend(handles=h, loc="best", fontsize=5)


fig, axes = create_figure(width="double", nrows=3, ncols=2, aspect_ratio=1.35)
(ax_fn_sr, ax_fn_lr), (ax_fp_sr, ax_fp_lr), (ax_rt_sr, ax_rt_lr) = axes

plot_fn(ax_fn_sr, SR, "Short read (ART), FN vs window", unit="pair")
plot_fn(ax_fn_lr, LR, "Long read (pbsim3), FN vs window", unit="read")
plot_fp(ax_fp_sr, SR, "Short read (ART), FP vs window", unit="pair")
plot_fp(ax_fp_lr, LR, "Long read (pbsim3), FP vs window", unit="read")
plot_runtime(ax_rt_sr, SR, "Short read (ART), runtime vs window", unit="pair")
plot_runtime(ax_rt_lr, LR, "Long read (pbsim3), runtime vs window", unit="read")

for ax in (ax_fn_sr, ax_fn_lr, ax_fp_sr, ax_fp_lr):
    ax.legend(loc="best", fontsize=5)

add_panel_labels(axes)
save_nature_figure(fig, str(FIGDIR / "fig_supp_window_sweep"))
print("wrote fig_supp_window_sweep")

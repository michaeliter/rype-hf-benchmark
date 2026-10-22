#!/usr/bin/env python3
"""Supplementary figure: rype decision-threshold sweep, short read + long read, both at w=20.

2-panel figure:
  (a) short read (ART holdout/wol2): rype FN/FP rate vs threshold, w=20, with
      Movi SOTA FN/FP as flat dashed references. Deployed classify
      threshold (1e-4) marked.
  (b) long read (pbsim3 holdout/wol2): rype FN/FP rate vs threshold, w=20, with
      Movi SOTA references. Deployed operating threshold (0.05, the FN=0 knee)
      marked.

Both panels share axes scaling for direct visual comparison of how far the
short-read curve sits from a corner (motivating the mm2 rescue step) versus how
the long-read curve reaches a good corner on its own.

Reads data/<mod>/threshold_sweep_w20.tsv and data/<mod>/movi_ref.tsv
(eval/eval_op.py).
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_style import (set_nature_style, create_figure, save_nature_figure,  # noqa: E402
                          add_panel_labels)
from figdata import FIGDIR, LR, SR  # noqa: E402

FN_COLOR, FP_COLOR = "#0077BB", "#CC3311"
W = 20


def load_sweep(mod):
    rows = list(csv.DictReader(open(mod.DATA / f"threshold_sweep_w{W}.tsv"), delimiter="\t"))
    thr = [float(r["threshold"]) for r in rows]
    fnr = [float(r["fnr"]) for r in rows]
    fpr = [float(r["fpr"]) for r in rows]
    return thr, fnr, fpr


def plot_panel(ax, mod, title, op_thr, op_label, floor, unit="pair"):
    thr, fnr, fpr = load_sweep(mod)
    movi = mod.movi_ref()
    # FN=0 can't be shown on a log axis; clip to the shared floor and mark it so
    # the point isn't silently dropped (matplotlib otherwise draws a spurious
    # vertical line down to the axis edge for non-positive log-scale values).
    fnr_plot = [max(v, floor) for v in fnr]
    is_zero = [v == 0 for v in fnr]

    ax.axvline(op_thr, color="0.35", ls="-", lw=0.7, alpha=0.8, zorder=1)
    ax.annotate(op_label, (op_thr, 1.0), xycoords=("data", "axes fraction"),
                textcoords="offset points", va="top",
                # left of the line when there is room, so it clears the FN curve
                xytext=(-3, -2) if op_thr > 1e-2 else (3, -2), ha="right" if op_thr > 1e-2 else "left",
                fontsize=7, color="0.35")
    ax.plot(thr, fnr_plot, marker="o", ms=3, lw=0.9, color=FN_COLOR, label="RYpe FN rate", zorder=3)
    for t, y, z in zip(thr, fnr_plot, is_zero):
        if z:
            ax.plot(t, y, marker="o", ms=5, mfc="none", mec=FN_COLOR, mew=1.0, zorder=4)
    if any(is_zero):
        # Label the run of zero-FN points from above, centred on it (log x), so
        # it stays clear of the x axis.
        zeros = [t for t, z in zip(thr, is_zero) if z]
        mid_t = (min(zeros) * max(zeros)) ** 0.5
        ax.annotate("RYpe FN = 0", (mid_t, floor), textcoords="offset points",
                    xytext=(0, 4), ha="center", va="bottom", fontsize=6.5, color=FN_COLOR)
    if movi["fnr"] > 0:
        ax.axhline(movi["fnr"], ls="--", lw=0.8, color=FN_COLOR, alpha=0.6, label="movi SOTA FN rate")
    else:
        ax.axhline(floor, ls="--", lw=0.8, color=FN_COLOR, alpha=0.6, label="movi SOTA FN rate = 0")
    ax.plot(thr, fpr, marker="s", ms=3, lw=0.9, color=FP_COLOR, label="RYpe FP rate", zorder=3)
    ax.axhline(movi["fpr"], ls="--", lw=0.8, color=FP_COLOR, alpha=0.6, label="movi SOTA FP rate")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("RYpe decision threshold")
    ax.set_ylabel(f"Per-{unit} error rate")
    ax.set_title(title)
    # Raised off the bottom so the movi SOTA FN line (and the zero-FN floor)
    # pass below it rather than through it.
    ax.legend(loc="lower left", bbox_to_anchor=(0, 0.13), fontsize=6.5)


set_nature_style()
fig, (ax_sr, ax_lr) = create_figure(width="double", ncols=2, aspect_ratio=0.62)

# shared floor for FN=0 points, well below the smallest genuinely-measured nonzero rate
# -- must also account for the SOTA reference rates (axhline'd in plot_panel), not just
# the sweep's own values, or a SOTA rate smaller than the sweep floor gets clipped off
# the bottom of the axis and silently disappears (this happened to short-read SOTA FN,
# 8.99e-7, well below rype's own smallest measured FN rate in the sweep, ~4.6e-5).
_, sr_fnr, sr_fpr = load_sweep(SR)
_, lr_fnr, lr_fpr = load_sweep(LR)
sr_movi, lr_movi = SR.movi_ref(), LR.movi_ref()
all_pos = [v for v in sr_fnr + sr_fpr + lr_fnr + lr_fpr if v > 0]
all_pos += [v for v in (sr_movi.get("fnr"), sr_movi.get("fpr"), lr_movi.get("fnr"), lr_movi.get("fpr")) if v]
FLOOR = min(all_pos) / 3.0

plot_panel(ax_sr, SR, f"Short read (ART), w={W}", SR.OP_THR, f"operating\nthreshold {SR.OP_THR:g}", FLOOR, unit="pair")
plot_panel(ax_lr, LR, f"Long read (pbsim3), w={W}", LR.OP_THR, f"operating\nthreshold {LR.OP_THR:g}", FLOOR, unit="read")

# shared y-limits so the panels are directly comparable
ax_sr.set_ylim(bottom=FLOOR / 1.5); ax_lr.set_ylim(bottom=FLOOR / 1.5)
ylo = min(ax_sr.get_ylim()[0], ax_lr.get_ylim()[0])
yhi = max(ax_sr.get_ylim()[1], ax_lr.get_ylim()[1])
ax_sr.set_ylim(ylo, yhi); ax_lr.set_ylim(ylo, yhi)

add_panel_labels([ax_sr, ax_lr])
save_nature_figure(fig, str(FIGDIR / "fig_supp_threshold_sweep"))
print("wrote fig_supp_threshold_sweep")

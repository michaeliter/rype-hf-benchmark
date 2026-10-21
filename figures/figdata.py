"""Paths and table loaders shared by the figure scripts.

Tables are read from data/<modality>/ (written by eval/eval_op.py); figures
are written to figures/output/. SR and LR expose the same interface: window
grid, operating threshold and loaders for the sweep tables. Approach order,
colors and labels come from series.py.
"""
import csv
import statistics
from collections import defaultdict
from pathlib import Path

from series import APPROACHES, COLOR, LABEL, MARKER, WINDOW_INDEPENDENT

REPO = Path(__file__).resolve().parent.parent
DATA_ROOT = REPO / "data"
FIGDIR = REPO / "figures" / "output"
FIGDIR.mkdir(parents=True, exist_ok=True)


def _rows(path):
    with open(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


class Modality:
    APPROACHES, COLOR, LABEL, MARKER, WINDOW_INDEPENDENT = APPROACHES, COLOR, LABEL, MARKER, WINDOW_INDEPENDENT

    def __init__(self, name, windows, op_thr):
        self.DATA = DATA_ROOT / name
        self.WINDOWS = windows      # window sweep x-axis; w=20 is the operating point
        self.OP_THR = op_thr        # operating decision threshold

    @staticmethod
    def mean_sd(xs):
        return (statistics.mean(xs) if xs else float("nan"),
                statistics.pstdev(xs) if len(xs) >= 2 else 0.0)

    def load_fn(self):      # -> {(approach, window): {"fnr": [per haplotype], "fn": [...]}}
        out = defaultdict(lambda: {"fnr": [], "fn": []})
        for r in _rows(self.DATA / "fn.tsv"):
            if r["approach"] in APPROACHES:
                k = (r["approach"], int(r["window"]))
                out[k]["fnr"].append(float(r["fnr"])); out[k]["fn"].append(int(r["fn"]))
        return out

    def load_fp(self):      # -> {(approach, window): {"fp": int, "fpr": float}}
        return {(r["approach"], int(r["window"])): {"fp": int(r["fp"]), "fpr": float(r["fpr"])}
                for r in _rows(self.DATA / "fp.tsv") if r["approach"] in APPROACHES}

    def load_runtime(self):     # -> {(dataset, approach, window): [s per million reads or pairs]}
        out = defaultdict(list)
        for r in _rows(self.DATA / "runtime.tsv"):
            if r["approach"] in APPROACHES and r.get("reads"):
                out[(r["dataset"], r["approach"], int(r["window"]))].append(
                    float(r["wall_seconds"]) / int(r["reads"]) * 1e6)
        return out

    def movi_ref(self):     # -> {"fnr": rate, "fpr": rate}, SOTA reference rates
        return {r["metric"]: float(r["rate"]) for r in _rows(self.DATA / "movi_ref.tsv")}


# Short read counts pairs; long read counts reads.
SR = Modality("short-read", [10, 20, 25, 35, 45], 0.0001)
LR = Modality("long-read", [20, 50, 100, 150, 200], 0.05)

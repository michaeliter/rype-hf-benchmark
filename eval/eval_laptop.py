#!/usr/bin/env python3
"""Laptop-tier tables (8 CPU, 16 or 8 GB): data/<mod>/laptop-<N>gb_runtime.tsv.

One row per dataset x approach:
  status            ok | failed (exit N) | missing | does_not_fit
  wall_seconds      tree_rss.tsv wall time of the timed call
  peak_rss_kb       peak summed RSS of the process tree (tree_rss.tsv)
  peak_proc_rss_kb  largest single process (/usr/bin/time -v)
  reads             reads (LR) or pairs (SR) in the dataset, from op_runtime.tsv
  server_wall_seconds, server_peak_rss_kb   the same approach in the 64 GB
                    op_runtime.tsv (run eval_op.py first)
  counts_match      rype/rype-mm2 only: stats.tsv identical to the 64 GB run
                    (the memory budget must not change which reads survive)

Tools whose 64 GB peak exceeds the tier's limit (Movi2 and mm2->mm2->Movi;
cleanifier exact at 8 GB) are listed as does_not_fit with their 64 GB peak. So is
rype-mm2 wherever minimap2 alone leaves no room for rype: its long-read T2T
map-hifi index peaks at 13.8 GB, the SR index at 11.4 GB (and the long-read
operating point is rype alone).

The 16 GB tier covers both modalities; the 8 GB tier is long read only.

Usage: source common.sh first (it exports the paths), then
       eval_laptop.py [--mem 16|8]      (default 16)
"""
import argparse
import csv
import os
import sys
from pathlib import Path

def env(name):
    try:
        return os.environ[name]
    except KeyError:
        sys.exit(f"{Path(__file__).name}: {name} is not set; source common.sh first")


REPO, RESULTS = Path(env("REPO")), Path(env("RESULTS_DIR"))
TIERS = {
    16: {"mods": ("short-read", "long-read"),
         "run": ("rype", "rype-mm2", "cleanifier", "cleanifier-filter"),
         "not_run": ("movi2-alone", "mm2-mm2-movi")},
    8: {"mods": ("long-read",), "run": ("rype", "cleanifier-filter"),
        "not_run": ("cleanifier", "movi2-alone", "mm2-mm2-movi")},
}
# minimap2's own footprint, not a 64 GB peak, is why rype-mm2 does not fit
NOT_RUN_MM2 = {(16, "long-read"): ("rype-mm2",), (8, "long-read"): ("rype-mm2",)}
FIELDS = ["dataset", "approach", "status", "wall_seconds", "peak_rss_kb", "peak_proc_rss_kb",
          "reads", "server_wall_seconds", "server_peak_rss_kb", "counts_match"]


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def time_v_rss(path):
    if not path.exists():
        return "NA"
    for line in path.read_text().splitlines():
        if "Maximum resident set size" in line:
            return line.rsplit(":", 1)[1].strip()
    return "NA"


def build(mem_gb, mod_dir):
    tier = TIERS[mem_gb]
    limit_kb = mem_gb * 1024 * 1024
    server, reads = {}, {}
    op = REPO / "data" / mod_dir / "op_runtime.tsv"
    if op.exists():
        for r in read_rows(op):
            server.setdefault((r["dataset"], r["approach"]), r)
            reads.setdefault(r["dataset"], r["reads"])
    not_run = tier["not_run"] + NOT_RUN_MM2.get((mem_gb, mod_dir), ())
    run = tuple(a for a in tier["run"] if a not in not_run)
    rows = []
    for dataset in ("holdout", "wol2"):
        for approach in run + not_run:
            srv = server.get((dataset, approach), {})
            row = dict.fromkeys(FIELDS, "NA")
            row.update(dataset=dataset, approach=approach, reads=reads.get(dataset, "NA"),
                       server_wall_seconds=srv.get("wall_seconds", "NA"),
                       server_peak_rss_kb=srv.get("peak_rss_kb", "NA"))
            if approach in not_run:
                row["status"] = "does_not_fit"
                rows.append(row)
                continue
            d = RESULTS / f"laptop-{mem_gb}gb" / mod_dir / dataset / approach
            tr = d / "tree_rss.tsv"
            if not tr.exists():
                row["status"] = "missing"
                rows.append(row)
                continue
            t = read_rows(tr)[0]
            row.update(status="ok" if t["exit_code"] == "0" else f"failed (exit {t['exit_code']})",
                       wall_seconds=t["wall_seconds"], peak_rss_kb=t["peak_tree_rss_kb"],
                       peak_proc_rss_kb=time_v_rss(d / "time-v.stderr"))
            if approach in ("rype", "rype-mm2"):
                ref = RESULTS / mod_dir / dataset / approach / "rep1" / "stats.tsv"
                if (d / "stats.tsv").exists() and ref.exists():
                    row["counts_match"] = str(read_rows(d / "stats.tsv") == read_rows(ref)).lower()
            rows.append(row)
    for r in rows:
        if (r["status"] == "does_not_fit" and r["approach"] in tier["not_run"]
                and r["server_peak_rss_kb"] not in ("NA", "")):
            assert int(r["server_peak_rss_kb"]) > limit_kb, f"{r['approach']} would fit in {mem_gb} GB"
    out = REPO / "data" / mod_dir / f"laptop-{mem_gb}gb_runtime.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mem", type=int, choices=sorted(TIERS), default=16)
    mem = ap.parse_args().mem
    for mod_dir in TIERS[mem]["mods"]:
        build(mem, mod_dir)

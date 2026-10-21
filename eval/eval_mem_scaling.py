#!/usr/bin/env python3
"""Long-read rype memory-scaling table: data/long-read/mem_scaling.tsv.

One row per job memory limit x dataset x replicate (run-mem-scaling.sh):
  status             ok | oom (exit 137: SIGKILL, the cgroup OOM killer) |
                     failed (exit N) | no record (job started, no tree_rss.tsv) |
                     incomplete (exit 0 but no stats.tsv: rype-hf did not run
                     to its end) | missing (never ran)
  wall_seconds       tree_rss.tsv wall time of the timed call
  peak_rss_kb        peak summed RSS of the process tree (tree_rss.tsv)
  limit_kb           the job memory limit
  rype_budget_mib    the --max-memory rype-hf gave rype
  passes             rype classification passes (rype's summary line)
  reads              long reads in the dataset (the 64 GB run's reads_in)
  counts_match       stats.tsv identical to the 64 GB operating-point run
                     (long-read/<set>/rype/rep1): the budget must not
                     change which reads survive

Usage: source common.sh first (it exports the paths), then
       eval_mem_scaling.py
"""
import csv
import os
import re
import sys
from pathlib import Path

def env(name):
    try:
        return os.environ[name]
    except KeyError:
        sys.exit(f"{Path(__file__).name}: {name} is not set; source common.sh first")


REPO, RESULTS = Path(env("REPO")), Path(env("RESULTS_DIR"))
TIERS = (64, 32, 16, 8, 4, 2, 1)
REPS = 3
FIELDS = ["tier_gb", "dataset", "rep", "status", "wall_seconds", "peak_rss_kb", "limit_kb",
          "rype_budget_mib", "passes", "reads", "counts_match"]


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def last_int(pattern, text):
    hits = re.findall(pattern, text)
    return hits[-1] if hits else "NA"


def build():
    rows = []
    for tier in TIERS:
        for dataset in ("holdout", "wol2"):
            ref = RESULTS / "long-read" / dataset / "rype" / "rep1" / "stats.tsv"
            ref_rows = read_rows(ref) if ref.exists() else None
            for rep in range(1, REPS + 1):
                d = RESULTS / "mem-scaling" / "long-read" / dataset / f"{tier}gb" / f"rep{rep}"
                row = dict.fromkeys(FIELDS, "NA")
                row.update(tier_gb=tier, dataset=dataset, rep=rep, limit_kb=tier * 1024 * 1024,
                           reads=ref_rows[0]["reads_in"] if ref_rows else "NA")
                if not d.exists():
                    row["status"] = "missing"
                    rows.append(row)
                    continue
                err = (d / "rype-hf.stderr").read_text(errors="replace") if (d / "rype-hf.stderr").exists() else ""
                row.update(rype_budget_mib=last_int(r"rype --max-memory (\d+) MiB", err),
                           passes=last_int(r"passes=(\d+)", err))
                if not (d / "tree_rss.tsv").exists():
                    row["status"] = "no record"
                    rows.append(row)
                    continue
                t = read_rows(d / "tree_rss.tsv")[0]
                code = int(t["exit_code"])
                row.update(status="ok" if code == 0 else "oom" if code in (137, -9) else f"failed (exit {code})",
                           wall_seconds=t["wall_seconds"], peak_rss_kb=t["peak_tree_rss_kb"])
                if code == 0 and not (d / "stats.tsv").exists():
                    row["status"] = "incomplete"
                elif code == 0 and ref_rows is not None:
                    row["counts_match"] = str(read_rows(d / "stats.tsv") == ref_rows).lower()
                rows.append(row)
    out = REPO / "data" / "long-read" / "mem_scaling.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}", file=sys.stderr)
    for r in rows:
        if r["status"] != "ok" or r["counts_match"] == "false":
            print(f"  {r['tier_gb']} GB {r['dataset']} rep{r['rep']}: {r['status']}"
                  f" counts_match={r['counts_match']}", file=sys.stderr)


if __name__ == "__main__":
    build()

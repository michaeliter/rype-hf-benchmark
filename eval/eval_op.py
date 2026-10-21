#!/usr/bin/env python3
"""Build the figure tables, data/<mod>/, for each modality.

Starts from eval_external.py's tables (${RESULTS_DIR}/<mod>/tables/: read
totals, Movi 2 alone and SOTA) and adds rype, rype-mm2 and both cleanifier
variants from their runs:

  Sweep tables (window and threshold supplementary figures):
    fn.tsv           base, window, approach, fn, total, fnr   (holdout, per haplotype)
    fp.tsv           window, approach, fp, total, fpr          (wol2, pooled)
    runtime.tsv      dataset, approach, window, unit, wall_seconds, reads
    memory.tsv       base, window, approach, peak_rss_kb
      rype / rype-mm2 per window from the window sweep (run-window-sweep.sh;
      w=20 is the operating-point rep1 run); window-independent tools repeat
      their single run at every window.
    storage.tsv      component, window, bytes
    threshold_sweep_w20.tsv  threshold, fn_mean, fnr, fp, fpr: rype alone at
      w=20 over 11 thresholds (run-threshold-sweep.sh), checked against the
      operating-point run at the operating threshold.
    movi_ref.tsv     SOTA FNR / FPR (from eval_external.py)

  Operating-point tables (main benchmarking figure):
    op_runtime.tsv   dataset, approach, rep, wall_seconds, peak_rss_kb, reads
    op_fn.tsv        approach, base, fn, total, fnr
    op_fp.tsv        approach, fp, total, fpr
      rype / rype-mm2: all replicates for timing, replicate 1's survivors for
      FN/FP; every other approach: its single run.

Usage: source common.sh first (it exports the paths), then
       eval_op.py [sr] [lr]      (default: both)
"""
import csv
import os
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import eval_common as E


def env(name):
    try:
        return os.environ[name]
    except KeyError:
        sys.exit(f"eval_op.py: {name} is not set; source common.sh first")


REPO, RESULTS = Path(env("REPO")), Path(env("RESULTS_DIR"))
PIGZ = env("PIGZ_BIN")
# approach -> (storage.tsv component, index files)
CLEANIFIER = {
    "cleanifier": ("cleanifier_index", [Path(env("CLEANIFIER_EXACT_IDX") + ext) for ext in (".hash", ".info")]),
    "cleanifier-filter": ("cleanifier_filter_index", [Path(env("CLEANIFIER_FILTER_IDX") + ext) for ext in (".filter", ".info")]),
}
# rype-mm2's minimap2 index: T2T-CHM13v2.0 alone in both modalities.
RESCUE_MMI = {"sr": Path(env("SR_RESCUE_MMI")), "lr": Path(env("LR_RESCUE_MMI"))}
OP_WINDOW = 20
OP_THRESHOLD = {"sr": "0.0001", "lr": "0.05"}
THR_GRID = ("0.0001", "0.001", "0.005", "0.01", "0.02", "0.05", "0.1", "0.2", "0.3", "0.5", "0.7")
NEW_RYPE = ("rype", "rype-mm2")

MODS = {
    "sr": {"dir": "short-read", "paired": True, "windows": [10, 20, 25, 35, 45]},
    "lr": {"dir": "long-read", "paired": False, "windows": [20, 50, 100, 150, 200]},
}


def read_tsv(path):
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        return rdr.fieldnames, list(rdr)


def write_tsv(path, fields, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def count_by_tag(path, paired):
    """Reads per source tag in a concatenated survivor FASTQ ('@<tag>::<id>').
    paired=True: interleaved pairs, records/2. Missing file -> {} (0 survivors)."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    awk = "NR % 4 == 1 { split($0, a, \"::\"); c[substr(a[1], 2)]++ } END { for (k in c) print k \"\\t\" c[k] }"
    out = subprocess.run(f"{PIGZ} -dc -- '{path}' | awk '{awk}'", shell=True, check=True,
                         capture_output=True, text=True).stdout
    counts = {}
    for line in out.splitlines():
        tag, n = line.split("\t")
        counts[tag] = int(n) // 2 if paired else int(n)
    return counts


def survivor_tags(run_dir, paired):
    """Per-tag survivor reads (pairs for short reads) of a rype-hf run: from
    survivors_R1.fastq.gz for split R1/R2 output (one record per pair), else
    from the single survivors.fastq.gz."""
    r1 = Path(run_dir) / "survivors_R1.fastq.gz"
    if r1.exists():
        return count_by_tag(r1, paired=False)
    return count_by_tag(Path(run_dir) / "survivors.fastq.gz", paired)


def fn_fp_rows(approach, holdout_counts, wol2_counts, base_totals, wol2_total):
    fn = [{"approach": approach, "base": b, "fn": holdout_counts.get(b, 0), "total": t,
           "fnr": f"{holdout_counts.get(b, 0) / t:.12f}"} for b, t in base_totals.items()]
    fp_n = wol2_total - sum(wol2_counts.values())
    fp = {"approach": approach, "fp": fp_n, "total": wol2_total, "fpr": f"{fp_n / wol2_total:.8e}"}
    return fn, fp


def build(mod_key):
    m = MODS[mod_key]
    results = RESULTS / m["dir"]
    tables, data = results / "tables", REPO / "data" / m["dir"]
    data.mkdir(parents=True, exist_ok=True)

    fn_fields, fn_rows = read_tsv(tables / "fn.tsv")
    fp_fields, fp_rows = read_tsv(tables / "fp.tsv")
    rt_fields, rt_rows = read_tsv(tables / "runtime.tsv")
    mem_fields, mem_rows = read_tsv(tables / "memory.tsv")
    st_fields, st_rows = read_tsv(tables / "storage.tsv")
    shutil.copy2(tables / "movi_ref.tsv", data / "movi_ref.tsv")

    _, totals = read_tsv(tables / "totals.tsv")
    base_totals = {r["tag"]: int(r["total"]) for r in totals if r["dataset"] == "holdout"}
    wol2_total = sum(int(r["total"]) for r in totals if r["dataset"] == "wol2")
    reads = {"holdout": sum(base_totals.values()), "wol2": wol2_total}
    rep1_counts = {}

    def survivor_counts(s, a, w):
        """(per-tag survivor counts, wall s, peak tree RSS kB) of a complete
        rype-hf run at window w, else None. w=20 is the operating-point rep1
        run (survivors counted here, cached for the op tables); other windows
        counted their survivors in the job (tag_counts.tsv)."""
        d = results / s / a / ("rep1" if w == OP_WINDOW else f"w{w}")
        if not (d / "tree_rss.tsv").exists():
            return None
        t = next(E.rows(d / "tree_rss.tsv"))
        if t["exit_code"] != "0":
            E.warn(f"{mod_key} {s} {a} w{w}: exit {t['exit_code']}")
            return None
        st = next(E.rows(d / "stats.tsv"))
        if int(st["reads_in"]) != reads[s]:
            E.warn(f"{mod_key} {s} {a} w{w}: reads_in {st['reads_in']} != {reads[s]}")
        if w == OP_WINDOW:
            if (s, a) not in rep1_counts:
                rep1_counts[(s, a)] = survivor_tags(d, m["paired"])
            counts = rep1_counts[(s, a)]
        elif (d / "tag_counts.tsv").exists():
            counts = {r["tag"]: int(r["reads"]) for r in E.rows(d / "tag_counts.tsv")}
        else:
            return None
        return counts, float(t["wall_seconds"]), int(t["peak_tree_rss_kb"])

    # -- cleanifier, exact and probabilistic (arms/arm-cleanifier.sbatch) -------------
    keep = "clean_1_keep.*" if m["paired"] else "clean_keep.*"
    for approach, (component, index_files) in CLEANIFIER.items():
        res = {}
        for s in ("holdout", "wol2"):
            d = results / s / approach
            if not (d / "tree_rss.tsv").exists():
                continue
            t = next(E.rows(d / "tree_rss.tsv"))
            kept = sorted(d.glob(keep))
            if t["exit_code"] == "0" and kept:
                # R1-only keep file: one record per pair already, so paired=False.
                res[s] = (count_by_tag(kept[0], paired=False), float(t["wall_seconds"]),
                          int(t["peak_tree_rss_kb"]))
        if len(res) == 2:
            (h_counts, h_wall, h_rss), (w_counts, w_wall, _) = res["holdout"], res["wol2"]
            fp_n = wol2_total - sum(w_counts.values())
            fn_rows += [{"base": b, "window": w, "approach": approach, "fn": h_counts.get(b, 0),
                         "total": t, "fnr": f"{h_counts.get(b, 0) / t:.12f}"}
                        for w in m["windows"] for b, t in base_totals.items()]
            fp_rows += [{"window": w, "approach": approach, "fp": fp_n, "total": wol2_total,
                         "fpr": f"{fp_n / wol2_total:.8e}"} for w in m["windows"]]
            rt_rows += [{"dataset": ds, "approach": approach, "window": w,
                         "unit": "base" if ds == "holdout" else "chunk", "wall_seconds": f"{wall:.2f}",
                         "reads": reads[ds]} for ds, wall in (("holdout", h_wall), ("wol2", w_wall))
                        for w in m["windows"]]
            mem_rows += [{"base": "ALL", "window": w, "approach": approach, "peak_rss_kb": h_rss}
                         for w in m["windows"]]
        else:
            E.warn(f"{mod_key}: {approach} results incomplete ({sorted(res)}), left out")
        if all(p.exists() for p in index_files):
            st_rows.append({"component": component, "window": "NA",
                            "bytes": sum(p.stat().st_size for p in index_files)})
        else:
            E.warn(f"{mod_key}: {approach} index files missing, no storage row")

    # -- rype / rype-mm2 window sweep (run-window-sweep.sh) ------------------------------
    for a in NEW_RYPE:
        runs = {(s, w): survivor_counts(s, a, w) for s in ("holdout", "wol2") for w in m["windows"]}
        missing = sorted(f"{s} w{w}" for (s, w), r in runs.items() if r is None)
        if missing:
            E.warn(f"{mod_key} {a}: window sweep incomplete ({', '.join(missing)}), left out of the sweep tables")
            continue
        for w in m["windows"]:
            (h_counts, h_wall, h_rss), (w_counts, w_wall, _) = runs[("holdout", w)], runs[("wol2", w)]
            fp_n = wol2_total - sum(w_counts.values())
            fn_rows += [{"base": b, "window": w, "approach": a, "fn": h_counts.get(b, 0),
                         "total": t, "fnr": f"{h_counts.get(b, 0) / t:.12f}"} for b, t in base_totals.items()]
            fp_rows.append({"window": w, "approach": a, "fp": fp_n, "total": wol2_total,
                            "fpr": f"{fp_n / wol2_total:.8e}"})
            rt_rows += [{"dataset": ds, "approach": a, "window": w, "unit": unit,
                         "wall_seconds": f"{wall:.2f}", "reads": reads[ds]}
                        for ds, unit, wall in (("holdout", "base", h_wall), ("wol2", "chunk", w_wall))]
            mem_rows.append({"base": "ALL", "window": w, "approach": a, "peak_rss_kb": h_rss})

    # -- rype threshold sweep at w=20 (run-threshold-sweep.sh): host reads per
    # threshold and source tag, from one `rype classify run` per set ------------------
    host = {}  # (set, threshold) -> {tag: host reads}
    done = []
    for s in ("holdout", "wol2"):
        f = results / s / "threshold-sweep" / "host_counts.tsv"
        if f.exists():
            done.append(s)
            for r in E.rows(f):
                host.setdefault((s, r["threshold"]), {})[r["tag"]] = int(r["host"])
    if len(done) < 2:
        E.warn(f"{mod_key}: threshold sweep incomplete (done: {done}), no threshold_sweep_w20.tsv")
    else:
        # The operating-point threshold must reproduce the rype-hf runs exactly.
        op_thr = OP_THRESHOLD[mod_key]
        for s in ("holdout", "wol2"):
            st = next(E.rows(results / s / "rype" / "rep1" / "stats.tsv"))
            removed = int(st["reads_in"]) - int(st["rype_kept"])
            swept = sum(host.get((s, op_thr), {}).values())
            if swept != removed:
                E.warn(f"{mod_key} {s}: threshold sweep host at t={op_thr} is {swept}, "
                       f"rype-hf rep1 removed {removed}")
        with open(data / "threshold_sweep_w20.tsv", "w") as fh:
            fh.write("threshold\tfn_mean\tfnr\tfp\tfpr\n")
            for t in THR_GRID:
                h = host.get(("holdout", t), {})
                fns = [(tot - h.get(b, 0), tot) for b, tot in base_totals.items()]
                fn_mean = statistics.mean(fn for fn, _ in fns)
                fnr = statistics.mean(fn / tot for fn, tot in fns)
                fp_n = sum(host.get(("wol2", t), {}).values())
                fh.write(f"{t}\t{fn_mean:.2f}\t{fnr:.6e}\t{fp_n}\t{fp_n / wol2_total:.6e}\n")

    write_tsv(data / "fn.tsv", fn_fields, sorted(fn_rows, key=lambda r: (r["base"], int(r["window"]), r["approach"])))
    write_tsv(data / "fp.tsv", fp_fields, sorted(fp_rows, key=lambda r: (int(r["window"]), r["approach"])))
    write_tsv(data / "runtime.tsv", rt_fields, sorted(rt_rows, key=lambda r: (r["dataset"], r["approach"], int(r["window"]))))
    write_tsv(data / "memory.tsv", mem_fields, sorted(mem_rows, key=lambda r: (int(r["window"]), r["approach"])))
    st_rows.append({"component": "rype_mm2_mmi", "window": "NA", "bytes": RESCUE_MMI[mod_key].stat().st_size})
    write_tsv(data / "storage.tsv", st_fields, st_rows)

    # -- operating-point tables ------------------------------------------------------
    op_rt, op_fn, op_fp = [], [], []
    op_mem = {r["approach"]: int(r["peak_rss_kb"]) for r in mem_rows if int(r["window"]) == OP_WINDOW}
    for a in sorted({r["approach"] for r in rt_rows} - set(NEW_RYPE)):
        for r in rt_rows:
            if r["approach"] == a and int(r["window"]) == OP_WINDOW:
                op_rt.append({"dataset": r["dataset"], "approach": a, "rep": 1,
                              "wall_seconds": r["wall_seconds"],
                              "peak_rss_kb": op_mem[a] if r["dataset"] == "holdout" else "NA",
                              "reads": r["reads"]})
        op_fn += [{"approach": a, "base": r["base"], "fn": r["fn"], "total": r["total"], "fnr": r["fnr"]}
                  for r in fn_rows if r["approach"] == a and int(r["window"]) == OP_WINDOW]
        op_fp += [{"approach": a, "fp": r["fp"], "total": r["total"], "fpr": r["fpr"]}
                  for r in fp_rows if r["approach"] == a and int(r["window"]) == OP_WINDOW]

    for a in NEW_RYPE:
        n_reps = {}
        for s in ("holdout", "wol2"):
            n_reps[s] = 0
            for rep_dir in sorted((results / s / a).glob("rep*")):
                tr = rep_dir / "tree_rss.tsv"
                if not rep_dir.name.removeprefix("rep").isdigit() or not tr.exists():
                    continue
                row = next(E.rows(tr))
                if int(row["exit_code"]) != 0:
                    E.warn(f"{mod_key} {s} {a} {rep_dir.name}: exit {row['exit_code']}, skipped")
                    continue
                st = next(E.rows(rep_dir / "stats.tsv"))
                if int(st["reads_in"]) != reads[s]:
                    E.warn(f"{mod_key} {s} {a} {rep_dir.name}: reads_in {st['reads_in']} != {reads[s]}")
                n_reps[s] += 1
                op_rt.append({"dataset": s, "approach": a, "rep": rep_dir.name.removeprefix("rep"),
                              "wall_seconds": row["wall_seconds"], "peak_rss_kb": row["peak_tree_rss_kb"],
                              "reads": reads[s]})
        if not n_reps["holdout"] or not n_reps["wol2"]:
            E.warn(f"{mod_key} {a}: operating-point runs incomplete, left out")
            continue
        counts = {s: rep1_counts[(s, a)] if (s, a) in rep1_counts
                  else survivor_tags(results / s / a / "rep1", m["paired"])
                  for s in ("holdout", "wol2")}
        fn, fp = fn_fp_rows(a, counts["holdout"], counts["wol2"], base_totals, wol2_total)
        op_fn += fn
        op_fp.append(fp)
        E.warn(f"{mod_key} {a}: FN {sum(r['fn'] for r in fn)}, FP {fp['fp']}")

    write_tsv(data / "op_runtime.tsv", ["dataset", "approach", "rep", "wall_seconds", "peak_rss_kb", "reads"], op_rt)
    write_tsv(data / "op_fn.tsv", ["approach", "base", "fn", "total", "fnr"], op_fn)
    write_tsv(data / "op_fp.tsv", ["approach", "fp", "total", "fpr"], op_fp)
    print(f"[{mod_key}] wrote {data}")


if __name__ == "__main__":
    for k in (sys.argv[1:] or ["sr", "lr"]):
        build(k)

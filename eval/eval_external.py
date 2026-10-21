#!/usr/bin/env python3
"""Read totals and the Movi 2 / SOTA tables, one modality at a time.

Counts the QC'd reads of every holdout haplotype and WoL2 chunk (the
denominators of every FN/FP rate) and the survivors of the two
window-independent arms that are split back per source tag:
  mm2-mm2-movi   arms/arm-sota.sbatch        (<mod>/<set>/sota/)
  movi2-alone    arms/arm-movi2-only.sbatch  (<mod>/<set>/movi2-only/)

FN (holdout, every read is host) = survivors; FP (wol2, every read is
non-host) = total - survivors, pooled over chunks. Short reads count pairs:
the survivor files are interleaved, so records // 2 (Movi filters per read and
can leave an orphan mate, hence the rounding down).

Writes ${RESULTS_DIR}/<mod>/tables/, which eval_op.py builds on:
  totals.tsv   dataset, tag, total        (reads, or pairs for short reads)
  fn.tsv       base, window, approach, fn, total, fnr   (repeated per window,
  fp.tsv       window, approach, fp, total, fpr          so they line up with
  runtime.tsv  dataset, approach, window, unit, wall_seconds, reads   the
  memory.tsv   base, window, approach, peak_rss_kb       window sweep)
  storage.tsv  component, window, bytes   (rype indices, SOTA minimap2
                                           indices, Movi 2 index)
  movi_ref.tsv metric, count, total, rate (SOTA FNR / FPR)
Runtime is the sum of an arm's timed stages, memory their maximum.

Usage: source common.sh first (it exports the paths), then
       eval_external.py <sr|lr>
"""
import os
import statistics
import sys
from pathlib import Path

import eval_common as E

MODS = {
    "sr": {"dir": "short-read", "paired": True, "windows": [10, 20, 25, 35, 45],
           "sota_mmi": [("sota_hg38_mmi", "QIITA_HG38_MMI"), ("sota_t2t_mmi", "QIITA_T2T_MMI")]},
    "lr": {"dir": "long-read", "paired": False, "windows": [20, 50, 100, 150, 200],
           "sota_mmi": [("sota_grch38_mmi", "LR_SOTA_HG38_MMI"), ("sota_t2t_mmi", "LR_RESCUE_MMI")]},
}
ARM_DIR = {"mm2-mm2-movi": "sota", "movi2-alone": "movi2-only"}


def env(name):
    try:
        return os.environ[name]
    except KeyError:
        sys.exit(f"eval_external.py: {name} is not set; source common.sh first")


def build(mod_key):
    m = MODS[mod_key]
    repo, reads_dir = Path(env("REPO")), Path(env("READS_DIR"))
    results = Path(env("RESULTS_DIR")) / m["dir"]
    out = results / "tables"
    out.mkdir(parents=True, exist_ok=True)
    holdout = [line.split("\t")[0] for line in (repo / "references" / "holdout10.tsv").read_text().splitlines()[1:]]
    qc = reads_dir / m["dir"]

    # -- denominators ---------------------------------------------------------------
    r1 = "_R1.fastq.gz" if m["paired"] else ".fastq.gz"
    hold_files = {b: qc / "holdout" / "qc" / f"{b}{r1}" for b in holdout}
    counts = E.count_many(hold_files.values())  # R1 records = pairs
    hold_total = {b: counts[str(f)] for b, f in hold_files.items()}
    chunk_lists = sorted((repo / "reads" / "wol2-chunks" / m["dir"]).glob("chunk_*.list"))
    wol2_total = E.chunk_totals_parallel(chunk_lists, qc / "wol2" / "qc", paired=m["paired"])
    with open(out / "totals.tsv", "w") as fh:
        fh.write("dataset\ttag\ttotal\n")
        for b, t in hold_total.items():
            fh.write(f"holdout\t{b}\t{t}\n")
        for c, t in wol2_total.items():
            fh.write(f"wol2\t{c}\t{t}\n")
    E.warn(f"{mod_key}: {len(hold_total)} holdout haplotypes, {len(wol2_total)} wol2 chunks")

    # -- survivors ------------------------------------------------------------------
    def surv(s, a, tag):
        return results / s / ARM_DIR[a] / "survivors" / f"{tag}.survivors.fastq.gz"

    paths = [surv("holdout", a, b) for a in ARM_DIR for b in hold_total] + \
            [surv("wol2", a, c) for a in ARM_DIR for c in wol2_total]
    # A missing survivor file means no read of that tag survived (split_survivors.py
    # only writes tags it sees), so .get(path, 0).
    n = E.count_many(paths, paired=m["paired"])
    fn = {a: {b: n.get(str(surv("holdout", a, b)), 0) for b in hold_total} for a in ARM_DIR}
    fp = {a: sum(t - n.get(str(surv("wol2", a, c)), 0) for c, t in wol2_total.items()) for a in ARM_DIR}
    tot_w = sum(wol2_total.values())

    with open(out / "fn.tsv", "w") as fh:
        fh.write("base\twindow\tapproach\tfn\ttotal\tfnr\n")
        for b, w, a in sorted((b, w, a) for a in ARM_DIR for b in hold_total for w in m["windows"]):
            fh.write(f"{b}\t{w}\t{a}\t{fn[a][b]}\t{hold_total[b]}\t{fn[a][b] / hold_total[b]:.12f}\n")
    with open(out / "fp.tsv", "w") as fh:
        fh.write("window\tapproach\tfp\ttotal\tfpr\n")
        for w in m["windows"]:
            for a in sorted(ARM_DIR):
                fh.write(f"{w}\t{a}\t{fp[a]}\t{tot_w}\t{fp[a] / tot_w:.8e}\n")
    sota_fn, sota_tot = sum(fn["mm2-mm2-movi"].values()), sum(hold_total.values())
    with open(out / "movi_ref.tsv", "w") as fh:
        fh.write("metric\tcount\ttotal\trate\n")
        fh.write(f"fnr\t{sota_fn}\t{sota_tot}\t{sota_fn / sota_tot:.6e}\n")
        fh.write(f"fpr\t{fp['mm2-mm2-movi']}\t{tot_w}\t{fp['mm2-mm2-movi'] / tot_w:.6e}\n")

    # -- runtime and memory ------------------------------------------------------------
    rt, mem = [], []
    for ds, totals, unit in (("holdout", hold_total, "base"), ("wol2", wol2_total, "chunk")):
        reads = sum(totals.values())
        for a in ARM_DIR:
            for base, (wall, rss) in E.sota_resource_totals(results / ds / ARM_DIR[a] / "resource").items():
                for w in m["windows"]:
                    rt.append((ds, a, w, unit, wall, reads))
                    if ds == "holdout":
                        mem.append((base, w, a, rss))
    with open(out / "runtime.tsv", "w") as fh:
        fh.write("dataset\tapproach\twindow\tunit\twall_seconds\treads\n")
        for ds, a, w, u, s, r in sorted(rt, key=lambda x: (x[0], x[1], x[2], x[3])):
            fh.write(f"{ds}\t{a}\t{w}\t{u}\t{s:.2f}\t{r}\n")
    with open(out / "memory.tsv", "w") as fh:
        fh.write("base\twindow\tapproach\tpeak_rss_kb\n")
        for b, w, a, k in sorted(mem):
            fh.write(f"{b}\t{w}\t{a}\t{k}\n")

    # -- storage --------------------------------------------------------------------
    sizes = [("rype_ryxdi", w, Path(env("INDEX_DIR")) / "rype" / f"human-all97-w{w}.ryxdi") for w in m["windows"]]
    sizes += [(label, "NA", Path(env(var))) for label, var in m["sota_mmi"]]
    sizes.append(("movi2_index", "NA", Path(env("MOVI_INDEX"))))
    with open(out / "storage.tsv", "w") as fh:
        fh.write("component\twindow\tbytes\n")
        for label, w, path in sizes:
            b = E.du_bytes(path) if path.exists() else None
            if b is None:
                E.warn(f"{mod_key}: {path} missing, no storage row")
                continue
            fh.write(f"{label}\t{w}\t{b}\n")
    E.warn(f"{mod_key}: wrote {out}")
    for a in ARM_DIR:
        E.warn(f"  {a}: FN {sum(fn[a].values())} (mean FNR "
               f"{statistics.mean(fn[a][b] / hold_total[b] for b in hold_total):.3e}), FP {fp[a]}")


if __name__ == "__main__":
    for k in (sys.argv[1:] or ["sr", "lr"]):
        build(k)

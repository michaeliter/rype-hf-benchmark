"""Shared counting helpers for eval_external.py and eval_op.py.

FN/FP are counted from survivor files: every arm writes the reads it keeps.
Interleaved short-read output holds both mates of a pair back to back, one
record each, so a pair is 2 records; a long read is 1 record.
"""
import csv
import gzip
import os
import shlex
import statistics
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

DEFAULT_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 4)


def warn(msg):
    import sys
    print(f"[eval] {msg}", file=sys.stderr)


def count_fastq_records(path):
    """Number of sequence records in a plain or .gz file, FASTQ (4 lines/record)
    or FASTA (2 lines/record) -- auto-detected from the first byte ('@' vs '>').
    0 if missing. movi's --filter emits FASTA despite arm-sota.sbatch naming its
    output *.fastq.gz, so this has to handle both.

    Streams through `wc -l` rather than materializing the decompressed file in
    Python (the original `zcat ... | .stdout` approach buffered the ENTIRE
    decompressed FASTQ as a bytes object -- fine serially, but count_many()
    runs many of these concurrently via a thread pool, and N full-file buffers
    alive at once OOM'd the eval job at -c16/32gb on real WGS-scale FASTQs).
    zcat and wc both stream page-by-page with a small fixed buffer, so memory
    stays flat regardless of file size or concurrency."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return 0
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            first_byte = fh.read(1)
        if not first_byte:
            return 0
        wc = subprocess.run(["bash", "-c", f"zcat -- {shlex.quote(str(path))} | wc -l"],
                             stdout=subprocess.PIPE, check=True)
    else:
        with open(path, "rb") as fh:
            first_byte = fh.read(1)
        if not first_byte:
            return 0
        wc = subprocess.run(["wc", "-l", "--", str(path)], stdout=subprocess.PIPE, check=True)
    n = int(wc.stdout.split()[0])
    lines_per_record = 2 if first_byte == b">" else 4
    assert n % lines_per_record == 0, f"{path}: {n} lines not a multiple of {lines_per_record}"
    return n // lines_per_record


def count_reads(path, paired):
    """Records -> reads: pairs (records/2) for paired interleaved output, else records.

    Movi filters per read, not per pair, so its survivor file can hold an
    orphaned mate (odd record count): round down rather than assert."""
    n = count_fastq_records(path)
    if paired:
        return n // 2
    return n


def count_many(paths, paired=False, workers=None):
    """Parallel count over many files -- each call is one subprocess (zcat)
    spent mostly blocked waiting on decompression, so a thread pool (GIL
    released during subprocess.run) gives near-linear wall-clock speedup
    over a serial loop. Missing/empty paths are simply omitted from the
    result dict (not zeroed). Use `.get(path, 0)` at call sites, not
    membership-as-skip: split_survivors.py only creates a file for a tag it
    actually sees a record for, so in this concatenated-input pipeline a
    missing survivor file means "genuinely 0 survivors" (once every arm has
    confirmed completed), not "didn't run yet" -- treating it as skip
    silently drops perfect-filtering (0-survivor) results.

    paired=True applies count_reads' interleaved pairs/2 division (for
    interleaved survivor files, R1 and R2 one record each back to back).
    paired=False (the default) returns
    raw record counts via count_fastq_records -- use this for R1-only
    per-mate files (e.g. wol2 QC accessions, see chunk_totals_parallel)
    where one record already equals one pair; applying the /2 there would
    silently halve the total.

    Returns {str(path): count}.
    """
    existing = [Path(p) for p in paths if Path(p).exists() and Path(p).stat().st_size > 0]
    if not existing:
        return {}
    fn = (lambda p: count_reads(p, True)) if paired else count_fastq_records
    with ThreadPoolExecutor(max_workers=workers or DEFAULT_WORKERS) as ex:
        counts = list(ex.map(fn, existing))
    return {str(p): c for p, c in zip(existing, counts)}


def chunk_totals_parallel(chunk_lists, qc_dir, paired, workers=None):
    """Parallel equivalent of calling chunk_total() once per chunk list --
    resolves every (chunk, accession) file path across ALL chunks up front,
    counts them in one shared thread pool, then sums per chunk. wol2 chunk
    lists run into the thousands of accessions combined, each needing its
    own zcat subprocess; doing that serially (the original chunk_total, one
    chunk at a time) is the single biggest wall-clock cost of the eval.

    Returns {chunk_id: total}, chunk_id = list_path.stem (matches the
    original chunk_total call sites' `cid = lf.stem` convention).
    """
    per_chunk_paths = {}
    for lf in chunk_lists:
        cid = Path(lf).stem
        paths = []
        for acc in Path(lf).read_text().split():
            if paired:
                paths.append(qc_dir / f"{acc}_R1.fastq.gz")
            else:
                p = Path(acc)
                paths.append(p if p.is_absolute() else qc_dir / f"{acc}.fastq.gz")
        per_chunk_paths[cid] = paths

    all_paths = [p for paths in per_chunk_paths.values() for p in paths]
    counts = count_many(all_paths, paired=False, workers=workers)  # R1-only files: 1 record = 1 pair already
    return {cid: sum(counts.get(str(p), 0) for p in paths)
            for cid, paths in per_chunk_paths.items()}


def parse_wall(s):
    """'/usr/bin/time -v' Elapsed (wall clock) time -> seconds. 'h:mm:ss' or 'm:ss.ss'."""
    if not s:
        return float("nan")
    parts = s.strip().split(":")
    parts = [float(p) for p in parts]
    if len(parts) == 3:
        h, m, sec = parts
    elif len(parts) == 2:
        h = 0.0; m, sec = parts
    else:
        return float("nan")
    return h * 3600 + m * 60 + sec


def rows(path):
    with open(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def sota_resource_totals(resource_dir):
    """sota's resource/<base>.tsv has one row per stage (align_hg38/align_t2t/
    movi2_filter, or GRCh38/T2T/movi2 for LR) -- sum wall_time, max peak_rss_kb,
    per base. Returns {base: (wall_seconds, peak_rss_kb)}."""
    out = {}
    for f in Path(resource_dir).glob("*.tsv"):
        wall = 0.0
        rss = 0
        base = None
        for r in rows(f):
            base = r["base"]
            wall += parse_wall(r["wall_time"])
            rss = max(rss, int(r["peak_rss_kb"] or 0))
        if base:
            out[base] = (wall, rss)
    return out


def du_bytes(path):
    try:
        return int(subprocess.check_output(["du", "-sbL", str(path)]).split()[0])
    except Exception:
        return None

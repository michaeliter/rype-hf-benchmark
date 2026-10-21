#!/usr/bin/env python3
"""Split a survivor FASTQ (or FASTA) produced from a concatenated,
ID-prefixed input (reads/build-concat-inputs.sbatch) back into one file per
source tag (holdout haplotype or WoL2 chunk), so FN/FP can be counted per tag.

Each record's header is "@<tag>::<original_id> ..." (FASTQ, 4 lines/record)
or ">..." (FASTA, 2 lines/record), auto-detected from the first byte: Movi's
--filter writes FASTA although the Movi arms name their output
*.survivors.fastq.gz. Records are grouped by <tag>, "<tag>::" is stripped
(restoring the original ID), and each group goes to <out_dir>/<tag><suffix>.

Compression goes through pigz (PIGZ_BIN, from config.sh) when available,
falling back to Python's gzip module, which is far slower on wol2-sized files.

Usage: split_survivors.py <concat_survivors.fastq.gz> <out_dir> <suffix>
"""
import gzip
import io
import os
import shutil
import subprocess
import sys
from pathlib import Path

PIGZ_BIN = os.environ.get("PIGZ_BIN") or shutil.which("pigz")
if PIGZ_BIN and not (os.path.isfile(PIGZ_BIN) and os.access(PIGZ_BIN, os.X_OK)):
    PIGZ_BIN = None
# Capped low (not pigz's auto-detected nproc) because up to ~16 tag writers
# run concurrently within one job's CPU allocation -- an uncapped -p would
# oversubscribe cores badly (16 tags x nproc threads each).
WRITER_THREADS = 2


def open_read_stream(path):
    """Returns (proc_or_None, text_file). proc is None when falling back to
    the gzip module (nothing to wait() on); non-None when pigz is used."""
    if PIGZ_BIN:
        proc = subprocess.Popen([PIGZ_BIN, "-dc", str(path)], stdout=subprocess.PIPE)
        return proc, io.TextIOWrapper(proc.stdout, encoding="utf-8")
    return None, gzip.open(path, "rt")


def close_read_stream(proc, fh):
    fh.close()
    if proc is not None:
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f"pigz -dc exited {rc}")


def open_write_stream(path):
    """Returns (proc_or_None, text_file, raw_file). raw_file is the
    underlying binary file handle pigz's stdout is wired to; kept around so
    it can be closed once the child has its own duplicated fd."""
    raw_f = open(path, "wb")
    if PIGZ_BIN:
        proc = subprocess.Popen(
            [PIGZ_BIN, "-c", "-p", str(WRITER_THREADS)],
            stdin=subprocess.PIPE, stdout=raw_f,
        )
        raw_f.close()  # child has its own dup'd fd now
        return proc, io.TextIOWrapper(proc.stdin, encoding="utf-8"), None
    return None, gzip.open(raw_f, "wt", compresslevel=6), raw_f


def close_write_stream(proc, fh, raw_f):
    fh.close()
    if proc is not None:
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f"pigz -c exited {rc}")
    elif raw_f is not None:
        raw_f.close()


def main():
    if len(sys.argv) != 4:
        sys.exit(f"usage: {sys.argv[0]} <concat_survivors.fastq.gz> <out_dir> <suffix>")
    src, out_dir, suffix = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
    out_dir.mkdir(parents=True, exist_ok=True)

    # Plain gzip module here, not open_read_stream/pigz -- a 1-byte peek gets
    # no benefit from parallel decompression, and deliberately stopping after
    # 1 byte would leave a pigz subprocess's write end of the pipe closed
    # under it, killing it with SIGPIPE (a real crash hit in production: the
    # subprocess's exit code -13 was wrongly treated as fatal).
    with gzip.open(src, "rb") as peek_fh:
        first_byte = peek_fh.read(1)
    if not first_byte:
        return  # empty input -- nothing to split
    is_fasta = first_byte == b">"
    lines_per_record = 2 if is_fasta else 4

    handles = {}
    counts = {}
    in_proc, fh = open_read_stream(src)
    try:
        while True:
            header = fh.readline()
            if not header:
                break
            lines = [fh.readline() for _ in range(lines_per_record - 1)]
            if not lines[-1]:
                raise ValueError(f"{src}: truncated record after {header!r}")

            assert header[0] in ("@", ">"), f"{src}: bad header {header!r}"
            tag, sep, rest = header[1:].partition("::")
            if not sep:
                raise ValueError(f"{src}: header missing tag delimiter '::': {header!r}")

            if tag not in handles:
                path = out_dir / f"{tag}{suffix}"
                tmp_path = path.with_suffix(path.suffix + ".tmp")
                w_proc, w_fh, w_raw = open_write_stream(tmp_path)
                handles[tag] = (w_proc, w_fh, w_raw, path, tmp_path)
                counts[tag] = 0

            _, fh_out, _, _, _ = handles[tag]
            fh_out.write(header[0] + rest)
            for line in lines:
                fh_out.write(line)
            counts[tag] += 1
    finally:
        close_read_stream(in_proc, fh)
        for tag, (w_proc, w_fh, w_raw, path, tmp_path) in handles.items():
            close_write_stream(w_proc, w_fh, w_raw)
            tmp_path.rename(path)
            print(f"[split_survivors] {tag}: {counts[tag]} records -> {path}", file=sys.stderr)


if __name__ == "__main__":
    main()

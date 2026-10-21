#!/usr/bin/env python3
"""Run a command; record its wall time and the peak SUMMED RSS of its whole
process tree.

`/usr/bin/time -v` reports the max RSS of the single largest process, which
under-counts a pipeline: rype-hf runs rype, minimap2, samtools and bgzip
concurrently, and their footprints add up. This samples
VmRSS of every descendant every INTERVAL seconds and keeps the largest total.

Usage: tree_rss.py OUT.tsv -- cmd [args...]
Writes one row: wall_seconds, peak_tree_rss_kb, exit_code. Exits with the
command's exit code.
"""
import os
import subprocess
import sys
import time

INTERVAL = 0.5


def children_map():
    kids = {}
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        try:
            with open(f"/proc/{name}/stat") as fh:
                stat = fh.read()
        except OSError:
            continue
        # ppid is the 2nd field after the parenthesised comm (which may hold spaces)
        ppid = int(stat.rsplit(")", 1)[1].split()[1])
        kids.setdefault(ppid, []).append(int(name))
    return kids


def tree_rss_kb(root):
    kids, total, stack = children_map(), 0, [root]
    while stack:
        pid = stack.pop()
        stack.extend(kids.get(pid, ()))
        try:
            with open(f"/proc/{pid}/status") as fh:
                for line in fh:
                    if line.startswith("VmRSS:"):
                        total += int(line.split()[1])
                        break
        except OSError:
            pass
    return total


def main():
    if len(sys.argv) < 4 or sys.argv[2] != "--":
        sys.exit(__doc__)
    out, cmd = sys.argv[1], sys.argv[3:]
    t0 = time.monotonic()
    proc = subprocess.Popen(cmd)
    peak = 0
    while proc.poll() is None:
        peak = max(peak, tree_rss_kb(proc.pid))
        time.sleep(INTERVAL)
    wall = time.monotonic() - t0
    with open(out, "w") as fh:
        fh.write("wall_seconds\tpeak_tree_rss_kb\texit_code\n")
        fh.write(f"{wall:.2f}\t{peak}\t{proc.returncode}\n")
    sys.exit(proc.returncode)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""BED of the chunks --ltr-timeout cut short, parsed from an LTRquest run log.

Reads the detector's [WARN] lines ("ltrharvest timeout for index .../ltr_runs/<chunk>/idx",
"ltr_finder timeout ... on .../ltr_runs/<chunk>/<chunk>...fa"), keeps the round from the
work-dir prefix (<prefix>_r<N>.work), and writes chrom start end tool round (0-based).
usage: capped_chunks.py run.log > capped.bed
"""
import re
import sys

pat_lh = re.compile(r"ltrharvest timeout for index (\S+?)_r(\d+)\.work/ltr_runs/([^/]+)/idx")
pat_lf = re.compile(r"ltr_finder timeout after \d+s on (\S+?)_r(\d+)\.work/ltr_runs/([^/]+)/")
pat_fl = re.compile(r"FastLTR timeout after \d+s on (\S+?)_r(\d+)\.work/ltr_runs/([^/]+)/")
chunk = re.compile(r"^(.+)\.chunk\d+\.(\d+)-(\d+)$")
seen = set()
for line in open(sys.argv[1], errors="replace"):
    for tool, pat in (("LTRharvest", pat_lh), ("LTR_FINDER", pat_lf), ("FastLTR", pat_fl)):
        for m in pat.finditer(line):
            c = chunk.match(m.group(3))
            if not c:
                continue
            key = (c.group(1), int(c.group(2)) - 1, int(c.group(3)), tool, int(m.group(2)))
            if key not in seen:
                seen.add(key)
for k in sorted(seen):
    print("\t".join(map(str, k)))

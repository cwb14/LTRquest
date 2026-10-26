#!/usr/bin/env python3
"""FasTAN .1aln (via ONEview) -> BED of tandem arrays, merged, with a unit class.

Each FasTAN alignment spans [b_start, a_end) of one array at period `unit`.
Overlapping/abutting (<= GAP bp) alignments of the same class are merged.
Classes by unit: micro (< 10 bp, includes telomeric repeats), satellite
(10 - 999 bp), long (>= 1000 bp, e.g. rDNA units). Output: chrom start end
class median_unit n_alignments (0-based half-open), arrays >= --min-len only.
usage: fastan_arrays.py ONEview in.1aln [--min-len 10000] [--gap 100] > arrays.bed
"""
import argparse
import subprocess
from collections import defaultdict
from statistics import median


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('oneview')
    ap.add_argument('aln')
    ap.add_argument('--min-len', type=int, default=10000)
    ap.add_argument('--gap', type=int, default=100)
    a = ap.parse_args()
    txt = subprocess.run([a.oneview, a.aln], capture_output=True, text=True, check=True).stdout
    contigs, scaf, off = [], None, 0
    rows, cur = [], None
    for line in txt.splitlines():
        f = line.split()
        if not f:
            continue
        t = f[0]
        if t == 'S':
            scaf, off = f[2], 0
        elif t == 'C':
            contigs.append((scaf, off))
            off += int(f[1])
        elif t == 'G':
            off += int(f[1])
        elif t == 'A':
            if cur:
                rows.append(cur)
            ac, as_, ae, bc, bs, be = map(int, f[1:7])
            name, o = contigs[ac]
            cur = [name, bs + o, ae + o, 0]
        elif t == 'U' and cur:
            cur[3] = int(f[1])
    if cur:
        rows.append(cur)
    cls = lambda u: 'micro' if u < 10 else ('satellite' if u < 1000 else 'long')  # noqa: E731
    by = defaultdict(list)
    for c, s, e, u in rows:
        by[(c, cls(u))].append((s, e, u))
    out = []
    for (c, k), v in by.items():
        v.sort()
        ms, me, us = v[0][0], v[0][1], [v[0][2]]
        for s, e, u in v[1:]:
            if s <= me + a.gap:
                me = max(me, e)
                us.append(u)
            else:
                out.append((c, ms, me, k, int(median(us)), len(us)))
                ms, me, us = s, e, [u]
        out.append((c, ms, me, k, int(median(us)), len(us)))
    for r in sorted(out):
        if r[2] - r[1] >= a.min_len:
            print('\t'.join(map(str, r)))


if __name__ == '__main__':
    main()

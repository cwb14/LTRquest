#!/usr/bin/env python3
"""RepeatMasker support for LTRquest elements in CHM13 (independent annotation).

An element is supported when each of its two LTRs is >= 50% covered by
RepeatMasker records of class LTR and the two LTRs share at least one
RepeatMasker family name (the two LTRs of one insertion are one family).
Writes name<TAB>0/1 for every element in the given clean tables.

usage: rm_support.py --rm chm13v2.0_RepeatMasker.bed --report GCA_..._assembly_report.txt
                     --out support.tsv run1/*_clean_ltr.tsv [...]
"""
import argparse
import bisect
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--rm', required=True)
    ap.add_argument('--report', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('tables', nargs='+')
    a = ap.parse_args()

    ucsc2acc = {}
    for line in open(a.report):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) > 9 and f[4] != 'na':
            ucsc2acc[f[9] if f[9] != 'na' else f[0]] = f[4]
            ucsc2acc[f[0]] = f[4]
            ucsc2acc['chr' + f[0]] = f[4]
    rm = defaultdict(list)
    for line in open(a.rm):
        f = line.split('\t')
        if f[6] != 'LTR':
            continue
        rm[ucsc2acc.get(f[0], f[0])].append((int(f[1]), int(f[2]), f[3]))
    for v in rm.values():
        v.sort()
    starts = {c: [s for s, _, _ in v] for c, v in rm.items()}

    def cover(c, s, e):
        v, k = rm.get(c, []), bisect.bisect_left(starts.get(c, []), s - 50000)
        segs, names = [], set()
        while k < len(v) and v[k][0] < e:
            rs, re_, n = v[k]
            if re_ > s:
                segs.append((max(s, rs), min(e, re_)))
                names.add(n)
            k += 1
        segs.sort()
        tot, cur = 0, None
        for x, y in segs:
            if cur is None or x > cur[1]:
                if cur:
                    tot += cur[1] - cur[0]
                cur = [x, y]
            else:
                cur[1] = max(cur[1], y)
        if cur:
            tot += cur[1] - cur[0]
        return tot / max(1, e - s), names

    with open(a.out, 'w') as out:
        for path in a.tables:
            hdr = None
            for line in open(path):
                if line.startswith('#'):
                    hdr = line[1:].rstrip('\n').split('\t')
                    continue
                f = dict(zip(hdr, line.rstrip('\n').split('\t')))
                chrom, span = f['seq_id'].split('#')[0].rsplit(':', 1)
                s1 = int(span.split('-')[0])
                c5, n5 = cover(chrom, s1 + int(f['ltr5_start']) - 2, s1 + int(f['ltr5_end']) - 1)
                c3, n3 = cover(chrom, s1 + int(f['ltr3_start']) - 2, s1 + int(f['ltr3_end']) - 1)
                ok = int(c5 >= 0.5 and c3 >= 0.5 and bool(n5 & n3))
                out.write(f"{f['seq_id']}\t{ok}\t{c5:.2f},{c3:.2f}\n")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Clean elements inside given regions (e.g. CENH3 centromeres), per run: count, confidence,
superfamily, and curated / RepeatMasker support when a support table is given.

usage: region_counts.py --bed regions.bed [--support support.tsv] RUN_DIR [RUN_DIR ...]
"""
import argparse
import glob
import os
from collections import Counter, defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--bed', required=True)
    ap.add_argument('--support', default=None)
    ap.add_argument('runs', nargs='+')
    a = ap.parse_args()
    reg = defaultdict(list)
    for line in open(a.bed):
        f = line.split()
        reg[f[0]].append((int(f[1]), int(f[2])))
    sup = {}
    if a.support:
        sup = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(a.support)}
    for run in a.runs:
        rows = []
        for p in glob.glob(os.path.join(run, '*_depth*_clean_ltr.tsv')):
            hdr = None
            for line in open(p):
                if line.startswith('#'):
                    hdr = line[1:].rstrip('\n').split('\t')
                    continue
                f = dict(zip(hdr, line.rstrip('\n').split('\t')))
                c, span = f['seq_id'].split('#')[0].rsplit(':', 1)
                s, e = (int(x) for x in span.split('-'))
                if not any(rs < e and s <= re_ for rs, re_ in reg.get(c, [])):
                    continue
                tsd = f.get('tsd', 'NA') not in ('NA', '.', '')
                tg = f.get('motif', '').lower() == 'tg...ca'
                conf = f.get('confidence') or ('high' if tsd and tg else 'medium' if tsd or tg else 'low')
                cls = f['seq_id'].split('#')[1] if '#' in f['seq_id'] else '?'
                rows.append((conf, '/'.join(cls.split('/')[1:3]), sup.get(f['seq_id'])))
        conf = Counter(r[0] for r in rows)
        cls = Counter(r[1] for r in rows).most_common(4)
        known = [r[2] for r in rows if r[2] is not None]
        s = f" supported {sum(known)}/{len(known)}" if known else ""
        print(f"{os.path.basename(run)}: {len(rows)} in regions (high {conf['high']}, medium {conf['medium']}, "
              f"low {conf['low']}){s}; top classes {cls}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Recall of planted LTR-RTs (spike_in.py truth) by each run's clean calls.

A planted element is recovered when one clean call spans it with >= --rec reciprocal
overlap (default 0.9). Recall is split by array class (satellite / long) and by divergence
(planted per-base substitution rate u; LTR-LTR divergence ~2u). The background is a real
assembly, so calls not on a planted element are native elements or errors and are only
counted, not scored.
With --failed-log, planted elements that no round-1 chunk could reach (every chunk holding
them failed in that run, e.g. gt rejecting an 'X') are reported separately and left out of
the second recall figure.
usage: spike_eval.py truth.bed RUN_DIR [RUN_DIR ...] [--failed-log run.log]
"""
import argparse
import glob
import os
from collections import defaultdict


def calls_of(run):
    out = defaultdict(list)
    for p in glob.glob(os.path.join(run, '*_depth*_clean_ltr.tsv')):
        for line in open(p):
            if line.startswith('#') or not line.strip():
                continue
            locus = line.split('\t', 1)[0].split('#')[0]
            c, span = locus.rsplit(':', 1)
            s, e = (int(x) for x in span.split('-'))
            out[c].append((s - 1, e))
    return out


def rec(a, b):
    o = max(0, min(a[1], b[1]) - max(a[0], b[0]))
    return min(o / (a[1] - a[0]), o / (b[1] - b[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('truth')
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--rec', type=float, default=0.9)
    ap.add_argument('--failed-log', default=None)
    a = ap.parse_args()
    truth = [l.rstrip('\n').split('\t') for l in open(a.truth)]
    failed = defaultdict(list)
    if a.failed_log:
        import re
        for line in open(a.failed_log, errors='replace'):
            m = re.search(r"chunk failed: (\S+)\.chunk\d+\.(\d+)-(\d+) ::", line)
            if m:
                failed[m.group(1)].append((int(m.group(2)) - 1, int(m.group(3))))
    def unreachable(c, s, e):
        # holding chunks are 500 kb with 25 kb overlap; reachable if any holding chunk did not fail
        step, size = 475000, 500000
        k0, holders = max(0, (e - size) // step), []
        for k in range(k0, s // step + 1):
            cs, ce = k * step, k * step + size
            if cs <= s and e <= ce:
                holders.append((cs, ce))
        return bool(holders) and all(any(fs == cs for fs, _ in failed.get(c, [])) for cs, _ in holders)
    for run in a.runs:
        calls = calls_of(run)
        n_calls = sum(len(v) for v in calls.values())
        hit = defaultdict(lambda: [0, 0])
        on_planted = 0
        for c, s, e, _name, cls, u in truth:
            s, e, u = int(s), int(e), float(u)
            ok = any(rec((s, e), x) >= a.rec for x in calls.get(c, []))
            if failed and unreachable(c, s, e):
                hit['unreachable'][0] += ok
                hit['unreachable'][1] += 1
            else:
                hit['reachable'][0] += ok
                hit['reachable'][1] += 1
            touch = any(min(e, x[1]) > max(s, x[0]) for x in calls.get(c, []))
            on_planted += touch
            for key in ('all', cls, 'u<0.01' if u < 0.01 else 'u0.01-0.02' if u < 0.02 else 'u0.02-0.03' if u < 0.03 else 'u>=0.03'):
                hit[key][0] += ok
                hit[key][1] += 1
        parts = '  '.join(f"{k} {v[0]}/{v[1]}={v[0]/v[1]:.2f}" for k, v in hit.items())
        print(f"{os.path.basename(run)}: clean calls {n_calls}; planted recovered: {parts}; "
              f"planted elements touched by any call: {on_planted}")


if __name__ == '__main__':
    main()

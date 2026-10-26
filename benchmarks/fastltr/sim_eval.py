#!/usr/bin/env python3
"""Score LTRquest clean calls on the PrinTE burn-in genome against its known intact LTR-RTs.

Truth: rows of the PrinTE BED whose name carries '#LTR' and not '_FRAG' (intact insertions).
A call is a true positive when its span and a truth span overlap reciprocally by >= --rec
(default 0.9, as in the README benchmark). Each truth element is matched at most once.
usage: sim_eval.py truth.bed RUN_DIR [RUN_DIR ...]
"""
import argparse
import glob
import os
from collections import defaultdict


def load_truth(path):
    t = defaultdict(list)
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if '#LTR' in f[3] and '_FRAG' not in f[3]:
            t[f[0]].append((int(f[1]), int(f[2])))
    return t


def load_calls(run):
    out = []
    for p in glob.glob(os.path.join(run, '*_depth*_clean_ltr.tsv')):
        hdr = None
        for line in open(p):
            if line.startswith('#'):
                hdr = line[1:].rstrip('\n').split('\t')
                continue
            f = dict(zip(hdr, line.rstrip('\n').split('\t')))
            chrom, span = f['seq_id'].split('#')[0].rsplit(':', 1)
            s, e = (int(x) for x in span.split('-'))
            tsd = f.get('tsd', 'NA') not in ('NA', '.', '')
            tg = f.get('motif', '').lower() == 'tg...ca'
            out.append((chrom, s - 1, e, f.get('confidence') or ('high' if tsd and tg else 'medium' if tsd or tg else 'low'),
                        f['seq_id']))
    return out


def rec(a, b):
    o = max(0, min(a[1], b[1]) - max(a[0], b[0]))
    return min(o / (a[1] - a[0]), o / (b[1] - b[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('truth')
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--rec', type=float, default=0.9)
    a = ap.parse_args()
    truth = load_truth(a.truth)
    nT = sum(len(v) for v in truth.values())
    tp_sets = {}
    for run in a.runs:
        calls = load_calls(run)
        used, tp, tp_names = set(), 0, set()
        byconf = defaultdict(lambda: [0, 0])
        for c, s, e, conf, name in calls:
            hit = None
            for i, t in enumerate(truth.get(c, [])):
                if (c, i) not in used and rec((s, e), t) >= a.rec:
                    hit = (c, i)
                    break
            byconf[conf][1] += 1
            if hit:
                used.add(hit)
                tp += 1
                byconf[conf][0] += 1
                tp_names.add(name)
        fp, fn = len(calls) - tp, nT - tp
        p, r = tp / max(1, len(calls)), tp / nT
        f1 = 2 * p * r / (p + r) if p + r else 0
        tp_sets[run] = used
        conf = '; '.join(f"{k} {v[0]}/{v[1]}" for k, v in sorted(byconf.items()))
        print(f"{os.path.basename(run)}: calls={len(calls)} TP={tp} FP={fp} FN={fn} "
              f"precision={p:.3f} recall={r:.3f} F1={f1:.3f}  (TP/calls by confidence: {conf})")
    if len(a.runs) == 2:
        x, y = (tp_sets[r] for r in a.runs)
        print(f"truth found by {os.path.basename(a.runs[0])} only: {len(x - y)}; "
              f"by {os.path.basename(a.runs[1])} only: {len(y - x)}; both: {len(x & y)}")


if __name__ == '__main__':
    main()

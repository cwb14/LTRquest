#!/usr/bin/env python3
"""Score LTRquest benchmark runs of one genome against its uncapped reference.

For every config the clean element set (all depths) is matched to the reference
config's: two elements match when each LTR overlaps its counterpart by >= 50% of
the shorter. Elements are stratified by what --ltr-timeout cost their location in
the capped run named by --capped-from (round 1 chunks, 500 kb, 25 kb overlap):
  lost2  every round-1 chunk holding the element was cut short for both callers
  lost1  ... for one of the two callers
  open   at least one holding chunk finished for both
and by tandem context (FasTAN arrays >= 10 kb; satellite = unit 10-999 bp,
long = unit >= 1 kb; telomeric = within --tel bp of a chromosome end).

usage: analyze.py --genome G --runs DIR --configs U T TR TA UA --ref U --capped-from T
                  --fai G.fa.fai --arrays G.arrays10k.bed [--support support.tsv] --out PREFIX
"""
import argparse
import bisect
import glob
import os
import re
from collections import Counter, defaultdict

CHUNK, OVERLAP = 500_000, 25_000


def load_elements(run_dir, genome):
    rows = []
    for path in sorted(glob.glob(os.path.join(run_dir, f"{genome}_LTRs_depth*_clean_ltr.tsv"))):
        depth = int(re.search(r"_depth(\d+)_clean", path).group(1))
        hdr = None
        for line in open(path):
            if line.startswith('#'):
                hdr = line[1:].rstrip('\n').split('\t')
                continue
            if not line.strip():
                continue
            f = dict(zip(hdr, line.rstrip('\n').split('\t')))
            chrom, span = f['seq_id'].split('#')[0].rsplit(':', 1)
            s1, e1 = (int(x) for x in span.split('-'))
            g5 = (s1 + int(f['ltr5_start']) - 2, s1 + int(f['ltr5_end']) - 1)
            g3 = (s1 + int(f['ltr3_start']) - 2, s1 + int(f['ltr3_end']) - 1)
            tsd = f.get('tsd', 'NA') not in ('NA', '.', '')
            tgca = f.get('motif', '').lower() == 'tg...ca'
            conf = f.get('confidence') or ('high' if tsd and tgca else 'medium' if tsd or tgca else 'low')
            try:
                ident = float(f['identity'])
            except (KeyError, ValueError):
                ident = float('nan')
            rows.append({"name": f['seq_id'], "chrom": chrom, "start": s1 - 1, "end": e1, "g5": g5, "g3": g3,
                             "depth": depth, "conf": conf, "dom": f.get('domains', '.') not in ('.', '', 'NA'),
                             "ident": ident, "cls": f['seq_id'].split('#')[1] if '#' in f['seq_id'] else '?'})
    return rows


def ov(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def same(x, y):
    return (ov(x['g5'], y['g5']) >= 0.5 * min(x['g5'][1] - x['g5'][0], y['g5'][1] - y['g5'][0]) and
            ov(x['g3'], y['g3']) >= 0.5 * min(x['g3'][1] - x['g3'][0], y['g3'][1] - y['g3'][0]))


def matched(A, B):
    """Indices of A with a match in B, and of B with a match in A."""
    idx = defaultdict(list)
    for j, y in enumerate(B):
        idx[y['chrom']].append((y['g5'][0], j))
    for v in idx.values():
        v.sort()
    keys = {c: [s for s, _ in v] for c, v in idx.items()}
    mA, mB = set(), set()
    for i, x in enumerate(A):
        v = idx.get(x['chrom'], [])
        k = bisect.bisect_left(keys.get(x['chrom'], []), x['g5'][0] - 20000)
        while k < len(v) and v[k][0] < x['g5'][1]:
            if same(x, B[v[k][1]]):
                mA.add(i)
                mB.add(v[k][1])
            k += 1
    return mA, mB


def chunks_of(lengths):
    out = defaultdict(list)
    for c, L in lengths.items():
        s = 0
        while True:
            e = min(L, s + CHUNK)
            out[c].append((s, e))
            if e >= L:
                break
            s += CHUNK - OVERLAP
    return out


def load_capped(path):
    capped = defaultdict(set)   # (chrom, start) -> tools, round 1 only
    for line in open(path):
        c, s, e, tool, rnd = line.split()[:5]
        if int(rnd) == 1:
            capped[(c, int(s))].add(tool)
    return capped


def stratum(el, chunks, capped):
    holding = [(s, e) for s, e in chunks[el['chrom']] if s <= el['start'] and el['end'] <= e]
    if not holding:
        return 'span'          # wider than any round-1 chunk overlap
    lost = [len(capped.get((el['chrom'], s), ())) for s, _ in holding]
    if min(lost) >= 2:
        return 'lost2'
    if min(lost) == 1:
        return 'lost1'
    return 'open'


def load_bed(path):
    d = defaultdict(list)
    for line in open(path):
        f = line.split()
        d[(f[0], f[3])].append((int(f[1]), int(f[2])))
    for v in d.values():
        v.sort()
    return d


def hits(intervals, c, s, e):
    v = intervals.get(c, [])
    k = bisect.bisect_right(v, (s, float('inf'))) - 1
    k = max(k - 50, 0)
    for a, b in v[k:]:
        if a >= e:
            break
        if b > s:
            return True
    return False


def context(el, arrays, lengths, tel):
    c, s, e = el['chrom'], el['start'], el['end']
    if hits(arrays.get('satellite', {}), c, s, e):
        return 'satellite'
    if hits(arrays.get('long', {}), c, s, e):
        return 'long'
    if s < tel or e > lengths[c] - tel:
        return 'telomeric'
    return 'other'


def describe(rows, support):
    n = len(rows)
    if not n:
        return {"n": 0}
    c = Counter(r['conf'] for r in rows)
    idents = sorted(r['ident'] for r in rows if r['ident'] == r['ident'])
    out = {"n": n, "high": c['high'] / n, "medium": c['medium'] / n, "low": c['low'] / n,
               "domains": sum(r['dom'] for r in rows) / n,
               "ident": idents[len(idents) // 2] if idents else float('nan')}
    if support is not None:
        known = [support[r['name']] for r in rows if r['name'] in support]
        out['curated'] = sum(known) / len(known) if known else float('nan')
        out['curated_n'] = len(known)
    return out


def fmt(d):
    if not d.get('n'):
        return "n=0"
    s = (f"n={d['n']} high={d['high']:.2f} med={d['medium']:.2f} low={d['low']:.2f} "
         f"dom={d['domains']:.2f} ident={d['ident']:.3f}")
    if 'curated' in d:
        s += f" curated={d['curated']:.2f}"
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--genome', required=True)
    ap.add_argument('--runs', required=True)
    ap.add_argument('--configs', nargs='+', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--capped-from', required=True)
    ap.add_argument('--capped-bed', required=True)
    ap.add_argument('--fai', required=True)
    ap.add_argument('--arrays', required=True)
    ap.add_argument('--tel', type=int, default=100_000)
    ap.add_argument('--support', default=None, help='TSV name<TAB>0/1 curated-library support')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    lengths = {l.split()[0]: int(l.split()[1]) for l in open(a.fai)}
    chunks = chunks_of(lengths)
    capped = load_capped(a.capped_bed)
    arrays = defaultdict(dict)
    for (c, k), v in load_bed(a.arrays).items():
        arrays[k][c] = v
    support = None
    if a.support:
        support = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(a.support)}

    els = {cfg: load_elements(os.path.join(a.runs, f"{a.genome}_{cfg}"), a.genome) for cfg in a.configs}
    ref = els[a.ref]
    for x in ref:
        x['stratum'] = stratum(x, chunks, capped)
        x['context'] = context(x, arrays, lengths, a.tel)
    out = []
    out.append(f"# {a.genome}: reference {a.ref} ({len(ref)} clean elements); strata from {a.capped_from} round-1 caps")
    sc = Counter(x['stratum'] for x in ref)
    out.append("reference strata: " + ", ".join(f"{k} {sc[k]}" for k in ('open', 'lost1', 'lost2', 'span')))
    cc = Counter(x['context'] for x in ref)
    out.append("reference contexts: " + ", ".join(f"{k} {cc[k]}" for k in ('satellite', 'long', 'telomeric', 'other')))
    out.append(f"reference all: {fmt(describe(ref, support))}")
    tsv = ["genome\tconfig\tgroup\tstratum\tn\tn_ref\trecall\thigh\tmedium\tlow\tdomains\tident\tcurated"]
    for cfg in a.configs:
        if cfg == a.ref:
            continue
        B = els[cfg]
        for y in B:
            y['stratum'] = stratum(y, chunks, capped)
            y['context'] = context(y, arrays, lengths, a.tel)
        mA, mB = matched(ref, B)
        out.append(f"\n## {cfg}: {len(B)} clean elements ({len(B) - len(ref):+d} vs {a.ref})")
        for st in ('all', 'open', 'lost1', 'lost2'):
            R = [x for x in ref if st == 'all' or x['stratum'] == st]
            got = [x for i, x in enumerate(ref) if i in mA and (st == 'all' or x['stratum'] == st)]
            rec = len(got) / len(R) if R else float('nan')
            byc = []
            for conf in ('high', 'medium', 'low'):
                Rc = [x for x in R if x['conf'] == conf]
                gc = [x for x in got if x['conf'] == conf]
                byc.append(f"{conf} {len(gc)}/{len(Rc)}")
            sup_txt = ''
            if support is not None:
                Rs = [x for x in R if support.get(x['name'])]
                gs = [x for x in got if support.get(x['name'])]
                sup_txt = f"; curated-supported {len(gs)}/{len(Rs)} = {len(gs)/len(Rs) if Rs else float('nan'):.3f}"
            out.append(f"  recall of {a.ref} [{st:5s}]: {len(got)}/{len(R)} = {rec:.3f}  ({'; '.join(byc)}){sup_txt}")
            tsv.append(f"{a.genome}\t{cfg}\trecall\t{st}\t{len(got)}\t{len(R)}\t{rec:.4f}\t\t\t\t\t\t")
        added = [y for j, y in enumerate(B) if j not in mB]
        lost = [x for i, x in enumerate(ref) if i not in mA]
        for label, rows in (("added (not in ref)", added), ("lost (ref only)", lost)):
            out.append(f"  {label}: {fmt(describe(rows, support))}")
            for st in ('open', 'lost1', 'lost2'):
                sub = [r for r in rows if r['stratum'] == st]
                if sub:
                    out.append(f"      [{st:5s}] {fmt(describe(sub, support))}")
            d = describe(rows, support)
            if d.get('n'):
                tsv.append(f"{a.genome}\t{cfg}\t{label.split()[0]}\tall\t{d['n']}\t\t\t{d['high']:.3f}\t{d['medium']:.3f}"
                           f"\t{d['low']:.3f}\t{d['domains']:.3f}\t{d['ident']:.4f}\t{d.get('curated', float('nan')):.3f}")
        ctx = Counter(y['context'] for y in B)
        out.append("  contexts: " + ", ".join(f"{k} {ctx[k]} ({ctx[k] - cc[k]:+d})"
                                              for k in ('satellite', 'long', 'telomeric', 'other')))
    open(a.out + '.txt', 'w').write('\n'.join(out) + '\n')
    open(a.out + '.tsv', 'w').write('\n'.join(tsv) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()

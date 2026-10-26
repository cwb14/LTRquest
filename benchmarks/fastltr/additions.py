#!/usr/bin/env python3
"""Elements a test run adds over a reference run, by array context, with support and confidence.

Contexts as in analyze.py: long-unit arrays (rDNA), satellite arrays, telomeric ends, other.
usage: additions.py --genome G --ref REF_DIR --test TEST_DIR --fai G.fa.fai
                    --arrays G.arrays10k.bed --support support.tsv > additions.txt
"""
import argparse
from collections import Counter

import analyze as A


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--genome', required=True)
    ap.add_argument('--ref', required=True, help='run directory of the reference')
    ap.add_argument('--test', required=True, help='run directory of the run tested')
    ap.add_argument('--fai', required=True)
    ap.add_argument('--arrays', required=True)
    ap.add_argument('--support', required=True, help='TSV name<TAB>0/1')
    ap.add_argument('--tel', type=int, default=100_000)
    a = ap.parse_args()
    lengths = {line.split()[0]: int(line.split()[1]) for line in open(a.fai)}
    arrays = {}
    for (c, k), v in A.load_bed(a.arrays).items():
        arrays.setdefault(k, {})[c] = v
    sup = {line.split('\t')[0]: int(line.split('\t')[1]) for line in open(a.support)}
    ref = A.load_elements(a.ref, a.genome)
    test = A.load_elements(a.test, a.genome)
    _, matched = A.matched(ref, test)
    add = [y for j, y in enumerate(test) if j not in matched]
    print(f"additions: {len(add)} of {len(test)} test elements")
    for ctx in ('long', 'satellite', 'telomeric', 'other'):
        sub = [y for y in add if A.context(y, arrays, lengths, a.tel) == ctx]
        if sub:
            s = sum(sup.get(y['name'], 0) for y in sub)
            print(f"  [{ctx}] {len(sub)}, supported {s} ({s / len(sub):.2f}); "
                  f"confidence {dict(Counter(y['conf'] for y in sub))}")


if __name__ == '__main__':
    main()

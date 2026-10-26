#!/usr/bin/env python3
"""Plant known LTR-RTs into the real tandem arrays of a T2T assembly (exact truth in hard regions).

Each listed region is cut out of the genome as its own record. Insertion sites are drawn
inside its tandem arrays (FasTAN BED, >= 10 kb), >= --spacing bp apart and >= 20 kb from the
region ends. At each site one library element (random, either orientation) is inserted with a
5-bp target-site duplication; every base of the inserted copy is substituted with probability
u ~ U(0, --max-div), so its two LTRs differ by ~2u. Writes the spiked FASTA and a truth BED
(record, start, end, element, array class, u) in spiked coordinates (0-based, end exclusive,
TSD excluded).
usage: spike_in.py --genome G.fa --arrays G.arrays10k.bed --lib elements.fa.gz
                   --regions chr:start-end [...] --out PREFIX [--spacing 30000] [--seed 1]
"""
import argparse
import bisect
import gzip
import random


def fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    name, buf = None, []
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    yield name, ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, ''.join(buf)


def revcomp(s):
    return s.translate(str.maketrans('ACGTacgtNn', 'TGCAtgcaNn'))[::-1]


def mutate(s, u, rng):
    out = list(s)
    for i, b in enumerate(out):
        if b in 'ACGT' and rng.random() < u:
            out[i] = rng.choice('ACGT'.replace(b, ''))
    return ''.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--genome', required=True)
    ap.add_argument('--arrays', required=True)
    ap.add_argument('--lib', required=True)
    ap.add_argument('--regions', nargs='+', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--spacing', type=int, default=30000)
    ap.add_argument('--max-div', type=float, default=0.04)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    lib = [(n, s.upper()) for n, s in fasta(a.lib)]
    want = {}
    for r in a.regions:
        c, span = r.rsplit(':', 1)
        s, e = (int(x) for x in span.split('-'))
        want.setdefault(c, []).append((s, e))
    arrays = {}
    for line in open(a.arrays):
        f = line.split()
        arrays.setdefault(f[0], []).append((int(f[1]), int(f[2]), f[3]))
    fa_out = open(a.out + '.fa', 'w')
    bed_out = open(a.out + '.truth.bed', 'w')
    n_total = 0
    for chrom, seq in fasta(a.genome):
        for s, e in want.get(chrom, []):
            region = seq[s:e].upper()
            rec = f"{chrom}_{s}_{e}"
            # candidate sites: inside arrays, region-relative
            sites = []
            for x, y, cls in arrays.get(chrom, []):
                x, y = max(x, s + 20000), min(y, e - 20000)
                for p in range(x + 2000, y - 2000, 1000):
                    sites.append((p - s, cls))
            rng.shuffle(sites)
            chosen = []
            for p, cls in sites:
                k = bisect.bisect_left([q for q, _ in chosen], p)
                if all(abs(p - q) >= a.spacing for q, _ in chosen[max(0, k - 1):k + 1]):
                    bisect.insort(chosen, (p, cls))
            pieces, prev, shift = [], 0, 0
            for p, cls in chosen:
                name, el = rng.choice(lib)
                if rng.random() < 0.5:
                    el = revcomp(el)
                u = rng.uniform(0, a.max_div)
                el = mutate(el, u, rng)
                tsd = region[p - 5:p]
                pieces.append(region[prev:p])         # up to and including the TSD's first copy
                start = p + shift
                pieces.append(el + tsd)               # element, then the TSD's second copy
                bed_out.write(f"{rec}\t{start}\t{start + len(el)}\t{name}\t{cls}\t{u:.4f}\n")
                shift += len(el) + 5
                prev = p
                n_total += 1
            pieces.append(region[prev:])
            out = ''.join(pieces)
            fa_out.write(f">{rec}\n")
            for i in range(0, len(out), 80):
                fa_out.write(out[i:i + 80] + '\n')
    fa_out.close()
    bed_out.close()
    print(f"planted {n_total} elements")


if __name__ == '__main__':
    main()

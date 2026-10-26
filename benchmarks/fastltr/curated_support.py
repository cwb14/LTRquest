#!/usr/bin/env python3
"""Does a curated TE library support each LTRquest element? (80-80-80 style)

Every element of the given clean FASTAs is searched (blastn) against the curated
library's LTR entries; an element is supported when HSPs of >= 80% identity and
>= 80 bp cover >= 80% of its A/C/G/T bases (masked nested inners, stored as runs
of IUPAC letters, are not counted). Writes name<TAB>0/1<TAB>covered_fraction.

usage: curated_support.py --lib lib.fa [--lib-filter '#LTR'] --threads 16 --out support.tsv clean1.fa [clean2.fa ...]
"""
import argparse
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict


def read_fasta(path):
    name, buf = None, []
    opener = open
    if path.endswith('.gz'):
        import gzip
        opener = gzip.open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, ''.join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--lib', required=True)
    ap.add_argument('--lib-filter', default='#LTR', help="keep library entries whose header contains this")
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--out', required=True)
    ap.add_argument('fastas', nargs='+')
    a = ap.parse_args()

    seqs = {}
    for fa in a.fastas:
        for n, s in read_fasta(fa):
            seqs.setdefault(n, s)
    tmp = tempfile.mkdtemp(prefix='cursup.', dir=os.path.dirname(os.path.abspath(a.out)))
    q = os.path.join(tmp, 'q.fa')
    lib = os.path.join(tmp, 'lib.fa')
    acgt = {}
    with open(q, 'w') as fh:
        for i, (n, s) in enumerate(seqs.items()):
            fh.write(f'>e{i}\n{s}\n')
            acgt[f'e{i}'] = (n, sum(1 for ch in s if ch in 'ACGTacgt'), s)
    with open(lib, 'w') as fh:
        for n, s in read_fasta(a.lib):
            if a.lib_filter in n:
                fh.write(f'>{n}\n{s}\n')
    subprocess.run(['makeblastdb', '-in', lib, '-dbtype', 'nucl', '-out', lib], check=True,
                   stdout=subprocess.DEVNULL)
    res = subprocess.run(['blastn', '-task', 'blastn', '-query', q, '-db', lib, '-evalue', '1e-10',
                          '-num_threads', str(a.threads), '-max_target_seqs', '200', '-max_hsps', '20',
                          '-outfmt', '6 qseqid qstart qend pident length'],
                         check=True, capture_output=True, text=True).stdout
    shutil.rmtree(tmp)
    cov = defaultdict(list)
    for line in res.splitlines():
        qid, qs, qe, pid, ln = line.split('\t')
        if float(pid) >= 80 and int(ln) >= 80:
            qs, qe = sorted((int(qs), int(qe)))
            cov[qid].append((qs - 1, qe))
    with open(a.out, 'w') as out:
        for qid, (name, n_acgt, s) in acgt.items():
            mask = bytearray(len(s))
            for qs, qe in cov.get(qid, []):
                mask[qs:qe] = b'\x01' * (qe - qs)
            covered = sum(1 for i, ch in enumerate(s) if mask[i] and ch in 'ACGTacgt')
            frac = covered / n_acgt if n_acgt else 0.0
            out.write(f'{name}\t{int(frac >= 0.8)}\t{frac:.3f}\n')


if __name__ == '__main__':
    main()

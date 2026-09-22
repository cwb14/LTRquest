"""A synthetic LTR-RT family with truncated calls, written as LTRquest writes a run.

One contig, `chrS`. Family `syn_fam00001` (LTR 400 bp, TG..CA, internal 2 kb) has
13 intact copies whose calls are right and six whose calls stop short at an
obstacle between their two LTRs -- the way LTRharvest / LTR_FINDER stop:

  del_left         60 bp deleted 50 bp into the left LTR; the call starts after it
  patch_right      30 bp at 30% divergence, 40-70 bp from the right end; the call ends before it
  ins_left         1.5 kb inserted 40 bp into the left LTR; the call starts after it
  del_left_minus   del_left on a minus-strand copy: the obstacle sits at the forward right end
  nested_del_left  del_left inside the internal region of a one-copy host (depth 1)
  conflict_left    del_left with an annotated neighbour overlapping the recoverable bases

Every copy has a 5 bp TSD at its true ends except `decayed`, an intact copy whose
TSD was broken. Tables carry LTRquest's 33 columns; FASTAs store minus-oriented
records reverse-complemented and mask the host's nested inner with `N` (depth 0).
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import List

from ltrquest.detect import revcomp as revcomp_record
from ltrquest.kmer2ltr import COLUMNS

CHROM = "chrS"
FAMILY = "syn_fam00001"
LTR_LEN, INT_LEN, SPACER = 400, 2000, 1500
TAIL = ["strand", "family", "domains", "nest_status"]
_COMP = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass
class Element:
    kind: str
    family: str
    strand: str
    true_start: int
    true_end: int
    start: int
    end: int
    l1: int
    r0: int
    depth: int = 0
    tsd: str = "."
    k2p: float = 0.01
    nest_status: str = "."

    @property
    def key(self) -> str:
        return f"{CHROM}:{self.start}-{self.end}"

    @property
    def name(self) -> str:
        return f"{self.key}#LTR/Gypsy/Synth"


@dataclass
class Fixture:
    genome: Path
    indir: Path
    prefix: str
    elements: List[Element]
    ltr: str

    def kind(self, kind: str) -> Element:
        (hit,) = [e for e in self.elements if e.kind == kind]
        return hit

    def family_members(self) -> List[Element]:
        return [e for e in self.elements if e.family == FAMILY]


class _Contig:
    def __init__(self, rng: random.Random):
        self.rng = rng
        self.parts: List[str] = []
        self.pos = 0

    def rand(self, n: int) -> str:
        return "".join(self.rng.choice("ACGT") for _ in range(n))

    def add(self, seq: str) -> int:
        """Append `seq`; return its 1-based start."""
        self.parts.append(seq)
        start = self.pos + 1
        self.pos += len(seq)
        return start

    def tsd(self) -> str:
        while True:
            t = self.rand(5)
            if len(set(t)) >= 2:
                return t


def _mutate(rng: random.Random, seq: str, rate: float,
            keep_ends: int = 2) -> str:
    out = list(seq)
    for i in range(keep_ends, len(out) - keep_ends):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def _patch(rng: random.Random, seq: str, n: int) -> str:
    out = list(seq)
    for i in rng.sample(range(len(out)), n):
        out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def build(root: Path, prefix: str = "syn_LTRs",
          seed: int = 7) -> Fixture:
    rng = random.Random(seed)
    g = _Contig(rng)
    ltr0 = "TG" + g.rand(LTR_LEN - 4) + "CA"
    internal0 = g.rand(INT_LEN)
    elements: List[Element] = []

    def copy_pair():
        base = _mutate(rng, ltr0, 0.02)
        return (_mutate(rng, base, 0.003),
                _mutate(rng, base, 0.003),
                _mutate(rng, internal0, 0.02, keep_ends=0))

    def place(kind, left, internal, right, strand, family=FAMILY,
              broken_tsd=False, k2p=0.01):
        """Insert one copy between TSDs. `left`/`right` are biological 5'/3' LTRs."""
        t = g.tsd()
        g.add(t)
        body = (left + internal + right if strand == "+"
                else rc(left + internal + right))
        start = g.add(body)
        end = start + len(body) - 1
        g.add(t if not broken_tsd
              else t[:2] + ("A" if t[2] != "A" else "C") + t[3:])
        forward_left = len(left) if strand == "+" else len(right)
        forward_right = len(right) if strand == "+" else len(left)
        e = Element(kind, family, strand, start, end, start, end,
                    start + forward_left - 1, end - forward_right + 1,
                    tsd="." if broken_tsd else t, k2p=k2p)
        elements.append(e)
        return e

    for i in range(12):
        g.add(g.rand(SPACER))
        a, b, internal = copy_pair()
        place("normal", a, internal, b, "+" if i % 3 else "-",
              k2p=0.004 + 0.001 * i)
    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    place("decayed", a, internal, b, "+", broken_tsd=True, k2p=0.02)

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("del_left", a[:50] + a[110:], internal, b, "+", k2p=0.03)
    e.start, e.r0, e.tsd = e.true_start + 50, e.r0 + 110, "."

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("patch_right", a, internal,
              b[:330] + _patch(rng, b[330:360], 9) + b[360:], "+",
              k2p=0.03)
    e.end, e.l1, e.tsd = e.true_end - 70, e.l1 - 70, "."

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("ins_left", a[:40] + g.rand(1500) + a[40:], internal, b, "+",
              k2p=0.03)
    e.start, e.r0, e.tsd = e.true_start + 40 + 1500, e.r0 + 40, "."

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("del_left_minus", a[:50] + a[110:], internal, b, "-",
              k2p=0.03)
    e.end, e.l1, e.tsd = e.true_end - 50, e.l1 - 110, "."

    g.add(g.rand(SPACER))
    host_ltr = "TG" + g.rand(296) + "CA"
    host_tsd = g.tsd()
    g.add(host_tsd)
    host_start = g.add(host_ltr + g.rand(1500))
    a, b, internal = copy_pair()
    inner = place("nested_del_left", a[:50] + a[110:], internal, b, "+",
                  k2p=0.03)
    inner.start, inner.r0, inner.tsd = (inner.true_start + 50,
                                        inner.r0 + 110, ".")
    g.add(g.rand(1500) + _mutate(rng, host_ltr, 0.01))
    host_end = g.pos
    g.add(host_tsd)
    host = Element("host", "syn_fam00002", "+", host_start, host_end,
                   host_start, host_end, host_start + 299,
                   host_end - 299, depth=1, tsd=host_tsd, k2p=0.01)
    elements.append(host)
    host.nest_status = f"nest-outer:{inner.key}"
    inner.nest_status = f"nest-inner:{host.key}"

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("conflict_left", a[:50] + a[110:], internal, b, "+", k2p=0.03)
    e.start, e.r0, e.tsd = e.true_start + 50, e.r0 + 110, "."
    xs, xe = e.true_start - 405, e.true_start + 20
    elements.append(Element("neighbour", "syn_fam00003", "+", xs, xe, xs, xe,
                           xs + 99, xe - 99))
    g.add(g.rand(SPACER))

    contig = "".join(g.parts)
    root.mkdir(parents=True, exist_ok=True)
    genome = root / "genome.fa"
    with open(genome, "w") as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, len(contig), 60):
            fh.write(contig[i:i + 60] + "\n")

    header = "#" + "\t".join(COLUMNS + TAIL)
    for depth in (0, 1):
        rows = [e for e in elements if e.depth == depth]
        for variant in ("_clean_ltr", "_ltr"):
            base = root / f"{prefix}_depth{depth}{variant}"
            with open(f"{base}.tsv", "w") as fh:
                fh.write(header + "\n")
                for e in rows:
                    fh.write("\t".join(_row(e)) + "\n")
            with open(f"{base}.fa", "w") as fh:
                for e in rows:
                    seq = contig[e.start - 1:e.end]
                    if e.kind == "host":
                        a0, a1 = inner.start - e.start, inner.end - e.start + 1
                        seq = seq[:a0] + "N" * (a1 - a0) + seq[a1:]
                    if e.strand == "-":
                        seq = revcomp_record(seq)
                    fh.write(f">{e.name}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
    return Fixture(genome, root, prefix, elements, ltr0)


def _row(e: Element) -> List[str]:
    n = e.end - e.start + 1
    l5, l3 = e.l1 - e.start + 1, e.end - e.r0 + 1
    vals = {
        "seq_id": e.name,
        "seq_len": str(n),
        "status": "pass",
        "ltr5_start": "1",
        "ltr5_end": str(l5),
        "ltr3_start": str(e.r0 - e.start + 1),
        "ltr3_end": str(n),
        "ltr5_len": str(l5),
        "ltr3_len": str(l3),
        "flank5_len": "0",
        "flank3_len": "0",
        "aln_len": str(max(l5, l3)),
        "n_sites": str(min(l5, l3)),
        "n_ts": "1",
        "n_tv": "1",
        "n_gapcols": "0",
        "identity": "0.99",
        "p_dist": "0.01",
        "k2p": f"{e.k2p:g}",
        "k2p_se": "0.001",
        "bitscore": "500",
        "flank_margin_bits": "NA",
        "cigar": ".",
        "motif": "tg...ca",
        "k2p_time": str(round(e.k2p / 6e-8)),
        "orientation": "-" if e.strand == "-" else "+",
        "tsd": e.tsd,
        "tsd_offset": "NA" if e.tsd == "." else "0,0",
        "tsd_input": e.tsd,
        "strand": e.strand,
        "family": e.family,
        "domains": ".",
        "nest_status": e.nest_status,
    }
    return [vals[c] for c in COLUMNS + TAIL]

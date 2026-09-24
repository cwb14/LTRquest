"""Genome access and alignment primitives for ltrquest.reboundary.

`Member` is one called element as the clean tables describe it; `Genomes` gives
random access to every genome of a run by prefix, never reading one whole; `glocal`
aligns a model LTR end to end inside a window.

Coordinates are genome-forward, 1-based, inclusive. `Member.l1` and
`Member.r0` are the inner ends of the genomic-left and genomic-right LTR,
which is how the depth tables store them even for minus-strand elements.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import parasail
import pyfaidx

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2, -3, 5, 2
TERM = 30         # outer bases whose identity gates a moved end
_COMP = str.maketrans("ACGTN", "TGCAN")
_MATRIX = None


def matrix():
    """The parasail scoring matrix, built once per process (it does not
    pickle)."""
    global _MATRIX
    if _MATRIX is None:
        _MATRIX = parasail.matrix_create("ACGTN", MATCH, MISMATCH)
    return _MATRIX


def rc(seq: str) -> str:
    """Reverse complement of genomic sequence (ACGTN). Records use
    detect.revcomp."""
    return seq.translate(_COMP)[::-1]


@dataclass(frozen=True)
class Member:
    prefix: str
    name: str
    chrom: str
    start: int
    end: int
    l1: int
    r0: int
    strand: str
    orientation: str
    family: str
    depth: int
    k2p: Optional[float]
    tsd: str
    nest_status: str = "."
    tsd_offset: str = ""        # Kmer2LTR's d5,d3: "0,0" = the duplication sits at the called ends

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def uid(self) -> str:
        """Unique across genomes: the same coordinates can occur in two
        assemblies."""
        return f"{self.prefix}\t{self.key}"

    @property
    def len_left(self) -> int:
        return self.l1 - self.start + 1

    @property
    def len_right(self) -> int:
        return self.end - self.r0 + 1

    @property
    def stranded(self) -> bool:
        return self.strand in ("+", "-")

    @property
    def has_tsd(self) -> bool:
        return self.tsd not in ("", ".", "NA")


class Genomes:
    """Random access to each genome by prefix. Nothing is ever read whole."""

    def __init__(self, paths: Dict[str, str]):
        self.paths = dict(paths)
        self._open: Dict[str, pyfaidx.Fasta] = {}

    def __getstate__(self):
        return {"paths": self.paths}

    def __setstate__(self, state):
        self.paths = state["paths"]
        self._open = {}

    def index(self) -> None:
        """Build any missing or stale .fai now, in this one process.

        Call before a pool starts. Workers open genomes lazily and all at once,
        and pyfaidx locks per process: left to them, every worker rebuilds the
        same .fai, and some read it while another has it truncated.
        """
        for path in self.paths.values():
            pyfaidx.Fasta(path).close()

    def _fasta(self, prefix: str) -> pyfaidx.Fasta:
        fa = self._open.get(prefix)
        if fa is None:
            fa = pyfaidx.Fasta(self.paths[prefix], as_raw=True,
                               sequence_always_upper=True)
            self._open[prefix] = fa
        return fa

    def length(self, prefix: str, chrom: str) -> int:
        return len(self._fasta(prefix)[chrom])

    def fetch(self, prefix: str, chrom: str, lo: int,
              hi: int) -> Tuple[str, int]:
        """Bases lo..hi (1-based, inclusive), clipped to the sequence.
        Returns (seq, first)."""
        n = self.length(prefix, chrom)
        a, b = max(1, lo), min(n, hi)
        if a > b:
            return "", a
        return self._fasta(prefix)[chrom][a - 1:b], a


@dataclass(frozen=True)
class Glocal:
    start: int        # 0-based first window base the model covers
    end: int          # 0-based last window base
    id_first: float   # identity over the model's first `term` bases
    id_last: float
    id_all: float
    score: int
    head: Tuple[int, ...] = ()   # window index of model bases 0..edge-1
    tail: Tuple[int, ...] = ()   # same for the last bases, outermost first


def _identity(q: str, t: str) -> float:
    n = sum(1 for a in q if a != "-")
    m = sum(1 for a, b in zip(q, t) if a == b and a != "-")
    return m / n if n else 0.0


def glocal(model: str, window: str, term: int = TERM,
           edge: int = 8) -> Glocal:
    """Align `model` end to end inside `window`; only the window's ends are
    free."""
    r = parasail.sg_dx_trace_scan_sat(model, window, GAP_OPEN, GAP_EXTEND,
                                      matrix())
    q, t = r.traceback.query, r.traceback.ref
    i0 = 0
    while i0 < len(q) and q[i0] == "-":
        i0 += 1
    i1 = len(q) - 1
    while i1 >= 0 and q[i1] == "-":
        i1 -= 1
    start = i0 - t[:i0].count("-")
    end = start + sum(1 for c in t[i0:i1 + 1] if c != "-") - 1
    q, t = q[i0:i1 + 1], t[i0:i1 + 1]
    k = j = 0
    while j < len(q) and k < term:
        k += q[j] != "-"
        j += 1
    k, j2 = 0, len(q)
    while j2 > 0 and k < term:
        j2 -= 1
        k += q[j2] != "-"
    qmap: List[int] = []
    wi = start
    for a, b in zip(q, t):
        if a != "-":
            qmap.append(wi if b != "-" else -1)
        if b != "-":
            wi += 1
    return Glocal(start, end, _identity(q[:j], t[:j]),
                  _identity(q[j2:], t[j2:]), _identity(q, t), r.score,
                  tuple(qmap[:edge]), tuple(reversed(qmap[-edge:])))

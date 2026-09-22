"""Family LTR models for ltrquest.reboundary: what a family's LTR is, end to
end.

A model's termini come only from sequence agreement -- where a family's 5' and
3' LTR copies stop matching each other -- never from a target-site duplication
or a terminal motif, so the TSD stays an independent check on what reboundary
does.

Coordinates are genome-forward, 1-based, inclusive. `Member.l1` and
`Member.r0` are the inner ends of the genomic-left and genomic-right LTR,
which is how the depth tables store them even for minus-strand elements.
"""

from __future__ import annotations

import math
import os
import random
import subprocess
import tempfile
import zlib
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

import numpy as np
import parasail
import pyfaidx

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2, -3, 5, 2
TERM = 30         # outer bases whose identity gates a moved end
PAD = 200         # flank kept around each reference LTR in the model alignment
MIN_REFS = 10     # fewer reference copies than this and no model is built
KMER = 11
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


def canonical_kmers(seq: str, k: int = KMER) -> FrozenSet[str]:
    """Strand-independent k-mer set: each k-mer or its reverse complement,
    whichever is smaller."""
    out = set()
    for i in range(len(seq) - k + 1):
        w = seq[i:i + k]
        if "N" in w:
            continue
        r = rc(w)
        out.add(w if w <= r else r)
    return frozenset(out)


def jaccard(a: FrozenSet[str], b: FrozenSet[str]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


@dataclass(frozen=True)
class Model:
    model_id: str
    family: str
    seq: str            # biological 5'->3'
    modal_len: float    # the family's modal called LTR length
    n_refs: int
    kmers: FrozenSet[str] = field(default=frozenset(), repr=False,
                                   compare=False)

    def __post_init__(self):
        if not self.kmers:
            object.__setattr__(self, "kmers", canonical_kmers(self.seq))


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


def oriented_ltrs(m: Member, g: Genomes, pad: int = 0) -> Tuple[str, str]:
    """(5' LTR, 3' LTR), biological 5'->3', each padded by `pad` bases of
    context."""
    left, _ = g.fetch(m.prefix, m.chrom, m.start - pad, m.l1 + pad)
    right, _ = g.fetch(m.prefix, m.chrom, m.r0 - pad, m.end + pad)
    return (left, right) if m.strand != "-" else (rc(right), rc(left))


def modal_length(lengths: Sequence[int]) -> float:
    """Mode of LTR length: densest 5%-wide log bin (smoothed by its
    neighbours), its median."""
    bins = Counter(int(math.log(x) / math.log(1.05)) for x in lengths if x > 0)
    top = max(bins, key=lambda b: bins[b] + 0.5 * (bins.get(b - 1, 0) +
                                                     bins.get(b + 1, 0)))
    inbin = [x for x in lengths if x > 0 and abs(
        int(math.log(x) / math.log(1.05)) - top) <= 1]
    return float(np.median(inbin))


def select_references(members: Sequence[Member], strategy: str,
                      seed_key: str, exclude: FrozenSet[str] = frozenset(),
                      n_young: int = 40,
                      n_random: int = 40) -> Tuple[List[Member],
                                                    Optional[float]]:
    """Reference copies for a model, and the family's modal called LTR length.

    `modal`: both called LTRs within 15% of the modal length. `tsd`: copies
    whose called ends carry a Kmer2LTR TSD (independent evidence the call is
    right), falling back to `modal` below MIN_REFS. The youngest first, then a
    seeded random draw of the rest, so the choice does not depend on table
    order.
    """
    stranded = [m for m in members if m.stranded and m.key not in exclude]
    if len(stranded) < MIN_REFS:
        return [], None
    modal = modal_length([m.len_left for m in stranded] +
                         [m.len_right for m in stranded])
    pool: List[Member] = []
    if strategy == "tsd":
        pool = [m for m in stranded if m.has_tsd]
    if len(pool) < MIN_REFS:
        pool = [m for m in stranded
                if abs(m.len_left - modal) <= 0.15 * modal
                and abs(m.len_right - modal) <= 0.15 * modal]
    if len(pool) < MIN_REFS:
        return [], modal
    age = lambda m: (math.inf if m.k2p is None else m.k2p, m.key)  # noqa: E731
    young = sorted(pool, key=age)[:n_young]
    taken = {m.key for m in young}
    rest = sorted((m for m in pool if m.key not in taken), key=lambda m: m.key)
    rng = random.Random(zlib.crc32(seed_key.encode()))
    return young + rng.sample(rest, min(n_random, len(rest))), modal


def mafft_align(seqs: Sequence[str], mafft: str = "mafft") -> List[str]:
    fd, path = tempfile.mkstemp(suffix=".fa")
    try:
        with os.fdopen(fd, "w") as fo:
            for i, s in enumerate(seqs):
                fo.write(f">{i}\n{s}\n")
        out = subprocess.run([mafft, "--auto", "--thread", "1", "--quiet",
                             path], capture_output=True, text=True,
                            check=True).stdout
    finally:
        os.unlink(path)
    rows: Dict[int, List[str]] = {}
    cur = -1
    for line in out.splitlines():
        if line.startswith(">"):
            cur = int(line[1:])
            rows[cur] = []
        elif cur >= 0:
            rows[cur].append(line.strip().upper())
    return ["".join(rows[i]) for i in range(len(seqs))]


def consensus_from_alignment(aln: Sequence[str],
                             sides: Sequence[int]) -> Optional[str]:
    """The LTR in an alignment of padded 5'-LTR (side 5) and 3'-LTR (side 3)
    rows.

    A column is LTR when both sides agree on a base carried by >= 40% of each.
    Beyond the outer end the 5' rows are flank; beyond the inner end the 3'
    rows are; either way the sides disagree. The span runs from the first to
    the last sustained block (>= 15 of 21 columns) so dips inside the LTR are
    ignored.
    """
    A = np.frombuffer("".join(aln).encode(), dtype="S1").reshape(len(aln), -1)
    side = np.asarray(sides)
    bases = np.array([b"A", b"C", b"G", b"T"])
    c5 = np.stack([(A[side == 5] == b).sum(0) for b in bases])
    c3 = np.stack([(A[side == 3] == b).sum(0) for b in bases])
    f5 = c5.max(0) / max(1, int((side == 5).sum()))
    f3 = c3.max(0) / max(1, int((side == 3).sum()))
    agree = (c5.argmax(0) == c3.argmax(0)) & (f5 >= 0.4) & (f3 >= 0.4)
    smooth = np.convolve(agree.astype(float), np.ones(21) / 21, mode="same")
    hit = np.flatnonzero(smooth >= 15 / 21)
    if len(hit) == 0:
        return None
    lo, hi, n = int(hit[0]), int(hit[-1]), len(agree)
    if hi - lo < 50:
        return None
    while lo > 0 and agree[lo - 1]:
        lo -= 1
    while hi < n - 1 and agree[hi + 1]:
        hi += 1
    while lo < hi and not agree[lo]:
        lo += 1
    while hi > lo and not agree[hi]:
        hi -= 1
    tot = c5 + c3
    occ = tot.sum(0) / len(aln)
    return "".join(bases[tot[:, j].argmax()].decode()
                   for j in range(lo, hi + 1) if occ[j] >= 0.5)


def polish_ends(seq: str, refs: Sequence[Member], g: Genomes, k: int = 8,
                thr: float = 0.7, rounds: int = 3) -> Tuple[str, Tuple[int,
                                                             int]]:
    """Move the model's ends to where each reference's two LTRs stop agreeing.

    That is the definition of an LTR end. The model is placed in both LTRs of
    every reference and the base pairs at and just beyond each end are
    compared: inside the LTR they agree (1 - divergence); beyond it one side
    is flank and the other internal region. Extend while >= `thr` of the
    references agree, trim while fewer do. A model base aligned to a gap
    counts as disagreement, so bases no reference carries are trimmed. Returns
    (seq, (bases added at the start, bases added at the end)); negative =
    trimmed.
    """
    total_s = total_e = 0
    for _ in range(rounds):
        placed = []
        for m in refs:
            s5, s3 = oriented_ltrs(m, g, PAD)
            g5, g3 = glocal(seq, s5, edge=k), glocal(seq, s3, edge=k)
            if min(g5.id_all, g3.id_all) >= 0.6:
                placed.append((s5, s3, g5, g3))
        if len(placed) < MIN_REFS:
            break

        def vote(pairs) -> Tuple[float, str]:
            votes: Counter = Counter()
            hits = 0
            for pair in pairs:
                if pair is not None and pair[0] == pair[1] and pair[0] != "N":
                    hits += 1
                    votes[pair[0]] += 1
            return hits / len(placed), (votes.most_common(1)[0][0]
                                        if votes else "N")

        def inside(i: int, at_start: bool) -> Tuple[float, str]:
            """Agreement at model base i counted from that end (0 = the
            outermost)."""
            pairs = []
            for s5, s3, g5, g3 in placed:
                a = (g5.head if at_start else g5.tail)[i]
                b = (g3.head if at_start else g3.tail)[i]
                pairs.append((s5[a], s3[b]) if a >= 0 and b >= 0 else None)
            return vote(pairs)

        def outside(j: int, at_start: bool) -> Tuple[float, str]:
            """Agreement j bases beyond that end of the model (1 = adjacent).
            """
            pairs = []
            for s5, s3, g5, g3 in placed:
                i5 = g5.start - j if at_start else g5.end + j
                i3 = g3.start - j if at_start else g3.end + j
                ok = 0 <= i5 < len(s5) and 0 <= i3 < len(s3)
                pairs.append((s5[i5], s3[i3]) if ok else None)
            return vote(pairs)

        def trim_to(at_start: bool) -> int:
            """Model bases to drop from this end: up to the first 3 agreeing
            in a row.

            A lone agreeing base is not trusted: the aligner parks an
            overhanging base on any matching base nearby, in both LTRs at once.
            """
            for t in range(k - 2):
                if all(inside(t + d, at_start)[0] >= thr for d in range(3)):
                    return t
            return 0

        ds = 0
        while ds < k and outside(ds + 1, True)[0] >= thr:
            ds += 1
        if ds == 0:
            ds = -trim_to(True)
        de = 0
        while de < k and outside(de + 1, False)[0] >= thr:
            de += 1
        if de == 0:
            de = -trim_to(False)
        if ds == 0 and de == 0:
            break
        pre = "".join(outside(j, True)[1] for j in range(ds, 0, -1)
                      ) if ds > 0 else ""
        post = "".join(outside(j, False)[1] for j in range(1, de + 1)
                       ) if de > 0 else ""
        seq = pre + seq[(-ds if ds < 0 else 0):(len(seq) + de if de < 0
                                                  else len(seq))] + post
        total_s += ds
        total_e += de
    return seq, (total_s, total_e)


def consensus_model(family: str, refs: Sequence[Member], g: Genomes,
                    modal: Optional[float], mafft: str = "mafft",
                    model_id: Optional[str] = None) -> Optional[Model]:
    """One model from `refs`: MAFFT their padded 5'/3' LTRs, take the agreed
    span, polish."""
    if len(refs) < MIN_REFS or modal is None:
        return None
    seqs: List[str] = []
    sides: List[int] = []
    for m in refs:
        s5, s3 = oriented_ltrs(m, g, PAD)
        seqs += [s5, s3]
        sides += [5, 3]
    seq = consensus_from_alignment(mafft_align(seqs, mafft), sides)
    if seq is None:
        return None
    seq, _ = polish_ends(seq, refs, g)
    return Model(model_id or f"{family}:consensus", family, seq, modal,
                 len(refs))


def ratio_ok(model: Model, max_ratio: float) -> bool:
    """A model much longer than the family's usual called LTR is not an LTR
    model."""
    return len(model.seq) <= max_ratio * model.modal_len


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p).

    Summed in log space: a family can put thousands of candidates through this
    test, and `math.comb(n, i)` for n in the thousands is an integer too large
    to multiply by a float.
    """
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    log_p, log_q, log_n = math.log(p), math.log1p(-p), math.lgamma(n + 1)
    terms = [log_n - math.lgamma(i + 1) - math.lgamma(n - i + 1)
             + i * log_p + (n - i) * log_q for i in range(k, n + 1)]
    hi = max(terms)
    if hi == -math.inf:
        return 0.0
    total = sum(math.exp(t - hi) for t in terms)
    return min(1.0, math.exp(hi + math.log(total)))


def tsd_enrichment(hits: int, n: int, p0: float, min_n: int,
                   alpha: float = 0.01) -> str:
    """Family QC on TSDs at proposed ends: 'pass', 'fail', or 'untested'
    (n < min_n).

    `p0` is the displaced-flank null rate, floored at 0.5% so one lucky family
    cannot pass on a null of zero.
    """
    if n < min_n:
        return "untested"
    return "pass" if binom_sf(hits, n, max(p0, 0.005)) < alpha else "fail"


def subfamily_models(family: str, refs: Sequence[Member], g: Genomes, modal: Optional[float],
                     jaccard_min: float, mafft: str = "mafft") -> List[Model]:
    """One consensus per cluster of similar reference LTRs (>= MIN_REFS copies), plus the rest.

    Clusters grow greedily around centroids taken youngest-first (select_references
    returns the youngest first), on canonical k-mer Jaccard of each reference's 5'
    LTR. With no cluster big enough, this is the family consensus.
    """
    profiles = [canonical_kmers(oriented_ltrs(m, g)[0]) for m in refs]
    centroids: List[int] = []
    clusters: List[List[int]] = []
    for i, prof in enumerate(profiles):
        for c, members in zip(centroids, clusters):
            if jaccard(prof, profiles[c]) >= jaccard_min:
                members.append(i)
                break
        else:
            centroids.append(i)
            clusters.append([i])
    models: List[Model] = []
    big = [cl for cl in clusters if len(cl) >= MIN_REFS]
    for n, cl in enumerate(big, start=1):
        mo = consensus_model(family, [refs[i] for i in cl[:80]], g, modal, mafft,
                             f"{family}:sub{n}")
        if mo is not None:
            models.append(mo)
    rest = [refs[i] for cl in clusters if len(cl) < MIN_REFS for i in cl]
    if big and len(rest) >= MIN_REFS:
        mo = consensus_model(family, rest[:80], g, modal, mafft, f"{family}:rest")
        if mo is not None:
            models.append(mo)
    if not models:
        mo = consensus_model(family, list(refs)[:80], g, modal, mafft)
        if mo is not None:
            models.append(mo)
    return models

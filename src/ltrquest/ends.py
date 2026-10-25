"""Canonical LTR-RT ends: TG..CA termini flanked by a target-site duplication.

An LTR retrotransposon is integrated into a staggered cut, so the few target bases
beside it end up on both sides (the TSD, 4-6 bp; 5 in most LTR-RTs), and integrase
leaves TG at the element's 5' end and CA at its 3' end. TG..CA reads the same on the
reverse strand, so both tests are strand-free and work in genome-forward coordinates.

Kmer2LTR records a TSD it finds after moving a boundary 1 bp (`tsd_offset`) but never
moves the boundary, and it does not look at the termini. Its first hit is not always
the element's ends: X1..X5|TG..CA|X1..X5 called 1 bp wide on one side carries an
exact 4-bp TSD at the called ends, and a 5-mer matches by chance at a 1-bp shift of a
correct call about one time in four. Two hallmarks that agree are hard to meet by
chance (TG..CA plus a 5-bp TSD: about 1 in 260,000 per boundary pair), so a boundary is
moved only onto ends that have both.

Coordinates follow Kmer2LTR's `find_tsd`: the element is `ctx[b5:b3]`, and a shift
(d5, d3) is positive into the element at either end.
"""

from __future__ import annotations

from typing import NamedTuple, Optional, Tuple

from .table import is_missing

MAX_SHIFT = 2               # bp either end may move; Kmer2LTR's own search tries 1
TSD_LENGTHS = (4, 5, 6)     # LTR-RT TSDs; 5 is by far the most common
MODAL_TSD = 5
FLANK = MAX_SHIFT + max(TSD_LENGTHS)     # bases each end's window needs either side


class Ends(NamedTuple):
    d5: int         # bp the left end moves right (into the element when positive)
    d3: int         # bp the right end moves left (into the element when positive)
    tsd: str


def tsd_at(ctx: str, b5: int, b3: int, ks: Tuple[int, ...] = TSD_LENGTHS) -> Optional[str]:
    """The exact duplication flanking `ctx[b5:b3]`, uppercase, or None.

    Lengths are tried nearest `MODAL_TSD` first. Two lengths fit at the same
    boundary only in periodic sequence, so the order rarely matters. Refused, as
    Kmer2LTR refuses them: a k-mer holding anything but ACGT (an assembly gap, a
    masked base), and one of fewer than two distinct bases (a homopolymer matches
    its own reflection at a large share of genomic positions).
    """
    if b3 - b5 < 4:
        return None
    for k in sorted(ks, key=lambda k: (abs(k - MODAL_TSD), k)):
        if b5 - k < 0 or b3 + k > len(ctx):
            continue
        left = ctx[b5 - k:b5].upper()
        if left != ctx[b3:b3 + k].upper():
            continue
        if left.strip("ACGT") or len(set(left)) < 2:
            continue
        return left
    return None


def canonical_tsd(ctx: str, b5: int, b3: int) -> Optional[str]:
    """The TSD when `ctx[b5:b3]` starts TG, ends CA and has an exact TSD, else None."""
    if ctx[b5:b5 + 2].upper() != "TG" or ctx[b3 - 2:b3].upper() != "CA":
        return None
    return tsd_at(ctx, b5, b3)


def snap(ctx: str, b5: int, b3: int, max_shift: int = MAX_SHIFT) -> Optional[Ends]:
    """The nearest canonical ends within `max_shift` bp of each called end, or None.

    Ranked by the total bases moved, then by TSD length nearest `MODAL_TSD`, then
    shorter, then by which end moves (deterministic). `Ends(0, 0, tsd)` means the
    called ends are already canonical.
    """
    best = None
    for d5 in range(-max_shift, max_shift + 1):
        for d3 in range(-max_shift, max_shift + 1):
            tsd = canonical_tsd(ctx, b5 + d5, b3 - d3)
            if tsd is None:
                continue
            key = (abs(d5) + abs(d3), abs(len(tsd) - MODAL_TSD), len(tsd), abs(d5), d5, d3)
            if best is None or key < best[0]:
                best = (key, Ends(d5, d3, tsd))
    return best[1] if best else None


def end_windows(seq: str, start: int, end: int, flank: int = FLANK) -> Tuple[str, int, int]:
    """`(ctx, b5, b3)` for an element at `start..end` (1-based, inclusive) of `seq`:
    `flank` bases either side of each end, joined, so `snap` sees both ends without
    the element's interior. An end near the contig edge gets a shorter window, and
    a candidate needing bases beyond it is simply not found."""
    lo = max(0, start - 1 - flank)
    left = seq[lo:start - 1 + flank]
    right = seq[end - flank:end + flank]
    b5 = start - 1 - lo
    return left + right, b5, len(left) + flank


def tsd_between(fetch, left: int, right: int) -> Optional[str]:
    """The TSD of an element at `left..right` (1-based, inclusive) whose ends are
    canonical, else None. `fetch(lo, hi)` returns (the bases lo..hi clipped to the
    contig, the position of the first), as ltr_model.Genomes.fetch does for a contig."""
    a, a0 = fetch(left - FLANK, left + FLANK - 1)
    b, _ = fetch(right - FLANK + 1, right + FLANK)
    return canonical_tsd(a + b, left - a0, len(a) + FLANK)

def hallmarks(tsd: str, tsd_offset: str, motif: str) -> Tuple[bool, bool]:
    """(an exact TSD at the call's own ends, TG..CA termini), from its table columns.
    A TSD Kmer2LTR found only after moving a boundary (`tsd_offset` other than 0,0) is
    evidence for other ends, not these."""
    has_tsd = not is_missing(tsd) and (tsd_offset or "").replace(" ", "") == "0,0"
    return has_tsd, (motif or "").lower() == "tg...ca"


def confidence(tsd: str, tsd_offset: str, motif: str) -> str:
    """`high` with both hallmarks at the reported ends, `medium` with one, `low` with
    neither.

    A real insertion is born with both. Its two TSD copies and its termini then
    mutate like the rest of the element, so an old real element can lose them, but
    a young one keeps them: in a young element their absence marks a wrong call.
    Species whose elements lack TG..CA reach at most `medium`.
    """
    has_tsd, tgca = hallmarks(tsd, tsd_offset, motif)
    return ("low", "medium", "high")[has_tsd + tgca]

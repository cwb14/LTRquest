"""One rule for calls that overlap: nested, or in conflict.

Real insertions are either disjoint or nested: an element inserted into another lies
strictly inside it, with an LTR pair of its own. Any other overlap -- a partial one
(two calls crossing), one call inside another but sharing an LTR with it, the same
span twice -- means at most one of the two calls describes a real insertion. The
calls are then ranked on structural evidence and the stronger one is kept:

  1. the hallmarks of an insertion at the call's own ends: an exact TSD (target-site
     duplication) and TG..CA termini. Both > TSD alone > TG..CA alone > neither. A TSD
     ranks above TG..CA because the LTRs of one family share their termini, so a
     pairing of LTRs from two different insertions keeps TG..CA but not a TSD;
  2. the LTR pair's alignment score (Kmer2LTR's `bitscore`), which grows with both
     LTR length and identity;
  3. how many distinct protein domains the call holds.

Calls are kept greedily, strongest first, skipping any that conflicts with one
already kept; ties fall to genomic order. This replaces the earlier rule, which
compared only calls overlapping by >= 80% of the shorter and kept every call that
took part in a nesting, so crossing calls survived it.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple

from .ends import hallmarks
from .table import as_float, is_missing

BIN = 100_000       # interval-index bin; an element spans at most a few


class Call(NamedTuple):
    key: str                                       # chrom:start-end
    chrom: str
    start: int                                     # 1-based, inclusive
    end: int
    ltrs: Optional[Tuple[int, int, int, int]]      # genomic left LTR start, end, right LTR start, end
    rank: Tuple                                    # `evidence`; larger is stronger


def n_domains(domains: str) -> int:
    """Distinct protein domains in a `domains` value (`GENE|clade@s-e;...`)."""
    if is_missing(domains):
        return 0
    return len({tok.split("|", 1)[0] for tok in domains.split(";") if tok})


def evidence(tsd: str, tsd_offset: str, motif: str, bitscore: str, domains: int) -> Tuple:
    has_tsd, tgca = hallmarks(tsd, tsd_offset, motif)
    return (2 * has_tsd + tgca, as_float(bitscore, 0.0), domains)


def _ltrs_overlap(a: Call, b: Call) -> bool:
    """Whether any LTR of one overlaps any LTR of the other. Unknown LTRs share nothing."""
    if a.ltrs is None or b.ltrs is None:
        return False
    for s1, e1 in ((a.ltrs[0], a.ltrs[1]), (a.ltrs[2], a.ltrs[3])):
        for s2, e2 in ((b.ltrs[0], b.ltrs[1]), (b.ltrs[2], b.ltrs[3])):
            if s1 <= e2 and s2 <= e1:
                return True
    return False


def relation(a: Call, b: Call) -> Optional[str]:
    """None (disjoint), 'nested' (one strictly inside the other, LTRs distinct) or
    'conflict' (any other overlap)."""
    if a.chrom != b.chrom or a.end < b.start or b.end < a.start:
        return None
    inside = (a.start >= b.start and a.end <= b.end) or (b.start >= a.start and b.end <= a.end)
    if inside and (a.end - a.start) != (b.end - b.start) and not _ltrs_overlap(a, b):
        return "nested"
    return "conflict"


def resolve(calls: Iterable[Call]) -> Tuple[List[Call], List[Tuple[Call, Call]]]:
    """(kept, [(dropped, the stronger call it conflicts with)]).

    Kept calls come back in input order.
    """
    calls = list(calls)
    order = sorted(range(len(calls)), key=lambda i: tuple(-x for x in calls[i].rank)
                   + (calls[i].start, calls[i].end, calls[i].key, i))
    bins: Dict[Tuple[str, int], List[Call]] = defaultdict(list)
    kept = [False] * len(calls)
    dropped: List[Tuple[Call, Call]] = []
    for i in order:
        c = calls[i]
        span = range(c.start // BIN, c.end // BIN + 1)
        winner = next((k for b in span for k in bins[(c.chrom, b)]
                       if relation(c, k) == "conflict"), None)
        if winner is not None:
            dropped.append((c, winner))
            continue
        kept[i] = True
        for b in span:
            bins[(c.chrom, b)].append(c)
    return [c for c, k in zip(calls, kept) if k], dropped


def nest_pairs(calls: Sequence[Call]) -> List[Tuple[Call, Call]]:
    """(inner, outer) for every nested pair among `calls`."""
    by_chrom: Dict[str, List[Call]] = defaultdict(list)
    for c in calls:
        by_chrom[c.chrom].append(c)
    out = []
    for group in by_chrom.values():
        group.sort(key=lambda c: (c.start, -c.end))
        for i, a in enumerate(group):
            for b in group[i + 1:]:
                if b.start > a.end:
                    break
                if relation(a, b) == "nested":
                    out.append((b, a) if b.end - b.start < a.end - a.start else (a, b))
    return out


CONFLICTS_HEADER = ("#dropped\tkept\tkind\tdropped_hallmarks\tdropped_bitscore\t"
                    "dropped_domains\tkept_hallmarks\tkept_bitscore\tkept_domains\n")
_HALLMARKS = ("none", "TG..CA", "TSD", "TSD+TG..CA")


def kind(a: Call, b: Call) -> str:
    """`same_span`, `shared_ltr` (one inside the other, an LTR in common) or `crossing`."""
    if (a.start, a.end) == (b.start, b.end):
        return "same_span"
    if (a.start >= b.start and a.end <= b.end) or (b.start >= a.start and b.end <= a.end):
        return "shared_ltr"
    return "crossing"


def conflict_row(loser: Call, winner: Call) -> List[str]:
    """One CONFLICTS_HEADER row: what was dropped, what it lost to, and on what evidence."""
    def fields(c):
        tier, bits, doms = c.rank
        return [_HALLMARKS[tier], f"{bits:g}", str(doms)]
    return [loser.key, winner.key, kind(loser, winner), *fields(loser), *fields(winner)]

"""Placing trusted templates on one element: where do its LTRs really end?

Each template contributes one model per side (`templates.models_for`): the LTR
whose outer end -- the end its own TSD verifies -- is outermost there. A model is
aligned end to end in a window around the called LTR on its side (the core). When
the core's outer bases do not match -- a large insertion sits between the true end
and the rest of the LTR -- the model's outermost `anchor_len` bases are searched
for locally up to `max_indel` beyond the call (the anchor), unless a complete copy
of the LTR starts (or ends) there: that is a neighbouring solo LTR or tandem copy,
not this element's terminus.

Per side, the placements that pass the identity gates vote: the low median, if at
least `min_support` of them agree with it within `agree` bp. An outer end then
moves outward by >= `min_ext` bp, or inward by at most `max_trim` bp (a small
over-call); a larger inward vote is ignored.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import parasail

from .ltr_model import (
    GAP_EXTEND,
    GAP_OPEN,
    Genomes,
    Member,
    glocal,
    matrix,
)


@dataclass(frozen=True)
class Params:
    min_identity: float = 0.8       # outer-TERM identity a placement needs to vote
    min_whole_identity: float = 0.6
    min_ext: int = 1                # an outward move must reach this far to count
    max_trim: int = 10              # the largest inward move (0 = never trim)
    anchor: bool = True
    anchor_len: int = 30
    max_indel: int = 5000
    n_templates: int = 5            # templates placed per element
    min_support: int = 2            # agreeing placements a moved end needs
    agree: int = 2                  # bp within which placements agree


@dataclass(frozen=True)
class Side:
    """One model placed on one side of an element."""
    outer: int          # this side's outer end (left: first base; right: last base)
    inner: int          # this side's inner end (left LTR's l1 / right LTR's r0)
    id_outer: float     # identity of the model's outer TERM bases
    whole: float        # identity of the whole model
    src: str = "core"   # core | anchor


@dataclass(frozen=True)
class Hit:
    q_start: int
    q_end: int
    w_start: int
    w_end: int
    identity: float
    score: int


@dataclass(frozen=True)
class Proposal:
    member: Member
    orient: str             # frame of ext5/ext3: the element's strand, else implied
    left: int               # committed outer ends (the call's own where nothing moved)
    right: int
    left_src: str           # core | anchor | '.' (did not move)
    right_src: str
    id_left: float          # median outer identity of the placements that carried it
    id_right: float
    sup_left: int           # placements agreeing with the voted end
    sup_right: int
    raw_left: int           # median of every placement, gates or not (reporting)
    raw_right: int
    inner_left: int         # median inner end on each side, for the obstacle report
    inner_right: int
    reason_left: str        # why a side did not move: '.', gate_identity, no_support,
    reason_right: str       # gap_in_added
    templates: Tuple[str, ...] = ()

    @property
    def gate_ok(self) -> bool:
        return self.left != self.member.start or self.right != self.member.end

    @property
    def ext_left(self) -> int:      # bp added at the left (negative: trimmed)
        return self.member.start - self.left

    @property
    def ext_right(self) -> int:
        return self.right - self.member.end

    @property
    def ext5(self) -> int:
        return self.ext_left if self.orient == "+" else self.ext_right

    @property
    def ext3(self) -> int:
        return self.ext_right if self.orient == "+" else self.ext_left


def local(query: str, window: str) -> Optional[Hit]:
    """Best local alignment of `query` in `window` (0-based, inclusive coordinates)."""
    if not query or not window:
        return None
    r = parasail.sw_trace_scan_sat(query, window, GAP_OPEN, GAP_EXTEND, matrix())
    if r.score <= 0:
        return None
    q, t = r.traceback.query, r.traceback.ref
    q_len = sum(1 for c in q if c != "-")
    t_len = sum(1 for c in t if c != "-")
    matches = sum(1 for a, b in zip(q, t) if a == b and a != "-")
    return Hit(r.end_query - q_len + 1, r.end_query, r.end_ref - t_len + 1, r.end_ref,
               matches / max(1, q_len), r.score)


def vote(positions: Sequence[int], agree: int, need: int) -> Tuple[Optional[int], int]:
    """(low median, support): support counts positions within `agree` of it; the
    median is returned only when at least `need` of them agree."""
    if not positions:
        return None, 0
    s = sorted(positions)
    med = s[(len(s) - 1) // 2]
    sup = sum(abs(x - med) <= agree for x in s)
    return (med if sup >= need else None), sup


def _extra(m: Member, model: str) -> int:
    """How far past the call a window reaches: the model may be longer than the
    shorter called LTR by any amount."""
    return max(100, len(model) - min(m.len_left, m.len_right) + 100)


def place_left(m: Member, model: str, g: Genomes) -> Optional[Side]:
    """The model aligned end to end around the called left LTR, if it lands on it.

    The window is widened by the truncation, so it can reach a neighbouring copy of
    the same family -- a solo LTR, an uncalled fragment, a tandem array member --
    which, being full length, would outscore the truncated call. A placement that
    does not overlap the LTR it was placed on is not a placement of this element.
    """
    extra = _extra(m, model)
    w, w0 = g.fetch(m.prefix, m.chrom, m.start - extra, m.l1 + extra)
    if len(w) < len(model) // 2:
        return None
    a = glocal(model, w)
    left, l1 = w0 + a.start, w0 + a.end
    if left > m.l1 or l1 < m.start:
        return None
    return Side(left, l1, a.id_first, a.id_all)


def place_right(m: Member, model: str, g: Genomes) -> Optional[Side]:
    """`place_left` for the called right LTR."""
    extra = _extra(m, model)
    w, w0 = g.fetch(m.prefix, m.chrom, m.r0 - extra, m.end + extra)
    if len(w) < len(model) // 2:
        return None
    b = glocal(model, w)
    r0, right = w0 + b.start, w0 + b.end
    if r0 > m.end or right < m.r0:
        return None
    return Side(right, r0, b.id_last, b.id_all)


def anchor_left(m: Member, fwd: str, g: Genomes, p: Params) -> Optional[Tuple[int, float]]:
    q = fwd[:p.anchor_len]
    w, w0 = g.fetch(m.prefix, m.chrom, m.start - p.max_indel - len(q), m.start + 200)
    h = local(q, w)
    if h is None or h.q_start > 2 or (h.q_end - h.q_start + 1) < 0.9 * len(q):
        return None
    pos = w0 + h.w_start - h.q_start       # the model's terminal base, extrapolated
    return (pos, h.identity) if pos >= 1 else None


def anchor_right(m: Member, fwd: str, g: Genomes, p: Params) -> Optional[Tuple[int, float]]:
    q = fwd[-p.anchor_len:]
    w, w0 = g.fetch(m.prefix, m.chrom, m.end - 200, m.end + p.max_indel + len(q))
    h = local(q, w)
    if h is None or h.q_end < len(q) - 3 or (h.q_end - h.q_start + 1) < 0.9 * len(q):
        return None
    pos = w0 + h.w_end + (len(q) - 1 - h.q_end)
    return (pos, h.identity) if pos <= g.length(m.prefix, m.chrom) else None


def complete_copy_from(m: Member, model: str, pos: int, g: Genomes, p: Params) -> bool:
    """Whether a complete copy of `model` starts at `pos` (a neighbour, not a remainder).

    Beyond a large insertion only the LTR's outer part remains before the insertion
    takes over, so the model cannot be placed there whole; a solo LTR or a tandem
    copy holds the whole model, inner end included.
    """
    w, w0 = g.fetch(m.prefix, m.chrom, pos - 2, pos + len(model) + 200)
    if len(w) < len(model):
        return False
    a = glocal(model, w)
    return (abs(w0 + a.start - pos) <= 2 and a.id_all >= p.min_whole_identity
            and a.id_last >= p.min_identity)


def complete_copy_to(m: Member, model: str, pos: int, g: Genomes, p: Params) -> bool:
    """`complete_copy_from` for a copy that ends at `pos`."""
    w, w0 = g.fetch(m.prefix, m.chrom, pos - len(model) - 200, pos + 2)
    if len(w) < len(model):
        return False
    a = glocal(model, w)
    return (abs(w0 + a.end - pos) <= 2 and a.id_all >= p.min_whole_identity
            and a.id_first >= p.min_identity)


def _with_anchor_left(m: Member, model: str, s: Side, g: Genomes, p: Params) -> Side:
    if not p.anchor or s.id_outer >= p.min_identity:
        return s
    a = anchor_left(m, model, g, p)
    if (a is None or a[1] < p.min_identity or a[0] > s.outer + 2
            or complete_copy_from(m, model, a[0], g, p)):
        return s
    return Side(a[0], s.inner, a[1], s.whole, "anchor")


def _with_anchor_right(m: Member, model: str, s: Side, g: Genomes, p: Params) -> Side:
    if not p.anchor or s.id_outer >= p.min_identity:
        return s
    a = anchor_right(m, model, g, p)
    if (a is None or a[1] < p.min_identity or a[0] < s.outer - 2
            or complete_copy_to(m, model, a[0], g, p)):
        return s
    return Side(a[0], s.inner, a[1], s.whole, "anchor")


def _median(xs: Sequence[float]) -> float:
    if not xs:
        return math.nan
    s = sorted(xs)
    return s[(len(s) - 1) // 2]


def _decide(m: Member, sides: List[Side], called: int, outward: int, g: Genomes,
            p: Params):
    """(end, src, id, support, raw, inner, reason, reportable) for one side.

    `outward` is -1 on the left (a smaller coordinate is further out), +1 on the right.
    """
    def moves(pos: Optional[int]) -> bool:
        if pos is None:
            return False
        d = called - pos if outward < 0 else pos - called      # + = further out
        return d >= p.min_ext or (d < 0 and -d <= p.max_trim)

    raw, _ = vote([s.outer for s in sides], p.agree, 1)
    passing = [s for s in sides if s.id_outer >= p.min_identity
               and s.whole >= p.min_whole_identity]
    voted, sup = vote([s.outer for s in passing], p.agree, p.min_support)
    inner = int(_median([s.inner for s in sides])) if sides else called
    if not moves(voted):
        reason = "."
        if moves(raw):
            reason = "gate_identity" if not passing else "no_support"
        return called, ".", math.nan, sup, raw if raw is not None else called, inner, \
            reason, moves(raw)
    carried = [s for s in passing if abs(s.outer - voted) <= p.agree]
    lo, hi = (voted, called - 1) if outward < 0 else (called + 1, voted)
    if hi >= lo and "N" in g.fetch(m.prefix, m.chrom, lo, hi)[0]:
        return called, ".", math.nan, sup, raw, inner, "gap_in_added", True
    srcs = sorted(s.src for s in carried)
    return (voted, srcs[len(srcs) // 2], _median([s.id_outer for s in carried]), sup, raw,
            int(_median([s.inner for s in carried])), ".", True)


def propose(m: Member, models: Sequence[Tuple[str, str, str]], g: Genomes, p: Params,
            orient: Optional[str] = None) -> Optional[Proposal]:
    """New outer ends for `m` from `(template, left_model, right_model)` triples.

    None when no side has anything to report: neither the vote nor the ungated
    median would move an end under the rules. Otherwise a Proposal whose `left`,
    `right` are the committed ends and whose reasons say why a side stayed.
    """
    lefts: List[Side] = []
    rights: List[Side] = []
    for _name, lmod, rmod in models:
        s = place_left(m, lmod, g)
        if s is not None:
            lefts.append(_with_anchor_left(m, lmod, s, g, p))
        s = place_right(m, rmod, g)
        if s is not None:
            rights.append(_with_anchor_right(m, rmod, s, g, p))
    L = _decide(m, lefts, m.start, -1, g, p)
    R = _decide(m, rights, m.end, +1, g, p)
    if not (L[7] or R[7]):
        return None
    frame = m.strand if m.stranded else (orient if orient in ("+", "-") else "+")
    return Proposal(member=m, orient=frame, left=L[0], right=R[0], left_src=L[1],
                    right_src=R[1], id_left=L[2], id_right=R[2], sup_left=L[3],
                    sup_right=R[3], raw_left=L[4], raw_right=R[4], inner_left=L[5],
                    inner_right=R[5], reason_left=L[6], reason_right=R[6],
                    templates=tuple(name for name, _, _ in models))


def _gap_near(q: str, t: str, c: int) -> int:
    best, n = 0, len(q)
    for i in range(max(0, c - 15), min(n, c + 16)):
        for s in (q, t):
            if s[i] == "-":
                j = i
                while j < n and s[j] == "-":
                    j += 1
                k = i
                while k > 0 and s[k - 1] == "-":
                    k -= 1
                best = max(best, j - k)
    return best


def _mismatch(q: str, t: str, cols) -> float:
    m = u = 0
    for j in cols:
        if 0 <= j < len(q) and q[j] != "-" and t[j] != "-":
            u += 1
            m += q[j] != t[j]
            if u == 30:
                break
    return m / u if u else 0.0


def obstacle(m: Member, prop: Proposal, g: Genomes) -> Tuple[str, str]:
    """What sits between the element's two LTRs at each old outer end that moved out.

    The new left LTR (outer end to the placements' inner end) is aligned to the new
    right LTR; the old end is mapped onto that alignment. Reported, never used
    to decide anything: 'gap:<bp>' for an indel >= 5 bp within 15 columns,
    'mm:<rate>' for >= 20% mismatches in the 30 bp just outside, else 'none'.
    """
    a, a0 = g.fetch(m.prefix, m.chrom, prop.left, prop.inner_left)
    b, b0 = g.fetch(m.prefix, m.chrom, prop.inner_right, prop.right)
    if not a or not b:
        return ".", "."
    r = parasail.nw_trace_scan_sat(a, b, 10, 1, matrix())
    q, t = r.traceback.query, r.traceback.ref
    col_a, col_b, ia, ib = {}, {}, a0, b0
    for j, (x, y) in enumerate(zip(q, t)):
        if x != "-":
            col_a[ia] = j
            ia += 1
        if y != "-":
            col_b[ib] = j
            ib += 1

    def label(col: Optional[int], outward: int) -> str:
        if col is None:
            return "."
        gap = _gap_near(q, t, col)
        if gap >= 5:
            return f"gap:{gap}"
        cols = range(col - 1, -1, -1) if outward < 0 else range(col + 1, len(q))
        mm = _mismatch(q, t, cols)
        return f"mm:{mm:.2f}" if mm >= 0.2 else "none"

    left = label(col_a.get(m.start), -1) if prop.left < m.start else "."
    right = label(col_b.get(m.end), +1) if prop.right > m.end else "."
    return left, right


def tsd_at(api, g: Genomes, m: Member, left: int, right: int, shift: int = 0) -> str:
    """Kmer2LTR's TSD call for an element spanning left..right; `shift` displaces
    the right flank."""
    a, a0 = g.fetch(m.prefix, m.chrom, left - 10, left + 9)
    b, _ = g.fetch(m.prefix, m.chrom, right - 9 + shift, right + 10 + shift)
    if a0 != left - 10 or len(a) != 20 or len(b) != 20:
        return "NA"
    hit = api.find_tsd(a + b, 10, 30, api.TSD_K, api.TSD_SHIFTS)
    return hit[0] if hit else "."

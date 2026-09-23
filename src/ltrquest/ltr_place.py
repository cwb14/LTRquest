"""Placing a family LTR model on one element: where would its LTRs really end?

The model is aligned end to end in a window around each called LTR (the core).
When the core's outer bases do not match -- a large insertion sits between the
true end and the rest of the LTR -- the model's outermost `anchor_len` bases are
searched for locally up to `max_indel` beyond the call (the anchor).

Only outward moves are ever proposed. An end that would move inward stays where
the caller put it: in the feasibility spike, trimming was wrong most of the time.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import parasail

from .ltr_model import (GAP_EXTEND, GAP_OPEN, Genomes, Member, Model, canonical_kmers, glocal,
                        jaccard, matrix, oriented_ltrs, rc)


@dataclass(frozen=True)
class Params:
    min_identity: float = 0.8       # outer-TERM identity a moved end must reach
    min_whole_identity: float = 0.6
    min_ext: int = 5                # an end must move at least this far to count
    anchor: bool = True
    anchor_len: int = 30
    max_indel: int = 5000
    n_models: int = 3               # models tried per element when a family has several
    combine: str = "best"      # best | median of the top n_models placements


@dataclass(frozen=True)
class Core:
    left: int
    l1: int
    r0: int
    right: int
    id_left: float
    id_right: float
    whole: float
    score: int


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
    model_id: str
    orient: str
    left: int
    right: int
    left_src: str       # core | anchor | '.' (did not move)
    right_src: str
    id_left: float
    id_right: float
    whole: float
    raw_left: int       # before the gates, for reporting
    raw_right: int
    inner_left: int     # the model's inner ends, for the obstacle report
    inner_right: int
    gate_ok: bool

    @property
    def ext_left(self) -> int:
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

    @property
    def raw_ext5(self) -> int:
        left, right = self.member.start - self.raw_left, self.raw_right - self.member.end
        return max(0, left if self.orient == "+" else right)

    @property
    def raw_ext3(self) -> int:
        left, right = self.member.start - self.raw_left, self.raw_right - self.member.end
        return max(0, right if self.orient == "+" else left)


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


def core(m: Member, fwd: str, g: Genomes) -> Optional[Core]:
    """The model (forward frame) glocal-aligned around each called LTR.

    Each window is widened by the truncation on either LTR, so it can reach a
    neighbouring copy of the same family -- a solo LTR, an uncalled fragment, a
    tandem array member -- which, being full length, outscores the truncated call
    site and passes every identity gate. A placement that does not overlap the LTR
    it was placed on is not a placement of this element, and nothing downstream
    would catch it: `conflict()` knows only about called elements.
    """
    lc = len(fwd)
    extra = max(100, lc - min(m.len_left, m.len_right) + 100)
    wl, wl0 = g.fetch(m.prefix, m.chrom, m.start - extra, m.l1 + extra)
    wr, wr0 = g.fetch(m.prefix, m.chrom, m.r0 - extra, m.end + extra)
    if len(wl) < lc // 2 or len(wr) < lc // 2:
        return None
    a, b = glocal(fwd, wl), glocal(fwd, wr)
    left, l1, r0, right = wl0 + a.start, wl0 + a.end, wr0 + b.start, wr0 + b.end
    if left > m.l1 or l1 < m.start or r0 > m.end or right < m.r0:
        return None
    return Core(left, l1, r0, right, a.id_first, b.id_last,
                min(a.id_all, b.id_all), a.score + b.score)


def anchor_left(m: Member, fwd: str, g: Genomes, p: Params) -> Optional[Tuple[int, float]]:
    q = fwd[:p.anchor_len]
    w, w0 = g.fetch(m.prefix, m.chrom, m.start - p.max_indel - len(q), m.start + 200)
    h = local(q, w)
    if h is None or h.q_start > 2 or (h.q_end - h.q_start + 1) < 0.9 * len(q):
        return None
    return w0 + h.w_start - h.q_start, h.identity


def anchor_right(m: Member, fwd: str, g: Genomes, p: Params) -> Optional[Tuple[int, float]]:
    q = fwd[-p.anchor_len:]
    w, w0 = g.fetch(m.prefix, m.chrom, m.end - 200, m.end + p.max_indel + len(q))
    h = local(q, w)
    if h is None or h.q_end < len(q) - 3 or (h.q_end - h.q_start + 1) < 0.9 * len(q):
        return None
    return w0 + h.w_end + (len(q) - 1 - h.q_end), h.identity


def candidate_models(m: Member, models: Sequence[Model], g: Genomes, n: int) -> List[Model]:
    """With several models, the `n` whose k-mers best match the element's called left LTR.

    A model `m` itself helped build is never one of them, whatever the method and
    whatever the caller passed: an element that is its own reference scores 1.0
    against itself, which would leave the identity gate -- the whole safeguard --
    inoperative for exactly the copies that built the model.
    """
    usable = [mo for mo in models if m.uid not in mo.refs]
    if len(usable) <= 1:
        return usable
    s, _ = g.fetch(m.prefix, m.chrom, m.start, m.l1)
    ks = canonical_kmers(s)
    return sorted(usable, key=lambda mo: (-jaccard(ks, mo.kmers), mo.model_id))[:n]


def _argmedian(cands, key):
    """The candidate holding the median coordinate -- a real placement, gates and all.

    Splicing a median coordinate into the best-scoring placement's `Core` (what
    `dataclasses.replace` used to do here) left the identity that authorises the
    move, the inner ends and the anchor comparison describing a different
    placement from the one committed.
    """
    order = sorted(cands, key=lambda t: (key(t[0]), t[1].model_id))
    return order[(len(order) - 1) // 2]


def propose(m: Member, models: Sequence[Model], g: Genomes, p: Params) -> Optional[Proposal]:
    cands = []
    for model in candidate_models(m, models, g, p.n_models):
        top = None
        for o in ((m.strand,) if m.stranded else ("+", "-")):
            fwd = model.seq if o == "+" else rc(model.seq)
            c = core(m, fwd, g)
            if c is not None and (top is None or c.score > top[0].score):
                top = (c, model, o, fwd)
        if top is not None:
            cands.append(top)
    if not cands:
        return None
    cands.sort(key=lambda t: -t[0].score)
    best = cands[0]
    at_left = at_right = best
    if p.combine == "median":
        # one orientation at a time: a median across two frames is not a coordinate
        same = [t for t in cands if t[2] == best[2]]
        if len(same) >= 3:
            at_left = _argmedian(same, lambda c: c.left)
            at_right = _argmedian(same, lambda c: c.right)
    cl, model_l, o, fwd_l = at_left        # the placement each committed end comes from,
    cr, model_r, _, fwd_r = at_right       # and whose gates therefore decide it

    left, left_id, left_src = cl.left, cl.id_left, "core"
    if p.anchor and left_id < p.min_identity:
        a = anchor_left(m, fwd_l, g, p)
        if a is not None and a[1] >= p.min_identity and a[0] <= cl.left + 2:
            left, left_id, left_src = a[0], a[1], "anchor"
    right, right_id, right_src = cr.right, cr.id_right, "core"
    if p.anchor and right_id < p.min_identity:
        a = anchor_right(m, fwd_r, g, p)
        if a is not None and a[1] >= p.min_identity and a[0] >= cr.right - 2:
            right, right_id, right_src = a[0], a[1], "anchor"

    moves_left = left <= m.start - p.min_ext
    moves_right = right >= m.end + p.min_ext
    if not (moves_left or moves_right):
        return None
    whole = min(cl.whole, cr.whole)
    whole_ok = whole >= p.min_whole_identity
    gate_left = moves_left and left_id >= p.min_identity and whole_ok
    gate_right = moves_right and right_id >= p.min_identity and whole_ok
    model_id = (model_l.model_id if model_l is model_r
                else f"{model_l.model_id},{model_r.model_id}")
    return Proposal(
        member=m, model_id=model_id, orient=o,
        left=left if gate_left else m.start, right=right if gate_right else m.end,
        left_src=left_src if gate_left else ".", right_src=right_src if gate_right else ".",
        id_left=left_id, id_right=right_id, whole=whole, raw_left=left, raw_right=right,
        inner_left=cl.l1, inner_right=cr.r0, gate_ok=gate_left or gate_right)


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
    """What sits between the element's two LTRs at each old outer end that moved.

    The new left LTR (outer end to the model's inner end) is aligned to the new
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


def nearest_templates(consensus: Model, refs: Sequence[Member], g: Genomes,
                      tol: int = 2) -> List[Model]:
    """References whose called ends the family consensus confirms, as models of their own.

    A template's termini are then checked by its family rather than taken on trust --
    both ends of both LTRs, not just the outer ones. A template is cut at the
    reference's *inner* ends (`m.start..m.l1`), and placing it on a target sets that
    target's opposite *outer* end: a reference over-called a few bases inward pushes
    every target it places that far past the true terminus, which the outer-30 gate
    tolerates (~6 bases) and which silently destroys the TSD the family QC looks for.
    Falls back to the consensus itself when no reference qualifies.
    """
    out: List[Model] = []
    for m in refs:
        fwd = consensus.seq if m.strand != "-" else rc(consensus.seq)
        c = core(m, fwd, g)
        if c is None or abs(c.left - m.start) > tol or abs(c.right - m.end) > tol:
            continue
        if abs(c.l1 - m.l1) > tol or abs(c.r0 - m.r0) > tol:
            continue
        five, _ = oriented_ltrs(m, g)
        out.append(Model(f"{consensus.family}:tpl:{m.key}", consensus.family, five,
                         consensus.modal_len, 1, frozenset({m.uid})))
    return out or [consensus]

"""Rejoining a host whose LTR pair a tenant split (a step of ltrquest.reconcile).

An element that inserts at position p of another element's LTR leaves
`[1..p] TENANT [p+1..L]`. LTRharvest and LTR_FINDER pair the near-side fragments,
the `[1..p]` of both LTRs, so the host is called truncated at the tenant, without a
TSD, and the tenant looks free beside it. The full host is out of the callers'
reach: its LTR alignment is broken, and with the tenant inside it is usually longer
than their maximum LTR distance. Masking between rounds does not help either, as
host and tenant are mostly called in the same round.

The rest of the split LTR is still there, just beyond the tenant. For two calls
that abut, this asks whether the bases beyond one continue the other's own LTR,
which they must from the junction on: the bases after the tenant against those
after the host's called left LTR (tenant in the right LTR), or, mirrored, the bases
before the tenant against those before the host's called right LTR (tenant in the
left LTR). Where one side and only one does, the host is rebuilt to the end of that
continuation with the tenant inside its LTR. Its divergence columns are Kmer2LTR's,
measured on the called fragment pair, which is a sample of the same LTR pair.

Coordinates are genome-forward, 1-based, inclusive, as in the depth tables.
"""

from __future__ import annotations

import math
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import ends
from .detect import _ltrs_shared
from .ltr_model import Member
from .ltr_place import local, tsd_call
from .reboundary_io import SpanIndex, conflict, forward, record_matches, stored
from .recover_strand import Genome
from .table import Columns

# Fixed by biology or convention, not tuned to any genome.
JUNCTION = 20           # bp: a TSD (4-6 bp) plus the few bases an x-drop LTR
                        # extension stops short of an insertion junction
MIN_IDENTITY = 0.70     # the callers' own LTR-similarity floor (ltrharvest -similar 70)
MAX_EVALUE = 1e-5       # the usual BLAST significance cutoff
LAMBDA, K = 0.55, 0.21  # Karlin-Altschul parameters, lambda kept below the +2/-3
                        # value (0.63 ungapped), which makes the E-value err high
MAX_LTR = 7000          # bp: the longest LTR the callers may report (-maxlenltr /
                        # -L 7000); a continuation past it is a tandem copy, not an LTR
FIRST_WINDOW = 256      # bp aligned on each side first, doubled while the alignment
MAX_WINDOW = 8192       # runs to the window's edge, up to past any LTR the callers
                        # report (-maxlenltr / -L 7000)


class Bases:
    """`ltr_model.Genomes.fetch` on one genome, gzip-safe, its index in a scratch dir."""

    def __init__(self, path: str, workdir: str):
        self._g = Genome(path, workdir)

    def fetch(self, _prefix: str, chrom: str, lo: int, hi: int) -> Tuple[str, int]:
        try:
            n = len(self._g.fa[chrom])
        except KeyError:
            return "", lo
        a, b = max(1, lo), min(n, hi)
        if a > b:
            return "", a
        return self._g.get(chrom, a, b) or "", a

    def names(self) -> List[str]:
        return list(self._g.fa.keys())

    def close(self) -> None:
        self._g.close()


@dataclass(frozen=True)
class Continuation:
    """The split LTR's remainder, aligned from the junction on."""
    q_far: int          # its outermost base, beyond the tenant: the host's true end
    t_far: int          # the host's LTR inner end it reaches on the other side
    identity: float
    score: int
    evalue: float


@dataclass(frozen=True)
class Rejoin:
    host: dict          # reconcile records, as called
    tenant: dict
    old_key: str        # the host's key as called
    ltr: str            # which of the host's LTRs holds the tenant: right | left
    start: int          # the rebuilt host
    end: int
    l1: int
    r0: int
    cont: Continuation

    @property
    def key(self) -> str:
        return f"{self.host['chrom']}:{self.start}-{self.end}"

    @property
    def remainder(self) -> int:
        return (self.end - self.tenant["e"] if self.ltr == "right"
                else self.tenant["s"] - self.start)

    @property
    def ltr_bases(self) -> Tuple[int, int]:
        """(left, right) LTR lengths of the rebuilt host, the tenant's bases left out."""
        t = self.tenant
        left, right = self.l1 - self.start + 1, self.end - self.r0 + 1
        lo, hi = (self.r0, self.end) if self.ltr == "right" else (self.start, self.l1)
        in_ltr = min(t["e"], hi) - max(t["s"], lo) + 1
        return (left, right - in_ltr) if self.ltr == "right" else (left - in_ltr, right)


def continuation(g: Bases, chrom: str, q_anchor: int, t_anchor: int, t_limit: int,
                 leftward: bool) -> Optional[Continuation]:
    """Does the sequence from `q_anchor` on continue the LTR from `t_anchor` on?

    Rightward, both run to higher coordinates and `t_limit` is the last base the
    LTR may reach (the host's right-LTR start, less one). Leftward, both run to
    lower coordinates and `t_limit` is the first. The alignment must start within
    JUNCTION of both anchors, on one diagonal.
    """
    w = FIRST_WINDOW
    while True:
        if leftward:
            q, q0 = g.fetch("", chrom, q_anchor - w + 1, q_anchor + JUNCTION)
            t, t0 = g.fetch("", chrom, max(t_limit, t_anchor - w + 1), t_anchor + JUNCTION)
            q_last, t_last = q0 + len(q) - 1, t0 + len(t) - 1
            q, t = q[::-1], t[::-1]
            qa, ta = q_last - q_anchor, t_last - t_anchor
            room_q, room_t = q0 > 1, t0 > t_limit
        else:
            q, q0 = g.fetch("", chrom, q_anchor - JUNCTION, q_anchor + w - 1)
            t, t0 = g.fetch("", chrom, t_anchor - JUNCTION, min(t_limit, t_anchor + w - 1))
            qa, ta = q_anchor - q0, t_anchor - t0
            room_q, room_t = len(q) - qa == w, t0 + len(t) - 1 < t_limit
        if len(q) <= qa + JUNCTION or len(t) <= ta + JUNCTION:
            return None
        h = local(q, t)
        if h is None:
            return None
        at_edge_q = h.q_end >= len(q) - JUNCTION
        at_edge_t = h.w_end >= len(t) - JUNCTION
        if w < MAX_WINDOW and ((at_edge_q and room_q) or (at_edge_t and room_t)):
            w *= 2
            continue
        break
    dq, dt = h.q_start - qa, h.w_start - ta
    if abs(dq) > JUNCTION or abs(dt) > JUNCTION or abs(dq - dt) > JUNCTION:
        return None
    evalue = K * len(q) * len(t) * math.exp(-LAMBDA * h.score)
    if h.identity < MIN_IDENTITY or evalue > MAX_EVALUE:
        return None
    if leftward:
        return Continuation(q_last - h.q_end, t_last - h.w_end, h.identity, h.score, evalue)
    return Continuation(q0 + h.q_end, t0 + h.w_end, h.identity, h.score, evalue)


def neighbours(survivors: Sequence[dict]) -> List[Tuple[dict, dict]]:
    """(left, right) calls that abut within JUNCTION, neither inside the other."""
    out = []
    by_chrom: Dict[str, List[dict]] = defaultdict(list)
    for e in survivors:
        by_chrom[e["chrom"]].append(e)
    for arr in by_chrom.values():
        arr.sort(key=lambda x: (x["s"], -x["e"]))
        for i in range(1, len(arr)):
            y = arr[i]
            j = i - 1
            while j >= 0 and arr[j]["e"] >= y["e"]:      # skip calls enclosing y
                j -= 1
            if j >= 0 and arr[j]["s"] < y["s"] and abs(y["s"] - arr[j]["e"] - 1) <= JUNCTION:
                out.append((arr[j], y))
    return out


def propose(x: dict, y: dict, bounds: Dict[str, Tuple[int, int, int, int]],
            g: Bases) -> Tuple[Optional[Rejoin], str]:
    """The rejoin one abutting pair supports, or None and why not."""
    xb, yb = bounds.get(x["key"]), bounds.get(y["key"])
    if xb is None or yb is None:
        return None, "no_bounds"
    chrom = x["chrom"]
    # y in x's right LTR: past y continues x's left LTR, from its called inner end
    right = continuation(g, chrom, y["e"] + 1, xb[1] + 1, xb[2] - 1, leftward=False)
    # x in y's left LTR: before x continues y's right LTR, back from its called inner end
    left = continuation(g, chrom, x["s"] - 1, yb[2] - 1, yb[1] + 1, leftward=True)
    if right and left:
        return None, "ambiguous"
    if right and right.q_far > y["e"] and right.t_far >= xb[1]:
        return Rejoin(x, y, x["key"], "right", x["s"], right.q_far, right.t_far, xb[2],
                      right), ""
    if left and left.q_far < x["s"] and left.t_far <= yb[2]:
        return Rejoin(y, x, y["key"], "left", left.q_far, y["e"], yb[1], left.t_far,
                      left), ""
    return None, "no_continuation"


def _member(rec: dict, b: Tuple[int, int, int, int]) -> Member:
    return Member("", rec["col1"], rec["chrom"], rec["s"], rec["e"], b[1], b[2], ".", ".",
                  ".", 0, None, ".")


def rejoin_split_hosts(survivors: List[dict], bounds: Dict[str, Tuple[int, int, int, int]],
                       records: Dict[str, str], cols: Columns, g: Bases,
                       genome_api) -> Tuple[List[Rejoin], Counter]:
    """Rebuild every host a tenant in one LTR split, in place.

    Rewrites the host's record dict (key, span, row), its FASTA record in `records`,
    and `bounds`. A rejoin is not made where it is ambiguous, where a rebuilt
    LTR would be longer than any the callers report, where an element would be in
    two, where it would swallow or cut through another call (the rules
    re-boundarying moves an end by) or another rebuilt host, where it would make
    the host share an LTR it did not before, or where the host's stored record is
    not its genome span. Returns the rejoins made, and why the others were not.
    """
    why: Counter = Counter()
    proposals = []
    for x, y in neighbours(survivors):
        rj, reason = propose(x, y, bounds, g)
        if rj is None:
            why[reason] += 1
        else:
            proposals.append(rj)
    members = [_member(e, bounds[e["key"]]) for e in survivors if e["key"] in bounds]
    index = SpanIndex(members)
    by_key = {m.key: m for m in members}
    used = set()
    rebuilt: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    done: List[Rejoin] = []
    for rj in sorted(proposals, key=lambda r: -r.cont.score):
        if id(rj.host) in used or id(rj.tenant) in used:      # keys change on a rejoin
            why["element_in_two"] += 1
            continue
        if max(rj.ltr_bases) > MAX_LTR:
            why["ltr_too_long"] += 1
            continue
        if any(s < rj.start <= e < rj.end or rj.start < s <= rj.end < e
               for s, e in rebuilt[rj.host["chrom"]]):
            why["overlaps_rejoin"] += 1
            continue
        reason = _refusal(rj, index, by_key, bounds)
        fwd = None if reason else _new_record(rj, records, cols, g)
        if fwd is None:
            why[reason or "record_mismatch"] += 1
            continue
        _apply(rj, fwd, bounds, records, cols, g, genome_api)
        used.update((id(rj.host), id(rj.tenant)))
        rebuilt[rj.host["chrom"]].append((rj.start, rj.end))
        done.append(rj)
    return done, why


def _refusal(rj: Rejoin, index: SpanIndex, by_key: Dict[str, Member], bounds) -> str:
    h, t = rj.host, rj.tenant
    inside = frozenset(k for s, e, k in index.overlapping("", t["chrom"], t["s"], t["e"])
                       if s >= t["s"] and e <= t["e"])
    clash = conflict(index, by_key[h["key"]], rj.start, rj.end, ignore=inside)
    if clash:
        return clash
    old, new = bounds[h["key"]], (rj.start, rj.l1, rj.r0, rj.end)
    for _s, _e, k in index.overlapping("", h["chrom"], rj.start, rj.end):
        if k == h["key"] or k in inside:
            continue
        pair = {"old": old, "new": new, "z": bounds[k]}
        if _ltrs_shared("new", "z", pair) and not _ltrs_shared("old", "z", pair):
            return "shares_ltr"
    return ""


def _orientation(rec: dict, cols: Columns) -> str:
    return cols.get(rec["cols"], "orientation")


def _new_record(rj: Rejoin, records: Dict[str, str], cols: Columns,
                g: Bases) -> Optional[str]:
    """The rebuilt host's record, forward: its stored record (round masks and all)
    with the genome spliced on; None if the stored record is not its span."""
    h, chrom = rj.host, rj.host["chrom"]
    rec = records.get(h["col1"])
    if rec is None:
        return None
    fwd = forward(rec, _orientation(h, cols))
    if not record_matches(fwd, g.fetch("", chrom, h["s"], h["e"])[0]):
        return None
    if rj.ltr == "right":
        fwd += g.fetch("", chrom, h["e"] + 1, rj.end)[0]
    else:
        fwd = g.fetch("", chrom, rj.start, h["s"] - 1)[0] + fwd
    return fwd if len(fwd) == rj.end - rj.start + 1 else None


def _apply(rj: Rejoin, fwd: str, bounds, records: Dict[str, str], cols: Columns,
           g: Bases, genome_api) -> None:
    h, chrom = rj.host, rj.host["chrom"]
    s, e = rj.start, rj.end
    left_bp, right_bp = rj.ltr_bases
    tsd, offset = tsd_call(genome_api, g, _member(h, bounds[h["key"]]), s, e)
    canon = ends.tsd_between(lambda lo, hi: g.fetch("", chrom, lo, hi), s, e)
    if canon:
        tsd, offset = canon, "0,0"
    suffix = h["col1"].split("#", 1)[1] if "#" in h["col1"] else None
    name = f"{chrom}:{s}-{e}" + (f"#{suffix}" if suffix is not None else "")
    row = list(h["cols"])

    def put(col: str, value) -> None:
        i = cols.index(col)
        if i is not None and i < len(row):
            row[i] = str(value)

    put("seq_id", name)
    put("seq_len", e - s + 1)
    put("ltr5_end", rj.l1 - s + 1)
    put("ltr3_start", rj.r0 - s + 1)
    put("ltr3_end", e - s + 1)
    put("ltr5_len", left_bp)
    put("ltr3_len", right_bp)
    put("motif", f"{fwd[:2]}...{fwd[-2:]}".lower())
    put("tsd", tsd)
    put("tsd_input", tsd)
    put("tsd_offset", offset)
    records[name] = stored(fwd, _orientation(h, cols))
    bounds[rj.key] = (s, rj.l1, rj.r0, e)
    h.update(key=rj.key, s=s, e=e, col1=name, cols=row, line="\t".join(row))


def sidecar_rows(rejoins: Sequence[Rejoin], cols: Columns) -> List[List[str]]:
    """annotate.REJOIN_HEADER rows, one per rebuilt host."""
    return [[rj.old_key, rj.key, rj.tenant["key"], rj.ltr, str(rj.remainder),
             f"{rj.cont.identity:.4f}", f"{rj.cont.evalue:.3g}",
             cols.get(rj.host["cols"], "tsd")] for rj in rejoins]


def run(survivors: List[dict], bounds, records: Dict[str, str], cols: Columns,
        genome: str, genome_api, scratch: str,
        log: Callable[[str], None]) -> List[Rejoin]:
    """The reconcile step: open the genome in a scratch dir, rejoin, report."""
    workdir = tempfile.mkdtemp(prefix=".rejoin_", dir=scratch)
    g = Bases(genome, workdir)
    try:
        # every call was read off this genome, so a sequence it lacks means the
        # wrong one was given, which would otherwise rebuild nothing, silently
        absent = sorted({e["chrom"] for e in survivors} - set(g.names()))
        if absent:
            raise RuntimeError(
                f"--genome {genome} has no sequence named {', '.join(absent[:5])}"
                f"{' ...' if len(absent) > 5 else ''}; give the genome the rounds were run on")
        done, why = rejoin_split_hosts(survivors, bounds, records, cols, g, genome_api)
    finally:
        g.close()
        shutil.rmtree(workdir, ignore_errors=True)
    n_right = sum(rj.ltr == "right" for rj in done)
    n_pairs = len(done) + sum(why.values())
    not_made = ", ".join(f"{k}={v}" for k, v in sorted(why.items())
                         if k not in ("no_continuation", "no_bounds"))
    log(f"{n_pairs} abutting pair(s); rejoined {len(done)} host(s) split by a tenant "
        f"in one LTR (right LTR {n_right}, left LTR {len(done) - n_right})"
        + (f"; not rejoined although the LTR continues: {not_made}" if not_made else ""))
    return done

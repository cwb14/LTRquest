"""Template-guided re-boundarying of LTR-RT calls (opt-in stage).

LTRharvest and LTR_FINDER extend an LTR pair outward from a seed and stop at the
first obstacle between the element's two LTRs -- an indel, a patch of mutations --
so such an element is called short; a few are called a few bases long. This stage
re-bounds every element whose called ends carry no target-site duplication (TSD).
Its templates are the elements whose called ends carry an exact one, pooled over
every genome and family: their outer ends are where the insertion put them. Each
target's nearest templates (BLASTN) are placed on it one LTR per side, and where at
least two agree an outer end moves out, or in by a few bases. A second call of the
same element that the move reaches is merged away. Kmer2LTR then re-measures the
pair the element was called with, on the new record, with the templates' ends as
external evidence: the outer ends stay where the templates put them, the other
LTR's inner end follows real homology, and bases with no partner are gaps, so the
divergence is not inflated. Only the `_clean_` tables and FASTAs are rewritten;
`<prefix>_reboundary.tsv` records every candidate and why it was or was not moved,
and doubles as the old->new key map for the GFF3.

Usage:
  ltrquest-reboundary --indir RUN --prefix P1 [P2 ...] --genome G1 [G2 ...] [options]
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from multiprocessing import Pool
from typing import Dict, List, Optional, Sequence, Tuple

from . import kmer2ltr as k2l
from . import templates as tp
from .kmer2ltr import COLUMNS
from .ltr_model import Genomes, Member
from .ltr_place import Params, Proposal, obstacle, propose, tsd_at
from .reboundary_io import (
    SIDECAR_SUFFIX,
    Accepted,
    Commit,
    SpanIndex,
    conflict,
    fetch_records,
    forward,
    genome_order,
    hosts_of,
    leftover_message,
    leftover_staging,
    load_clean_tables,
    members_from,
    merge_partner,
    mutual_conflicts,
    nested_index,
    rewrite,
    sanitize,
    sidecar_text,
    stored,
    sync_dir,
    trim_restore,
)
from .reconcile import IUPAC_DEPTH_SEQ

# Bits of external evidence handed to Kmer2LTR that the record's termini are the
# element's ends. With the homology-paired credit it can only authorize those outer
# ends -- the templates chose them -- so it is a constant, not a knob.
CREDIT_BITS = 1e6
EXIT_INCOMPATIBLE_KMER2LTR = 3
_I = {c: i for i, c in enumerate(COLUMNS)}
_CIGAR = re.compile(r"(\d+)([=XID])")


def log(msg: str) -> None:
    print(f"[reboundary] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[reboundary] WARNING: {msg}", file=sys.stderr, flush=True)


@dataclass(frozen=True)
class Settings:
    mutation_rate: float = 3e-8
    threads: int = 8
    blastn: str = "blastn"
    blast_task: str = "blastn"
    max_templates_per_family: int = 100     # database cap; only `place.n_templates` are used
    merge_bp: int = 5                        # a merge partner's end within this of the new end
    place: Params = field(default_factory=Params)


_G: Optional[Genomes] = None
_API = None


def _init(paths: Dict[str, str], tools_dir: str) -> None:
    """Worker set-up: genome handles and the Kmer2LTR API, once per process."""
    global _G, _API
    _G = Genomes(paths)
    _API = k2l.api(tools_dir)


def _label(t: Member) -> str:
    return f"{t.prefix}:{t.name}"


def _flip(strand: str) -> str:
    return {"+": "-", "-": "+"}.get(strand, strand)


def implied_orient(m: Member, hits: Sequence[Tuple[Member, str]]) -> str:
    """The frame ext5/ext3 are reported in: `m`'s strand, else the one its stranded
    templates imply (majority; a tie reads `+`), else `+`."""
    if m.stranded:
        return m.strand
    votes = Counter(t.strand if rel == "+" else _flip(t.strand) for t, rel in hits if t.stranded)
    if not votes:
        return "+"
    return sorted(votes.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def place_job(args):
    """Place `m`'s templates on it; for a move, the TSD at the new ends and its null."""
    m, hits, p = args
    models = [(_label(t), *tp.models_for(t, rel, _G)) for t, rel in hits]
    prop = propose(m, models, _G, p, orient=implied_orient(m, hits))
    if prop is None:
        return m.uid, None, ".", ".", (".", ".")
    if not prop.gate_ok:
        return m.uid, prop, ".", ".", (".", ".")
    return (m.uid, prop, tsd_at(_API, _G, m, prop.left, prop.right),
            tsd_at(_API, _G, m, prop.left, prop.right, shift=1000), obstacle(m, prop, _G))


@dataclass(frozen=True)
class Verdict:
    uid: str
    status: str
    accepted: Optional[Accepted]
    k2l_status: str
    tsd: str
    k2p: str
    unpaired: Tuple[int, int] = (0, 0)      # (genomic left, genomic right) unpaired bases


def window(m: Member, left: int, right: int):
    """Kmer2LTR's four reference cuts around a record (genome._cuts), from pyfaidx."""
    pad, probe = _API.PAD, _API.PROBE
    return _API.Window(_G.fetch(m.prefix, m.chrom, left - pad, left - 1)[0],
                       _G.fetch(m.prefix, m.chrom, left, left + probe - 1)[0],
                       _G.fetch(m.prefix, m.chrom, right - probe + 1, right)[0],
                       _G.fetch(m.prefix, m.chrom, right + 1, right + pad)[0])


def rebase(fields: List[str], f5: int, name: str) -> List[str]:
    """Kmer2LTR's row restated against the record cut to its own bounds
    (detect.rebase_to_trimmed)."""
    out = list(fields)
    out[_I["seq_id"]] = name
    for c in ("ltr5_end", "ltr3_start", "ltr3_end"):
        out[_I[c]] = str(int(out[_I[c]]) - f5)
    out[_I["ltr5_start"]] = "1"
    out[_I["seq_len"]] = out[_I["ltr3_end"]]
    out[_I["flank5_len"]] = out[_I["flank3_len"]] = "0"
    return out


def length_ok(l5, l3, aln, left: int, right: int) -> bool:
    """Detection's Step 8a size floors on the re-measured pair: LTRs >= 100, alignment
    >= 90, end - start >= 300 (detection's own element test, so >= 301 bp). Its 0.65
    length-ratio term is left out: LTRs of unequal length are what a deletion in one
    or an insertion in the other looks like, and here the templates, not the pair,
    set the ends."""
    if l5 is None or l3 is None or aln is None:
        return False
    return min(l5, l3) >= 100 and aln >= 90 and right - left >= 300


def trim_spans(cigar: str, spans: Tuple[int, int, int, int], cut5: int,
               cut3: int) -> Tuple[int, int, int, int]:
    """The called pair after a trim. Cutting `cut5` bases off the left LTR's outer end
    takes the partner bases the called alignment paired with them off the right LTR's
    inner end; a `cut3` cut on the right takes theirs off the left LTR's inner end.
    `cigar` is the called pair's (query = left LTR, ref = right LTR)."""
    l5b, l5e, l3b, l3e = spans
    ops = [op for n, op in _CIGAR.findall(cigar or "") for _ in range(int(n))]
    if cut5 > 0:
        q = r = 0
        for op in ops:
            if q >= cut5:
                break
            q += op in "=XI"
            r += op in "=XD"
        l3b += r
    if cut3 > 0:
        q = r = 0
        for op in reversed(ops):
            if r >= cut3:
                break
            r += op in "=XD"
            q += op in "=XI"
        l5e -= q
    return l5b, l5e, l3b, l3e


def _unpaired(cigar: str) -> Tuple[int, int]:
    """(leading I run, trailing D run): left-LTR bases at its outer end, right-LTR
    bases at its outer end, with no partner."""
    ops = _CIGAR.findall(cigar or "")
    lead = int(ops[0][0]) if ops and ops[0][1] == "I" else 0
    tail = int(ops[-1][0]) if ops and ops[-1][1] == "D" else 0
    return lead, tail


def record_matches(fwd: str, genomic: str) -> bool:
    """Whether a stored record, turned to the forward frame, is the genome at its called
    span wherever it is not masked (a depth letter stands for a nested element). Case
    is ignored: detection keeps a soft-masked genome's lowercase."""
    if fwd == genomic:
        return True
    return len(fwd) == len(genomic) and all(
        a == b for a, b in zip(fwd.upper(), genomic.upper()) if a in "ACGT")


def arbitrate(args) -> Verdict:
    """Kmer2LTR re-measures the called pair on the new record; the templates' ends are
    its credited termini."""
    p, record, cigar, mu = args
    m = p.member
    fwd = forward(record, m.orientation)
    if not record_matches(fwd, _G.fetch(m.prefix, m.chrom, m.start, m.end)[0]):
        # e.g. a record stored reverse-complemented under an orientation `+` row:
        # splicing genome segments onto it would build a chimera
        return Verdict(m.uid, "record_mismatch", None, "NA", "NA", "NA")
    cut5, cut3 = max(0, p.left - m.start), max(0, m.end - p.right)
    left_seg = _G.fetch(m.prefix, m.chrom, p.left, m.start - 1)[0] if p.left < m.start else ""
    right_seg = _G.fetch(m.prefix, m.chrom, m.end + 1, p.right)[0] if p.right > m.end else ""
    new = left_seg + fwd[cut5:len(fwd) - cut3] + right_seg
    if len(new) != p.right - p.left + 1:
        return Verdict(m.uid, "record_mismatch", None, "NA", "NA", "NA")
    n = len(new)
    spans = trim_spans(cigar, (max(0, m.start - p.left), m.l1 - p.left, m.r0 - p.left,
                               min(n - 1, m.end - p.left)), cut5, cut3)
    suffix = "#" + m.name.split("#", 1)[1] if "#" in m.name else ""
    name = f"{m.chrom}:{p.left}-{p.right}{suffix}"
    seq = sanitize(new)
    ctx = _API.orient(seq, window(m, p.left, p.right))
    res = _API.classify(name, seq, period_rule="outermost", mutation_rate=mu,
                        tsd_credit=CREDIT_BITS, spans=spans)
    if ctx is not None:
        res = _API.annotate(res, seq, ctx, _API.Options())
    if res.status != "pass":
        return Verdict(m.uid, "kmer2ltr_not_pass", None, res.status, "NA", "NA")
    row = _API.format_row(res).split("\t")
    tsd, k2p = row[_I["tsd"]], row[_I["k2p"]]
    if res.flank5_len or res.flank3_len:       # cannot happen with credit + known spans
        return Verdict(m.uid, "kmer2ltr_moved_ends", None, res.status, tsd, k2p)
    if not length_ok(res.ltr5_len, res.ltr3_len, res.aln_len, p.left, p.right):
        return Verdict(m.uid, "length_filter", None, res.status, tsd, k2p)
    fields = rebase(row, 0, name)
    fields[_I["orientation"]] = m.orientation   # a storage fact: the record stays stored as it was
    acc = Accepted(m, fields[0], fields, stored(new, m.orientation), p.left, p.right)
    return Verdict(m.uid, "moved", acc, res.status, tsd, k2p, _unpaired(res.cigar))


def _fmt_id(x: float) -> str:
    return "." if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.3f}"


def sidecar_row(p: Proposal, tsd_new: str, tsd_null: str, obst: Tuple[str, str]) -> Dict[str, str]:
    m = p.member
    plus = p.orient == "+"

    def bio(left, right):
        return (left, right) if plus else (right, left)

    ext_l = p.ext_left if p.left != m.start else m.start - p.raw_left
    ext_r = p.ext_right if p.right != m.end else p.raw_right - m.end
    reason = next((r for r in (p.reason_left, p.reason_right) if r != "."), ".")
    ext5, ext3 = bio(ext_l, ext_r)
    src5, src3 = bio(p.left_src, p.right_src)
    id5, id3 = bio(p.id_left, p.id_right)
    sup5, sup3 = bio(p.sup_left, p.sup_right)
    ob5, ob3 = bio(*obst)
    return {
        "old_seq_id": m.name, "new_seq_id": ".", "family": m.family, "method": "templates",
        "templates": ",".join(p.templates) or ".", "decision": "rejected",
        "reason": "." if p.gate_ok else reason, "ext5": str(ext5), "ext3": str(ext3),
        "end_source5": src5, "end_source3": src3, "id_outer5": _fmt_id(id5),
        "id_outer3": _fmt_id(id3), "support5": str(sup5), "support3": str(sup3),
        "k2l_status": ".", "tsd_called": m.tsd or ".", "tsd_new": tsd_new, "tsd_null": tsd_null,
        "k2p_called": "NA" if m.k2p is None else f"{m.k2p:g}", "k2p_new": ".",
        "unpaired5": ".", "unpaired3": ".", "obstacle5": ob5, "obstacle3": ob3,
        "merged_into": ".",
    }


def _partner_row(x: Member, into: str) -> Dict[str, str]:
    row = dict.fromkeys(("new_seq_id", "templates", "ext5", "ext3", "end_source5", "end_source3",
                         "id_outer5", "id_outer3", "support5", "support3", "k2l_status",
                         "tsd_new", "tsd_null", "k2p_new", "unpaired5", "unpaired3",
                         "obstacle5", "obstacle3"), ".")
    row.update(old_seq_id=x.name, family=x.family, method="templates", decision="merged",
               reason="split_call", new_seq_id=into, merged_into=into, tsd_called=x.tsd or ".",
               k2p_called="NA" if x.k2p is None else f"{x.k2p:g}")
    return row


def _unmerge(rows: Dict[str, Dict[str, str]], pu: Optional[str]) -> None:
    """A survivor dropped out: the call it would have absorbed keeps its own call, and a
    row it has (it was a candidate itself) says so."""
    if pu is not None and pu in rows:
        rows[pu]["reason"] = "merge_partner_failed"


def _row_order(row: Dict[str, str]):
    chrom, span = row["old_seq_id"].split("#", 1)[0].rsplit(":", 1)
    start, end = (int(x) for x in span.split("-"))
    return chrom, start, end, row["old_seq_id"]


def check_prefixes(prefixes: Sequence[str]) -> None:
    """One prefix may be named once: twice would stage every one of its files twice.

    Each file's staging name is derived from the file, so a second rewrite of
    the same path swaps that path's first rewrite into the `.old` file holding
    its original. `Commit.open` refuses that too; this is the boundary where the
    user finds out, before a post-hoc backup has moved anything.
    """
    repeated = sorted(p for p, n in Counter(prefixes).items() if n > 1)
    if repeated:
        raise SystemExit(f"reboundary: --prefix {', '.join(repeated)} given more than once; "
                         f"name each genome once (they are pooled anyway). Nothing has "
                         f"been changed.")


def check_indir(indir: str) -> None:
    """Refuse to start while an interrupted rewrite's `.new`/`.old` files are in `indir`."""
    leftovers = leftover_staging(indir)
    if leftovers:
        raise SystemExit(leftover_message(indir, leftovers))


def check_blastn(blastn: str) -> str:
    """The blastn to run, with makeblastdb beside it (or on PATH); fail before any work."""
    path = shutil.which(blastn) or (blastn if os.path.isfile(blastn) else None)
    if path is None:
        raise SystemExit(f"reboundary: blastn not found ({blastn}); it is in environment.yml")
    here = os.path.dirname(path)
    if not (os.path.isfile(os.path.join(here, "makeblastdb")) or shutil.which("makeblastdb")):
        raise SystemExit("reboundary: makeblastdb not found beside blastn or on PATH")
    return path


def _dump(path: str, props: Dict[str, Proposal], rows: Dict[str, Dict[str, str]]) -> None:
    with open(path, "w") as fh:
        fh.write("uid\tleft\tright\traw_left\traw_right\tsup_left\tsup_right\t"
                 "reason_left\treason_right\tdecision\treason\n")
        for u, p in sorted(props.items()):
            r = rows[u]
            fh.write("\t".join(str(x) for x in (
                u.replace("\t", "|"), p.left, p.right, p.raw_left, p.raw_right, p.sup_left,
                p.sup_right, p.reason_left, p.reason_right, r["decision"], r["reason"])) + "\n")


def run(indir: str, prefixes: Sequence[str], genomes: Sequence[str], s: Settings,
        tools_dir: str, dump: Optional[str] = None, verbose: bool = False) -> Dict[str, int]:
    t0 = time.time()
    vlog = log if verbose else (lambda msg: None)
    check_prefixes(prefixes)
    check_indir(indir)
    if len(prefixes) != len(genomes):
        raise SystemExit("reboundary: --prefix and --genome must pair up one to one")
    blastn = check_blastn(s.blastn)
    for gpath in genomes:
        if not os.path.isfile(gpath):
            raise SystemExit(f"reboundary: genome not found: {gpath}")
    k2l.api(tools_dir)     # fail here, and clone at most once, before any worker starts
    tables_by_prefix = {p: load_clean_tables(indir, p) for p in prefixes}
    for p in prefixes:
        if not tables_by_prefix[p]:
            warn(f"no {p}_depth<N>_clean_ltr.tsv in {indir}; skipping this genome")
    kept = [(p, g) for p, g in zip(prefixes, genomes) if tables_by_prefix[p]]
    if not kept:
        raise SystemExit(f"reboundary: no _clean_ depth tables for any prefix in {indir}")
    prefixes, genomes = [p for p, _ in kept], [g for _, g in kept]
    tables = {p: tables_by_prefix[p] for p in prefixes}
    members = [m for p in prefixes for m in members_from(p, tables[p])]
    by_uid = {m.uid: m for m in members}
    paths = dict(zip(prefixes, genomes))
    g = Genomes(paths)
    g.index()              # once, here: workers left to it race to build one .fai
    trusted = tp.trusted(members, g)
    templates = tp.cap_per_family(trusted, s.max_templates_per_family)
    targets = tp.targets(members)
    log(f"start: {len(members)} elements in {len(prefixes)} genome(s); {len(templates)} "
        f"templates (exact TSD), {len(targets)} targets (no TSD)")
    for p in prefixes:
        vlog(f"  {p}: {sum(m.prefix == p for m in members)} elements, "
             f"{sum(t.prefix == p for t in trusted)} trusted templates, "
             f"{sum(t.prefix == p for t in templates)} templates (after the per-family cap), "
             f"{sum(m.prefix == p for m in targets)} targets")

    rows: Dict[str, Dict[str, str]] = {}
    props: Dict[str, Proposal] = {}
    verdicts: List[Verdict] = []
    partner_of: Dict[str, str] = {}              # survivor uid -> the uid it absorbs
    if not templates:
        warn("no element carries an exact TSD at its called ends, so there are no templates; "
             "nothing was moved")
    elif targets:
        work = os.path.join(indir, prefixes[0] + "_reboundary.tmp")
        try:
            near = tp.nearest(targets, templates, g, k=s.place.n_templates, workdir=work,
                              blastn=blastn, threads=s.threads, task=s.blast_task)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        jobs = [(m, near[m.uid], s.place) for m in targets if m.uid in near]
        log(f"templates found for {len(jobs)} of {len(targets)} targets "
            f"({time.time() - t0:.0f} s)")
        per_target = Counter(len(h) for h in near.values())
        vlog(f"  templates per target: {dict(sorted(per_target.items()))}")
        with Pool(s.threads, initializer=_init, initargs=(paths, tools_dir)) as pool:
            for u, prop, tsd_new, tsd_null, obst in pool.imap_unordered(place_job, jobs,
                                                                        chunksize=8):
                if prop is not None:
                    props[u] = prop
                    rows[u] = sidecar_row(prop, tsd_new, tsd_null, obst)
            order = sorted((p for p in props.values() if p.gate_ok),
                           key=lambda p: genome_order(p.member))
            log(f"placement: {len(props)} candidate(s), {len(order)} with a supported move")
            if verbose:
                sides = Counter(r for p in props.values() for r in (p.reason_left, p.reason_right))
                sources = Counter(src for p in order
                                  for src, moved in ((p.left_src, p.left != p.member.start),
                                                     (p.right_src, p.right != p.member.end))
                                  if moved)
                vlog(f"  side outcomes: {dict(sides)}; end sources of supported moves: "
                     f"{dict(sources)} ({time.time() - t0:.0f} s)")
            index = SpanIndex(members)
            clear: List[Proposal] = []
            absorbed = set()
            for p in order:
                u = p.member.uid
                if u in absorbed:
                    continue
                part = merge_partner(index, by_uid, p.member, p.left, p.right, s.merge_bp)
                pu = f"{p.member.prefix}\t{part}" if part else None
                if pu is not None and (pu in absorbed or pu in partner_of):
                    part, pu = None, None       # taken, or itself absorbing another call
                why = conflict(index, p.member, p.left, p.right,
                               ignore=frozenset({part}) if part else frozenset())
                if why:
                    rows[u]["reason"] = why
                    continue
                if pu is not None:
                    partner_of[u] = pu
                    absorbed.add(pu)
                    clear = [c for c in clear if c.member.uid != pu]
                clear.append(p)
            for u in mutual_conflicts([(p.member, p.left, p.right) for p in clear]):
                rows[u]["reason"] = "overlaps_element"
            final = [p for p in clear if rows[p.member.uid]["reason"] == "."]
            for u, pu in list(partner_of.items()):
                if rows[u]["reason"] != ".":
                    _unmerge(rows, partner_of.pop(u))
            cigars = {}
            for prefix in prefixes:
                for t in tables[prefix]:
                    for row in t.rows:
                        cigars[(prefix, row[0])] = t.cols.get(row, "cigar")
            arb = []
            for prefix in prefixes:
                mine = [p for p in final if p.member.prefix == prefix]
                recs = fetch_records([t.fasta for t in tables[prefix]],
                                     {p.member.name for p in mine})
                for p in mine:
                    rec = recs.get(p.member.name)
                    if rec is None:
                        rows[p.member.uid]["reason"] = "record_missing"
                        _unmerge(rows, partner_of.pop(p.member.uid, None))
                    else:
                        arb.append((p, rec, cigars.get((prefix, p.member.name), "."),
                                    s.mutation_rate))
            if verbose:
                lost = Counter(rows[p.member.uid]["reason"] for p in order
                               if rows[p.member.uid]["reason"] != ".")
                vlog(f"  conflicts: {dict(lost)}; split calls to merge: {len(partner_of)}")
            log(f"re-measuring {len(arb)} element(s) with Kmer2LTR")
            verdicts = pool.map(arbitrate, arb, chunksize=8) if arb else []
            vlog(f"  Kmer2LTR status: {dict(Counter(v.k2l_status for v in verdicts))}; "
                 f"verdicts: {dict(Counter(v.status for v in verdicts))} "
                 f"({time.time() - t0:.0f} s)")

    accepted: Dict[str, Dict[str, Accepted]] = defaultdict(dict)
    retired: Dict[str, set] = defaultdict(set)
    fetch = (lambda prefix, chrom, a, b: g.fetch(prefix, chrom, a, b)[0])
    nested = None             # host -> nested elements, built on the first trim that needs it
    for v in verdicts:
        row = rows[v.uid]
        row["k2l_status"] = v.k2l_status
        pu = partner_of.get(v.uid)
        if v.accepted is None:
            row["reason"] = v.status
            _unmerge(rows, pu)
            continue
        a = v.accepted
        m = a.member
        trimmed = ([(m.start, a.left - 1)] if a.left > m.start else []) + \
                  ([(a.right + 1, m.end)] if a.right < m.end else [])
        if trimmed and hosts_of(m.nest_status):
            if nested is None:
                nested = nested_index(members)
            a = replace(a, restore=trim_restore(m, trimmed, by_uid, IUPAC_DEPTH_SEQ, fetch,
                                                nested=nested))
        plus = props[v.uid].orient == "+"
        u5, u3 = v.unpaired if plus else v.unpaired[::-1]
        e5, e3 = ((m.start - a.left, a.right - m.end) if plus
                  else (a.right - m.end, m.start - a.left))
        row.update(decision="moved", reason=".", new_seq_id=a.new_name, tsd_new=v.tsd,
                   k2p_new=v.k2p, ext5=str(e5), ext3=str(e3), unpaired5=str(u5),
                   unpaired3=str(u3))
        accepted[m.prefix][m.key] = a
        if pu is not None:
            x = by_uid[pu]
            retired[x.prefix].add(x.key)
            rows[pu] = _partner_row(x, a.new_name)

    commit = Commit()
    try:
        for prefix in prefixes:
            if accepted[prefix] or retired[prefix]:
                rewrite(tables[prefix], accepted[prefix], IUPAC_DEPTH_SEQ, commit,
                        retired=frozenset(retired[prefix]))
            mine = sorted((r for u, r in rows.items() if u.split("\t", 1)[0] == prefix),
                          key=_row_order)
            with commit.open(os.path.join(indir, prefix + SIDECAR_SUFFIX)) as fh:
                fh.write(sidecar_text(mine))
        commit.commit()
    except BaseException:
        commit.abort()
        raise
    if dump:
        _dump(dump, props, rows)
    counts = Counter(r["decision"] if r["decision"] in ("moved", "merged") else r["reason"]
                     for r in rows.values())
    log(f"done in {time.time() - t0:.0f} s: {len(rows)} candidates, "
        f"{counts.get('moved', 0)} moved, {counts.get('merged', 0)} merged; "
        + ", ".join(f"{k} {n}" for k, n in sorted(counts.items())
                    if k not in ("moved", "merged")))
    return dict(counts)


def resolve_mutation_rate(indir: str, prefixes: Sequence[str], given: Optional[float]) -> float:
    if given is not None:
        return given
    for p in prefixes:
        try:
            with open(os.path.join(indir, p + ".detect.json")) as fh:
                return float(json.load(fh)["settings"]["mutation_rate"])
        except (OSError, KeyError, ValueError, TypeError):
            continue
    warn("no --mutation-rate and no <prefix>.detect.json to read one from; using 3e-8")
    return 3e-8


BACKUP_SUFFIX = "_pre_reboundary"


def _clean_files(directory: str, prefix: str) -> List[str]:
    return sorted(glob.glob(os.path.join(directory, f"{prefix}_depth*_clean_ltr.tsv"))
                  + glob.glob(os.path.join(directory, f"{prefix}_depth*_clean_ltr.fa")))


def _derived(indir: str, prefix: str) -> List[str]:
    names = (f"{prefix}_all_depth_LTR_cleaned.gff3", f"{prefix}_all_depth_protein_LTR_cleaned.gff3",
             f"{prefix}_plots", prefix + SIDECAR_SUFFIX)
    return [os.path.join(indir, n) for n in names if os.path.lexists(os.path.join(indir, n))]


def _remove(path: str) -> None:
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif os.path.lexists(path):
        os.remove(path)


def _check_clean_set(directory: str, prefix: str) -> List[str]:
    """`directory`'s clean tables and FASTAs, refusing an incomplete set.

    Called before anything is deleted, so that a backup which lost files (an
    interrupted restore, a failed rmdir on NFS, a pruned directory) stops the
    run while the working copies are still the ones on disk. An incomplete
    backup that was merely non-empty used to go unnoticed and take a whole depth
    out of the annotation.
    """
    tables = sorted(glob.glob(os.path.join(directory, f"{prefix}_depth*_clean_ltr.tsv")))
    if not tables:
        raise SystemExit(f"reboundary: {directory} holds no {prefix}_depth<N>_clean_ltr.tsv "
                         f"tables; nothing has been removed")
    orphans = [t for t in tables if not os.path.isfile(t[:-len(".tsv")] + ".fa")]
    if orphans:
        names = ", ".join(os.path.basename(p) for p in orphans)
        raise SystemExit(f"reboundary: {directory} is incomplete: no .fa beside {names}; "
                         f"nothing has been removed")
    return _clean_files(directory, prefix)


def prepare_posthoc(indir: str, prefix: str) -> str:
    """Move the originals into <prefix>_pre_reboundary/ once; put fresh copies back to work on.

    If the backup already exists this run starts again from it, so post-hoc
    runs with different settings never stack on each other. The backup is
    checked for a complete set of clean tables before any working copy is
    removed, so a run that cannot be restarted from it removes nothing.
    """
    bdir = os.path.join(indir, prefix + BACKUP_SUFFIX)
    if os.path.isdir(bdir + ".partial"):
        raise SystemExit(f"reboundary: {bdir}.partial exists: an earlier backup was interrupted. "
                         f"Move its files back into {indir}, delete it, then re-run.")
    if not os.path.isdir(bdir):
        if os.path.lexists(os.path.join(indir, prefix + SIDECAR_SUFFIX)):
            # A pipeline run re-boundaries in place and keeps no backup: these tables
            # are already moved calls, and a second pass is not the first one again.
            raise SystemExit(
                f"reboundary: {prefix} was already re-boundaried in {indir} "
                f"({prefix}{SIDECAR_SUFFIX} exists, {os.path.basename(bdir)}/ does not); a "
                f"post-hoc run would re-boundary the moved calls a second time. Start from "
                f"a run made without --reboundary.")
        _check_clean_set(indir, prefix)
        tmp = bdir + ".partial"
        os.makedirs(tmp)
        for path in _clean_files(indir, prefix) + _derived(indir, prefix):
            shutil.move(path, os.path.join(tmp, os.path.basename(path)))
        os.rename(tmp, bdir)
        originals = _clean_files(bdir, prefix)
        log(f"{prefix}: originals moved to {os.path.basename(bdir)}/")
    else:
        originals = _check_clean_set(bdir, prefix)     # before a single working copy goes
        for path in _clean_files(indir, prefix) + _derived(indir, prefix):
            _remove(path)
        log(f"{prefix}: starting again from {os.path.basename(bdir)}/")
    for path in originals:
        shutil.copy2(path, os.path.join(indir, os.path.basename(path)))
    return bdir


def _fsync(path: str) -> None:
    """One copied file's bytes on disk before it is renamed into place, as `Commit` does."""
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _put_back(src: str, dst: str) -> None:
    """Copy one backed-up entry over its place in the run, leaving the backup untouched.

    Durable the way the commit path is: every byte is fsynced before the rename
    and the directory after it, so a power loss cannot leave a truncated
    original behind once the backup has been deleted.
    """
    tmp = dst + ".partial"
    here = os.path.dirname(os.path.abspath(dst))
    if os.path.isdir(src) and not os.path.islink(src):
        _remove(tmp)
        shutil.copytree(src, tmp, symlinks=True)
        for root, _dirs, names in os.walk(tmp):
            for name in names:
                path = os.path.join(root, name)
                if not os.path.islink(path):
                    _fsync(path)
        _remove(dst)
        os.rename(tmp, dst)
        sync_dir(here)
        return
    shutil.copy2(src, tmp)
    _fsync(tmp)
    os.replace(tmp, dst)
    sync_dir(here)


def restore(indir: str, prefix: str) -> None:
    """Put <prefix>_pre_reboundary/ back in place, copying before it deletes anything.

    The backup stays whole until every one of its entries is in place, so an
    interrupted restore loses nothing and re-running it finishes the job. A
    moving restore could split the originals between the two directories, and
    the next attempt would then delete the half it had already put back.

    A derived output the backup has no counterpart for -- a GFF3 or a plots
    directory `regenerate()` wrote where the run had none -- is kept, because
    nothing here can rebuild it, and named in a warning, because it describes
    the extended coordinates the restored tables no longer carry.
    """
    bdir = os.path.join(indir, prefix + BACKUP_SUFFIX)
    if not os.path.isdir(bdir):
        raise SystemExit(f"reboundary: nothing to restore: no {bdir}")
    held = sorted(os.listdir(bdir))
    for name in held:
        _put_back(os.path.join(bdir, name), os.path.join(indir, name))
    if prefix + SIDECAR_SUFFIX not in held:
        _remove(os.path.join(indir, prefix + SIDECAR_SUFFIX))   # this run's own output
    orphans = [os.path.basename(p) for p in _derived(indir, prefix)
               if os.path.basename(p) not in held]
    shutil.rmtree(bdir)
    log(f"{prefix}: originals restored from {os.path.basename(bdir)}/")
    if orphans:
        warn(f"{prefix}: left in place and now disagreeing with the restored tables: "
             f"{', '.join(orphans)}. Written from the re-boundaried coordinates, with no "
             f"earlier copy in the backup to put back; rebuild from the restored tables "
             f"with ltrquest.gff3 (and the plots).")


def regenerate(indir: str, prefix: str, genome: str, plots: bool) -> None:
    """Rewrite the GFF3s (through the key map) and, unless told not to, the plots."""
    from . import gff3
    cons = sorted(glob.glob(os.path.join(indir, "*_all_ltr.consensus_id*_cluster.tsv")))
    if len(cons) > 1:
        raise SystemExit(f"reboundary: {len(cons)} consensus cluster tables in {indir}; "
                         f"expected at most one")
    fam = {}
    if cons:
        fam = {"consensus_cluster": cons[0],
               "family_prefix": os.path.basename(cons[0]).split("_all_ltr.consensus_id")[0]}
    rec = os.path.join(indir, prefix + "_strand_recovery.tsv")
    side = os.path.join(indir, prefix + SIDECAR_SUFFIX)
    if gff3.convert(prefix, indir, genome, recovered_strands=rec if os.path.isfile(rec) else None,
                    reboundary_map=side if os.path.isfile(side) else None, **fam) != 0:
        raise SystemExit(f"reboundary: GFF3 regeneration failed for {prefix}")
    if plots:
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts", "plots.sh")
        env = dict(os.environ, LTRQUEST_PYTHON=sys.executable)
        done = subprocess.run(["bash", script, "--prefix", prefix, "--genome", genome,
                               "--indir", indir], env=env)
        if done.returncode != 0:
            warn(f"plotting reported failures for {prefix}; the tables and GFF3 are unaffected")


def _at_least(lo: int):
    def parse(value: str) -> int:
        x = int(value)
        if x < lo:
            raise argparse.ArgumentTypeError(f"must be >= {lo}, got {value}")
        return x
    return parse


def _fraction(value: str) -> float:
    x = float(value)
    if not 0.0 < x <= 1.0:
        raise argparse.ArgumentTypeError(f"must be in (0, 1], got {value}")
    return x


def _positive(value: str) -> float:
    x = float(value)
    if not x > 0.0:
        raise argparse.ArgumentTypeError(f"must be > 0, got {value}")
    return x


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="ltrquest-reboundary",
        description="Re-bound LTR-RT calls whose ends carry no TSD, using the calls whose ends "
                    "carry an exact one as templates. Rewrites the _clean_ depth tables and "
                    "FASTAs in place and writes <prefix>_reboundary.tsv.")
    ap.add_argument("--indir", default=".", help="directory holding the run (default: .)")
    ap.add_argument("--prefix", nargs="+", required=True,
                    help="genome prefix(es); all are pooled, since templates span genomes")
    ap.add_argument("--genome", nargs="+",
                    help="original (unmasked) genome FASTA per prefix, in the same order")
    ap.add_argument("--templates", type=_at_least(1), default=5,
                    help="nearest templates placed per element (default 5)")
    ap.add_argument("--min-support", type=_at_least(1), default=2,
                    help="agreeing template placements a moved end needs (default 2)")
    ap.add_argument("--min-ext", type=_at_least(1), default=1,
                    help="the smallest outward move, bp (default 1)")
    ap.add_argument("--max-trim", type=_at_least(0), default=10,
                    help="the largest inward move, bp; 0 never trims (default 10)")
    ap.add_argument("--min-identity", type=_fraction, default=0.8,
                    help="identity over a template's outer 30 bp a placement needs (default 0.8)")
    ap.add_argument("--anchor-len", type=_at_least(1), default=30,
                    help="bp of a template's end searched for past a large indel (default 30)")
    ap.add_argument("--no-anchor", action="store_true", help="do not search past large indels")
    ap.add_argument("--max-indel", type=_at_least(0), default=5000,
                    help="how far past the call the anchor looks, bp (default 5000)")
    ap.add_argument("--mutation-rate", type=_positive, default=None,
                    help="per site per year, for k2p_time (default: the run's own, else 3e-8)")
    ap.add_argument("--tools-dir", default=os.environ.get("LTRQUEST_TOOLS_DIR", "ltrquest_tools"),
                    help="Kmer2LTR checkout location, cloned there if absent "
                         "(default: $LTRQUEST_TOOLS_DIR or ./ltrquest_tools)")
    ap.add_argument("--blastn", default="blastn",
                    help="blastn executable; makeblastdb is taken from beside it or PATH "
                         "(default: blastn on PATH)")
    ap.add_argument("-t", "--threads", type=_at_least(1), default=8,
                    help="worker processes (default 8)")
    ap.add_argument("--dump-proposals", default=None,
                    help="also write every proposal to this TSV (benchmarks)")
    ap.add_argument("-v", "--verbose", action="store_true", help="more progress")
    ap.add_argument("--posthoc", action="store_true",
                    help="update a finished run in place: originals go to "
                         "<prefix>_pre_reboundary/ (once; later runs start from it), then the "
                         "GFF3s and plots are rewritten")
    ap.add_argument("--no-plots", action="store_true", help="with --posthoc: skip the plots")
    ap.add_argument("--restore", action="store_true",
                    help="put <prefix>_pre_reboundary/ back in place and stop")
    return ap


def settings_from(args) -> Settings:
    return Settings(
        mutation_rate=resolve_mutation_rate(args.indir, args.prefix, args.mutation_rate),
        threads=args.threads, blastn=args.blastn,
        place=Params(min_identity=args.min_identity, anchor=not args.no_anchor,
                     anchor_len=args.anchor_len, max_indel=args.max_indel,
                     n_templates=args.templates, min_support=args.min_support,
                     min_ext=args.min_ext, max_trim=args.max_trim))


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    if args.min_support > args.templates:
        ap.error(f"--min-support {args.min_support} exceeds --templates {args.templates}: "
                 f"no end could ever move")
    check_prefixes(args.prefix)
    if args.restore:
        for p in args.prefix:
            restore(args.indir, p)
        return 0
    check_indir(args.indir)        # before --posthoc moves a single original
    if not args.genome:
        raise SystemExit("reboundary: --genome is required")
    if len(args.genome) != len(args.prefix):
        raise SystemExit("reboundary: --prefix and --genome must pair up one to one")
    s = settings_from(args)
    check_blastn(s.blastn)         # both before --posthoc backs anything up
    try:
        k2l.api(args.tools_dir)
    except k2l.IncompatibleKmer2LTR as exc:
        print(f"reboundary: {exc}", file=sys.stderr, flush=True)
        return EXIT_INCOMPATIBLE_KMER2LTR
    if args.posthoc:
        for p in args.prefix:
            prepare_posthoc(args.indir, p)
    run(args.indir, args.prefix, args.genome, s, args.tools_dir, args.dump_proposals,
        args.verbose)
    if args.posthoc:
        for p, gpath in zip(args.prefix, args.genome):
            regenerate(args.indir, p, gpath, plots=not args.no_plots)
    return 0


if __name__ == "__main__":
    sys.exit(main())

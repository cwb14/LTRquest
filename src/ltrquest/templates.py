"""Trusted templates for ltrquest.reboundary: the calls that vouch for their own ends.

An element whose called ends carry an exact target-site duplication (Kmer2LTR's
`tsd` found at offset `0,0`) has both outer ends where the insertion put them:
the duplication is read from the host flanks at that element's own locus. Such
elements are the templates; every element with no duplication at its called ends
is a target. A template's TSD is a different insertion event from any target's,
so a target's own TSD stays an independent check on where it is re-bounded.

Templates are found per target by one BLASTN of every target's LTRs against every
template's genomic-left LTR, pooled over all genomes and families: no family
label, copy number or consensus is involved. Each template's relation to a target
(`+` same genomic frame, `-` opposite) is the strand of its best hit.
"""

from __future__ import annotations

import math
import os
import random
import re
import shutil
import subprocess
import zlib
from collections import defaultdict
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

from .ltr_model import Genomes, Member, rc

NEST_OUTER = re.compile(r"^nest-outer:(.+):(\d+)-(\d+)$")
EVALUE = "1e-10"


def _found(tsd: str) -> bool:
    return tsd not in ("", ".", "NA")


def _nested_spans(m: Member) -> List[Tuple[str, int, int]]:
    """(chrom, start, end) of every element the table says is nested in `m`."""
    out = []
    for tok in (m.nest_status or ".").split(";"):
        hit = NEST_OUTER.match(tok)
        if hit:
            out.append((hit.group(1), int(hit.group(2)), int(hit.group(3))))
    return out


def ltr_seqs(m: Member, g: Genomes) -> Tuple[str, str]:
    """(genomic-left LTR, genomic-right LTR), forward frame, as called."""
    return (g.fetch(m.prefix, m.chrom, m.start, m.l1)[0],
            g.fetch(m.prefix, m.chrom, m.r0, m.end)[0])


def trusted(members: Iterable[Member], g: Genomes) -> List[Member]:
    """The templates: calls with an exact TSD at their called ends.

    Skipped even so: an LTR holding an assembly gap, or overlapping an element the
    table says is nested in the call -- its sequence would carry the other TE into
    the search and the placement.
    """
    out: List[Member] = []
    for m in members:
        if not _found(m.tsd) or m.tsd_offset != "0,0":
            continue
        ltrs = ((m.start, m.l1), (m.r0, m.end))
        if any(c == m.chrom and a <= hi and b >= lo
               for c, a, b in _nested_spans(m) for lo, hi in ltrs):
            continue
        left, right = ltr_seqs(m, g)
        if "N" in left or "N" in right:
            continue
        out.append(m)
    return out


def cap_per_family(templates: Sequence[Member], cap: int) -> List[Member]:
    """At most `cap` templates per family: its youngest half by K2P, then a draw seeded
    by the family name (so the choice does not depend on table order). Only `k`
    templates are placed per target, but every template in the database is aligned
    against every query that hits its family, so a 20,000-copy family would dominate
    the search time. Unlabelled templates and families under the cap are kept
    whole; `cap <= 0` keeps everything. Input order is kept for what survives.
    """
    if cap <= 0:
        return list(templates)
    by_family: Dict[str, List[Member]] = defaultdict(list)
    for t in templates:
        if t.family not in ("", ".", "NA"):
            by_family[t.family].append(t)
    keep = set()
    for fam, ms in by_family.items():
        if len(ms) <= cap:
            keep.update(m.uid for m in ms)
            continue
        age = sorted(ms, key=lambda m: (math.inf if m.k2p is None else m.k2p, m.uid))
        young = age[:cap // 2]
        rest = sorted(age[cap // 2:], key=lambda m: m.uid)
        drawn = random.Random(zlib.crc32(fam.encode())).sample(rest, cap - len(young))
        keep.update(m.uid for m in young + drawn)
    return [t for t in templates if t.family in ("", ".", "NA") or t.uid in keep]


def targets(members: Iterable[Member]) -> List[Member]:
    """Calls with no TSD at their called ends: the ones re-boundarying may move."""
    return [m for m in members if not _found(m.tsd)]


def models_for(t: Member, rel: str, g: Genomes) -> Tuple[str, str]:
    """(left_model, right_model) in the target's forward frame.

    Each side is placed with the template LTR whose OUTER end is its outermost base
    there -- the end the template's TSD verifies. Same frame: the template's own
    left and right LTRs. Opposite frame: its right LTR, reverse-complemented, on the
    target's left, and its left LTR on the target's right. (One LTR placed on both
    sides would let a template's unverified inner end set the target's opposite
    outer end.)
    """
    left, right = ltr_seqs(t, g)
    return (left, right) if rel == "+" else (rc(right), rc(left))


def _implied(target: Member, template: Member) -> Optional[str]:
    """The relation two known strands imply, or None when either is unknown."""
    if target.stranded and template.stranded:
        return "+" if target.strand == template.strand else "-"
    return None


def rank(hits: Iterable[Tuple[Member, Member, float, str]], k: int,
         exclude: FrozenSet[str] = frozenset()) -> Dict[str, List[Tuple[Member, str]]]:
    """target uid -> up to `k` (template, relation), best first.

    `hits` are (target, template, bitscore, relation) rows. A template counts once,
    at its best row; an exact tie between relations goes to the one the two
    strands imply, else `+`. A template whose relation contradicts two known
    strands is dropped. Ordered by (-bitscore, template uid), so ties are
    deterministic.
    """
    best: Dict[Tuple[str, str], Tuple[float, str, Member, Member]] = {}
    for target, tpl, bits, rel in hits:
        if tpl.uid in exclude or tpl.uid == target.uid:
            continue
        key = (target.uid, tpl.uid)
        prev = best.get(key)
        if prev is None or bits > prev[0]:
            best[key] = (bits, rel, target, tpl)
        elif bits == prev[0] and rel != prev[1]:
            best[key] = (bits, _implied(target, tpl) or "+", target, tpl)
    per: Dict[str, List[Tuple[float, str, Member, str]]] = defaultdict(list)
    for (tu, pu), (bits, rel, target, tpl) in best.items():
        implied = _implied(target, tpl)
        if implied is not None and implied != rel:
            continue
        per[tu].append((-bits, pu, tpl, rel))
    return {tu: [(tpl, rel) for _, _, tpl, rel in sorted(v, key=lambda x: (x[0], x[1]))[:k]]
            for tu, v in per.items()}


def _run(argv: Sequence[str]) -> None:
    done = subprocess.run(list(argv), capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError(f"{os.path.basename(argv[0])} failed ({done.returncode}): "
                           f"{done.stderr.strip()[-2000:]}")


def nearest(targets: Sequence[Member], templates: Sequence[Member], g: Genomes, k: int,
            workdir: str, blastn: str = "blastn", threads: int = 1, task: str = "blastn",
            exclude: FrozenSet[str] = frozenset(),
            makeblastdb: Optional[str] = None) -> Dict[str, List[Tuple[Member, str]]]:
    """target uid -> up to `k` (template, relation), best first (see `rank`).

    One BLASTN of both LTRs of every target (genomic forward) against every
    template's genomic-left LTR, both strands, in `workdir`. Sequences are named by
    index: a Member's uid carries a TAB and element names repeat across genomes.
    """
    if not targets or not templates:
        return {}
    os.makedirs(workdir, exist_ok=True)
    if makeblastdb is None:        # beside blastn, else on PATH (what check_blastn accepts)
        beside = os.path.join(os.path.dirname(blastn), "makeblastdb")
        makeblastdb = (beside if os.path.dirname(blastn) and os.path.isfile(beside)
                       else shutil.which("makeblastdb") or "makeblastdb")
    db = os.path.join(workdir, "templates")
    with open(db + ".fa", "w") as fh:
        for i, t in enumerate(templates):
            fh.write(f">t{i}\n{ltr_seqs(t, g)[0]}\n")
    _run([makeblastdb, "-in", db + ".fa", "-dbtype", "nucl", "-out", db])
    query = os.path.join(workdir, "targets.fa")
    with open(query, "w") as fh:
        for i, m in enumerate(targets):
            left, right = ltr_seqs(m, g)
            fh.write(f">q{i}L\n{left}\n>q{i}R\n{right}\n")
    out = os.path.join(workdir, "hits.tsv")
    _run([blastn, "-task", task, "-query", query, "-db", db, "-evalue", EVALUE,
          "-max_hsps", "1", "-max_target_seqs", str(max(50, 10 * k)),
          "-outfmt", "6 qseqid sseqid bitscore sstrand", "-num_threads", str(threads),
          "-out", out])
    # blastn writes each query's rows together, in query order (L then R per target),
    # so a target is ranked as soon as its rows end: only its top k is ever held.
    found: Dict[str, List[Tuple[Member, str]]] = {}
    group: List[Tuple[Member, Member, float, str]] = []
    cur, seen = -1, set()
    with open(out) as fh:
        for line in fh:
            qid, sid, bits, strand = line.rstrip("\n").split("\t")
            ti = int(qid[1:-1])
            if ti != cur:
                if ti in seen:
                    raise RuntimeError(f"{os.path.basename(blastn)} output is not grouped by "
                                       f"query ({qid} after other queries' rows)")
                found.update(rank(group, k, exclude))
                group, cur = [], ti
                seen.add(ti)
            group.append((targets[ti], templates[int(sid[1:])], float(bits),
                          "+" if strand == "plus" else "-"))
    found.update(rank(group, k, exclude))
    return found

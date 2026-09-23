"""Family-guided re-boundarying of truncated LTR-RT calls (opt-in stage).

LTRharvest and LTR_FINDER extend an LTR pair outward from a seed and stop at the
first obstacle between the element's two LTRs -- an indel, a patch of mutations.
Kmer2LTR can only trim. So an element whose LTR carries an obstacle near its end
is called short. This stage builds each family's LTR model from its full-length
copies (pooled over every genome), places it on every member, and proposes
outward-only ends. Candidates pass family QC and conflict checks, then Kmer2LTR
re-scores each widened record with the family's evidence as `tsd_credit`: it
settles the final ends and every Kmer2LTR column. Only the `_clean_` tables and
FASTAs are rewritten; `<prefix>_reboundary.tsv` records every candidate and why
it was or was not extended, and doubles as the old->new key map for the GFF3.

Usage:
  ltrquest-reboundary --indir RUN --prefix P1 [P2 ...] --genome G1 [G2 ...] [options]
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import pickle
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field, replace
from multiprocessing import Pool
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

from . import kmer2ltr as k2l
from .kmer2ltr import COLUMNS
from .ltr_model import (
    Genomes,
    Member,
    Model,
    consensus_model,
    ratio_ok,
    select_references,
    subfamily_models,
    tsd_enrichment,
)
from .ltr_place import Params, Proposal, nearest_templates, obstacle, propose, tsd_at
from .reboundary_io import (
    SIDECAR_SUFFIX,
    Accepted,
    Commit,
    SpanIndex,
    conflict,
    fetch_records,
    forward,
    leftover_message,
    leftover_staging,
    load_clean_tables,
    members_from,
    mutual_conflicts,
    rewrite,
    sanitize,
    sidecar_text,
    stored,
    sync_dir,
)
from .reconcile import IUPAC_DEPTH_SEQ

METHODS = ("consensus", "subfamily", "nearest")
_I = {c: i for i, c in enumerate(COLUMNS)}
# Bump when a pickled FamilyResult changes shape, so `--cache` files written by an
# older build are recomputed instead of being loaded into the new code.
# 2: Model carries `refs` and no longer pickles its k-mer index.
CACHE_FORMAT = 2


def log(msg: str) -> None:
    print(f"[reboundary] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[reboundary] WARNING: {msg}", file=sys.stderr, flush=True)


@dataclass(frozen=True)
class Settings:
    method: str = "nearest"
    references: str = "modal"
    min_copies: int = 10
    credit: str = "5000"
    max_ratio: float = 1.15
    qc_min_n: int = 5
    untested: str = "accept"
    subfamily_jaccard: float = 0.5
    mutation_rate: float = 3e-8
    mafft: str = "mafft"
    threads: int = 8
    place: Params = field(default_factory=Params)

    def family_key(self) -> str:
        """Everything the family phase depends on, so a cached phase is reused only when valid."""
        # min_copies is left out on purpose: it only chooses which families run, and a
        # cached phase covering more families is reused for fewer (see _family_phase).
        keep = ("method", "references", "subfamily_jaccard", "max_ratio", "place")
        d = {k: v for k, v in asdict(self).items() if k in keep}
        d["cache_format"] = CACHE_FORMAT
        return hashlib.sha1(json.dumps(d, sort_keys=True).encode()).hexdigest()


_G: Optional[Genomes] = None
_API = None


def _init(paths: Dict[str, str], tools_dir: str) -> None:
    """Worker set-up: genome handles and the Kmer2LTR API, once per process."""
    global _G, _API
    _G = Genomes(paths)
    _API = k2l.api(tools_dir)


def build_models(family: str, members: Sequence[Member], g: Genomes, s: Settings,
                 exclude: FrozenSet[str] = frozenset()) -> Tuple[List[Model], Optional[float], str]:
    """The family's model(s), the modal called LTR length, and what stopped a model.

    Leave-one-out does not depend on this call: every model records the copies whose
    sequence went into it (`Model.refs`), and `ltr_place.candidate_models` refuses a
    model for the element that built it. `exclude` is the stronger, caller-driven
    form -- the copies named are kept out of the reference pool altogether, which the
    benchmark uses to keep a planted element out of its family's model. It can drop a
    family below MIN_REFS; that returns 'too_few_references' and no model, which
    `family_job` and the sidecar already report.
    """
    if s.method not in METHODS:
        raise ValueError(f"unknown method {s.method!r}")
    many = s.method == "subfamily"
    refs, modal = select_references(members, s.references, family, exclude,
                                    n_young=150 if many else 40, n_random=150 if many else 40)
    if not refs:
        return [], modal, "too_few_references"
    if s.method == "subfamily":
        models = subfamily_models(family, refs, g, modal, s.subfamily_jaccard, s.mafft)
        return models, modal, ("ok" if models else "no_ltr_span")
    model = consensus_model(family, refs, g, modal, s.mafft)
    if model is None:
        return [], modal, "no_ltr_span"
    if s.method == "nearest":
        return nearest_templates(model, refs, g), modal, "ok"
    return [model], modal, "ok"


def place_params(s: Settings) -> Params:
    """`nearest` takes the median of its top templates' placements (spec 5.1)."""
    return replace(s.place, combine="median") if s.method == "nearest" else s.place


@dataclass
class FamilyResult:
    family: str
    n_members: int
    models: List[Model]
    modal: Optional[float]
    status: str
    proposals: List[Proposal] = field(default_factory=list)
    tsd_new: Dict[str, str] = field(default_factory=dict)
    tsd_null: Dict[str, str] = field(default_factory=dict)
    obstacles: Dict[str, Tuple[str, str]] = field(default_factory=dict)


def family_job(args) -> FamilyResult:
    family, members, s = args
    models, modal, status = build_models(family, members, _G, s)
    fr = FamilyResult(family, len(members), models, modal, status)
    if not models:
        return fr
    kept = [mo for mo in models if ratio_ok(mo, s.max_ratio)]
    if not kept:
        fr.status = "ratio_failed"   # proposals still computed, so the sidecar says what failed
    for m in members:
        p = propose(m, kept or models, _G, place_params(s))
        if p is None:
            continue
        fr.proposals.append(p)
        if p.gate_ok:
            fr.tsd_new[m.uid] = tsd_at(_API, _G, m, p.left, p.right)
            fr.tsd_null[m.uid] = tsd_at(_API, _G, m, p.left, p.right, shift=1000)
            fr.obstacles[m.uid] = obstacle(m, p, _G)
    return fr


def _found(tsd: str) -> bool:
    return tsd not in ("", ".", "NA")


def qc(results: Sequence[FamilyResult], s: Settings) -> Tuple[Dict[str, str], float]:
    """Per-family QC status and the pooled displaced-flank null rate it was tested against.

    Evidence: TSDs at the proposed ends of candidates whose called ends had none
    (TSDs are a credibility signal the finders did not select on; TG..CA is not
    used because they did).
    """
    null_hits = null_n = 0
    for fr in results:
        for p in fr.proposals:
            u = p.member.uid
            if p.gate_ok and not p.member.has_tsd and u in fr.tsd_null:
                null_n += 1
                null_hits += _found(fr.tsd_null[u])
    p0 = null_hits / null_n if null_n else 0.0
    status: Dict[str, str] = {}
    for fr in results:
        if fr.status != "ok":
            status[fr.family] = fr.status
            continue
        ev = [fr.tsd_new[p.member.uid] for p in fr.proposals
              if p.gate_ok and not p.member.has_tsd and p.member.uid in fr.tsd_new]
        status[fr.family] = tsd_enrichment(sum(map(_found, ev)), len(ev), p0, s.qc_min_n)
    return status, p0


def family_accepts(status: str, s: Settings) -> bool:
    return status == "pass" or (status == "untested" and s.untested == "accept")


def credit_for(p: Proposal, s: Settings) -> float:
    """Bits of family evidence handed to Kmer2LTR that the proposed termini are the ends."""
    if s.credit == "model":
        ids = (([p.id_left] if p.left < p.member.start else [])
               + ([p.id_right] if p.right > p.member.end else []))
        return 2.0 * s.place.anchor_len * min(ids) if ids else 0.0
    return float(s.credit)


@dataclass(frozen=True)
class Verdict:
    uid: str
    status: str
    accepted: Optional[Accepted]
    k2l_status: str
    tsd: str
    k2p: str
    credit: float


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
    """Detection's Step 8a gate (detect.filter_kmer2ltr_in_place) on the settled pair."""
    if l5 is None or l3 is None or aln is None:
        return False
    hi = max(l5, l3, aln)
    return (hi > 0 and min(l5, l3) >= 100 and aln >= 90
            and min(l5, l3, aln) / hi >= 0.65 and right - left >= 300)


def arbitrate(args) -> Verdict:
    """Kmer2LTR re-scores the widened record; it settles the final ends and every column."""
    p, record, credit, mu, min_ext = args
    m = p.member
    left_seg = _G.fetch(m.prefix, m.chrom, p.left, m.start - 1)[0] if p.left < m.start else ""
    right_seg = _G.fetch(m.prefix, m.chrom, m.end + 1, p.right)[0] if p.right > m.end else ""
    fwd = left_seg + forward(record, m.orientation) + right_seg
    if len(fwd) != p.right - p.left + 1:
        return Verdict(m.uid, "record_mismatch", None, "NA", "NA", "NA", credit)
    suffix = "#" + m.name.split("#", 1)[1] if "#" in m.name else ""
    seq = sanitize(fwd)
    ctx = _API.orient(seq, window(m, p.left, p.right))
    res = _API.classify(f"{m.chrom}:{p.left}-{p.right}{suffix}", seq, period_rule="outermost",
                        mutation_rate=mu, tsd_credit=credit)
    if ctx is not None:
        res = _API.annotate(res, seq, ctx, _API.Options())
    if res.status != "pass":
        return Verdict(m.uid, "kmer2ltr_not_pass", None, res.status, "NA", "NA", credit)
    row = _API.format_row(res).split("\t")
    f5, f3 = res.flank5_len, res.flank3_len
    fl, fr = p.left + f5, p.right - f3
    tsd, k2p = row[_I["tsd"]], row[_I["k2p"]]
    if fl > m.start or fr < m.end or (m.start - fl < min_ext and fr - m.end < min_ext):
        return Verdict(m.uid, "kmer2ltr_reverted", None, res.status, tsd, k2p, credit)
    if not length_ok(res.ltr5_len, res.ltr3_len, res.aln_len, fl, fr):
        return Verdict(m.uid, "length_filter", None, res.status, tsd, k2p, credit)
    fields = rebase(row, f5, f"{m.chrom}:{fl}-{fr}{suffix}")
    fields[_I["orientation"]] = m.orientation   # a storage fact: the record stays stored as it was
    acc = Accepted(m, fields[0], fields, stored(fwd[f5:len(fwd) - f3], m.orientation), fl, fr)
    return Verdict(m.uid, "extended", acc, res.status, tsd, k2p, credit)


def sidecar_row(p: Proposal, fr: FamilyResult, s: Settings) -> Dict[str, str]:
    m, u = p.member, p.member.uid
    ob_l, ob_r = fr.obstacles.get(u, (".", "."))
    plus = p.orient == "+"
    return {
        "old_seq_id": m.name, "new_seq_id": ".", "family": m.family, "method": s.method,
        "model_id": p.model_id, "decision": "rejected", "reason": ".",
        "ext5": str(p.ext5 if p.gate_ok else p.raw_ext5),
        "ext3": str(p.ext3 if p.gate_ok else p.raw_ext3),
        "end_source5": p.left_src if plus else p.right_src,
        "end_source3": p.right_src if plus else p.left_src,
        "id_outer5": f"{(p.id_left if plus else p.id_right):.3f}",
        "id_outer3": f"{(p.id_right if plus else p.id_left):.3f}",
        "credit_bits": ".", "k2l_status": ".",
        "tsd_called": m.tsd or ".", "tsd_new": fr.tsd_new.get(u, "."),
        "tsd_null": fr.tsd_null.get(u, "."),
        "k2p_called": "NA" if m.k2p is None else f"{m.k2p:g}", "k2p_new": ".",
        "obstacle5": ob_l if plus else ob_r, "obstacle3": ob_r if plus else ob_l,
    }


def _row_order(row: Dict[str, str]):
    chrom, span = row["old_seq_id"].split("#", 1)[0].rsplit(":", 1)
    return chrom, int(span.split("-")[0])


def _family_phase(pool, eligible, s: Settings, cache: Optional[str], verbose: bool):
    key = s.family_key()
    if cache and os.path.isfile(cache):
        with open(cache, "rb") as fh:
            saved_key, results = pickle.load(fh)
        wanted = {f for f, _ in eligible}
        if saved_key == key and wanted <= {r.family for r in results}:
            log(f"family phase reused from {cache}")
            return [r for r in results if r.family in wanted]
    results = []
    for fr in pool.imap_unordered(family_job, [(f, ms, s) for f, ms in eligible]):
        results.append(fr)
        if verbose:
            lens = ",".join(str(len(mo.seq)) for mo in fr.models) or "-"
            log(f"{fr.family}: {fr.n_members} copies, model {lens} bp, modal {fr.modal}, "
                f"{fr.status}, {len(fr.proposals)} proposal(s)")
    results.sort(key=lambda fr: fr.family)
    if cache:
        with open(cache, "wb") as fh:
            pickle.dump((key, results), fh)
    return results


def _dump(path: str, results, status, rows) -> None:
    with open(path, "w") as fh:
        fh.write("uid\tfamily\tfamily_status\tgate_ok\tcalled_tsd\ttsd_new\ttsd_null\t"
                 "decision\treason\n")
        for fr in results:
            for p in fr.proposals:
                u = p.member.uid
                r = rows[u]
                fh.write("\t".join([u.replace("\t", "|"), fr.family, status[fr.family],
                                    str(int(p.gate_ok)), str(int(p.member.has_tsd)),
                                    fr.tsd_new.get(u, "."), fr.tsd_null.get(u, "."),
                                    r["decision"], r["reason"]]) + "\n")


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


def run(indir: str, prefixes: Sequence[str], genomes: Sequence[str], s: Settings,
        tools_dir: str, dump: Optional[str] = None, cache: Optional[str] = None,
        verbose: bool = False) -> Dict[str, int]:
    t0 = time.time()
    check_prefixes(prefixes)
    check_indir(indir)
    if len(prefixes) != len(genomes):
        raise SystemExit("reboundary: --prefix and --genome must pair up one to one")
    if s.method not in METHODS:
        raise SystemExit(f"reboundary: unknown method {s.method!r}")
    if shutil.which(s.mafft) is None:
        raise SystemExit(f"reboundary: mafft not found ({s.mafft}); it is in environment.yml")
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
    families: Dict[str, List[Member]] = defaultdict(list)
    for m in members:
        if m.family not in ("", ".", "NA"):
            families[m.family].append(m)
    spaces = {f.rsplit("_fam", 1)[0] for f in families if "_fam" in f}
    if len(spaces) > 1:
        raise SystemExit(f"reboundary: the prefixes carry {len(spaces)} family namespaces "
                         f"({', '.join(sorted(spaces))}); pool only genomes whose families "
                         f"were clustered together in one run")
    eligible = sorted(((f, ms) for f, ms in families.items() if len(ms) >= s.min_copies),
                      key=lambda kv: (-len(kv[1]), kv[0]))
    log(f"start: {len(members)} elements in {len(prefixes)} genome(s), {len(families)} "
        f"families, {len(eligible)} with >= {s.min_copies} copies; method {s.method}")

    rows: Dict[str, Dict[str, str]] = {}
    by_uid: Dict[str, Proposal] = {}
    with Pool(s.threads, initializer=_init, initargs=(dict(zip(prefixes, genomes)), tools_dir)) \
            as pool:
        results = _family_phase(pool, eligible, s, cache, verbose)
        status, p0 = qc(results, s)
        log(f"family phase: {sum(len(fr.proposals) for fr in results)} candidates; "
            f"QC {dict(Counter(status.values()))}; displaced-flank null {p0:.3f}")
        survivors: List[Proposal] = []
        for fr in results:
            accept = family_accepts(status[fr.family], s)
            for p in fr.proposals:
                u = p.member.uid
                by_uid[u] = p
                rows[u] = sidecar_row(p, fr, s)
                if not p.gate_ok:
                    rows[u]["reason"] = "gate_identity"
                elif not accept:
                    rows[u]["reason"] = ("family_untested" if status[fr.family] == "untested"
                                         else "family_qc_failed")
                else:
                    survivors.append(p)
        index = SpanIndex(members)
        clear: List[Proposal] = []
        for p in survivors:
            why = conflict(index, p.member, p.left, p.right)
            if why:
                rows[p.member.uid]["reason"] = why
            else:
                clear.append(p)
        lost = mutual_conflicts([(p.member, p.left, p.right) for p in clear])
        for u in lost:
            rows[u]["reason"] = "overlaps_element"
        final = [p for p in clear if p.member.uid not in lost]
        jobs = []
        for prefix in prefixes:
            mine = [p for p in final if p.member.prefix == prefix]
            recs = fetch_records([t.fasta for t in tables[prefix]], {p.member.name for p in mine})
            for p in mine:
                rec = recs.get(p.member.name)
                if rec is None:
                    rows[p.member.uid]["reason"] = "record_missing"
                else:
                    jobs.append((p, rec, credit_for(p, s), s.mutation_rate, s.place.min_ext))
        log(f"arbitration: {len(jobs)} candidate(s) to Kmer2LTR")
        verdicts = pool.map(arbitrate, jobs, chunksize=8) if jobs else []

    accepted: Dict[str, Dict[str, Accepted]] = defaultdict(dict)
    for v in verdicts:
        row = rows[v.uid]
        row["credit_bits"], row["k2l_status"] = f"{v.credit:g}", v.k2l_status
        if v.accepted is None:
            row["reason"] = v.status
            continue
        a = v.accepted
        ext_l, ext_r = a.member.start - a.left, a.right - a.member.end
        plus = by_uid[v.uid].orient == "+"
        row.update(decision="extended", reason=".", new_seq_id=a.new_name, tsd_new=v.tsd,
                   k2p_new=v.k2p, ext5=str(ext_l if plus else ext_r),
                   ext3=str(ext_r if plus else ext_l))
        accepted[a.member.prefix][a.member.key] = a

    commit = Commit()
    try:
        for prefix in prefixes:
            if accepted[prefix]:
                rewrite(tables[prefix], accepted[prefix], IUPAC_DEPTH_SEQ, commit)
            mine = sorted((r for u, r in rows.items() if u.split("\t", 1)[0] == prefix),
                          key=_row_order)
            with commit.open(os.path.join(indir, prefix + SIDECAR_SUFFIX)) as fh:
                fh.write(sidecar_text(mine))
        commit.commit()
    except BaseException:
        commit.abort()
        raise
    if dump:
        _dump(dump, results, status, rows)
    counts = Counter(r["decision"] if r["decision"] == "extended" else r["reason"]
                     for r in rows.values())
    log(f"done in {time.time() - t0:.0f} s: {len(rows)} candidates, "
        f"{counts.get('extended', 0)} extended; "
        + ", ".join(f"{k} {n}" for k, n in sorted(counts.items()) if k != "extended"))
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


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="ltrquest-reboundary",
        description="Extend truncated LTR-RT calls to the ends their family's LTR model "
                    "supports. Rewrites the _clean_ depth tables and FASTAs in place and "
                    "writes <prefix>_reboundary.tsv.")
    ap.add_argument("--indir", default=".", help="directory holding the run (default: .)")
    ap.add_argument("--prefix", nargs="+", required=True,
                    help="genome prefix(es); all are pooled, since families span genomes")
    ap.add_argument("--genome", nargs="+",
                    help="original (unmasked) genome FASTA per prefix, in the same order")
    ap.add_argument("--method", choices=METHODS, default="nearest",
                    help="family model (default: nearest)")
    ap.add_argument("--references", choices=("modal", "tsd"), default="modal",
                    help="copies models are built from: modal-length (default) or TSD-bearing")
    ap.add_argument("--subfamily-jaccard", type=float, default=0.5,
                    help="with --method subfamily: LTR 11-mer Jaccard to join a cluster "
                         "(default 0.5)")
    ap.add_argument("--min-copies", type=int, default=10,
                    help="families with fewer copies are left alone (default 10)")
    ap.add_argument("--min-identity", type=float, default=0.8,
                    help="identity over the outer 30 bp a moved end needs (default 0.8)")
    ap.add_argument("--anchor-len", type=int, default=30,
                    help="bp of the model's end searched for past a large indel (default 30)")
    ap.add_argument("--no-anchor", action="store_true", help="do not search past large indels")
    ap.add_argument("--max-indel", type=int, default=5000,
                    help="how far past the call the anchor looks, bp (default 5000)")
    ap.add_argument("--credit", default="5000",
                    help="bits of family evidence given to Kmer2LTR: a number, or 'model' "
                         "(default 5000)")
    ap.add_argument("--max-ratio", type=float, default=1.15,
                    help="family QC: model length / modal called LTR length ceiling (default 1.15)")
    ap.add_argument("--qc-min-n", type=int, default=5,
                    help="family QC: candidates needed to test TSD enrichment (default 5)")
    ap.add_argument("--untested", choices=("accept", "skip"), default="accept",
                    help="families with too few candidates to test (default: accept)")
    ap.add_argument("--mutation-rate", type=float, default=None,
                    help="per site per year, for k2p_time (default: the run's own, else 3e-8)")
    ap.add_argument("--tools-dir", default=os.environ.get("LTRQUEST_TOOLS_DIR", "ltrquest_tools"),
                    help="Kmer2LTR checkout location, cloned there if absent "
                         "(default: $LTRQUEST_TOOLS_DIR or ./ltrquest_tools)")
    ap.add_argument("--mafft", default="mafft", help="mafft executable (default: mafft on PATH)")
    ap.add_argument("-t", "--threads", type=int, default=8, help="worker processes (default 8)")
    ap.add_argument("--dump-proposals", default=None,
                    help="also write every proposal to this TSV (benchmarks)")
    ap.add_argument("--cache", default=None,
                    help="pickle the family phase here; reuse it when the settings allow")
    ap.add_argument("-v", "--verbose", action="store_true", help="one line per family")
    ap.add_argument("--posthoc", action="store_true",
                    help="update a finished run in place: originals go to "
                         "<prefix>_pre_reboundary/ (once; later runs start from it), then the "
                         "GFF3s and plots are rewritten")
    ap.add_argument("--no-plots", action="store_true", help="with --posthoc: skip the plots")
    ap.add_argument("--restore", action="store_true",
                    help="put <prefix>_pre_reboundary/ back in place and stop")
    return ap


def settings_from(args) -> Settings:
    if args.credit != "model":
        try:
            float(args.credit)
        except ValueError:
            raise SystemExit(f"reboundary: --credit must be a number or 'model' "
                             f"(got {args.credit!r})") from None
    return Settings(
        method=args.method, references=args.references, min_copies=args.min_copies,
        credit=args.credit, max_ratio=args.max_ratio, qc_min_n=args.qc_min_n,
        untested=args.untested, subfamily_jaccard=args.subfamily_jaccard,
        mutation_rate=resolve_mutation_rate(args.indir, args.prefix, args.mutation_rate),
        mafft=args.mafft, threads=args.threads,
        place=Params(min_identity=args.min_identity, anchor=not args.no_anchor,
                     anchor_len=args.anchor_len, max_indel=args.max_indel))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    check_prefixes(args.prefix)
    if args.restore:
        for p in args.prefix:
            restore(args.indir, p)
        return 0
    check_indir(args.indir)        # before --posthoc moves a single original
    s = settings_from(args)
    if not args.genome:
        raise SystemExit("reboundary: --genome is required")
    if len(args.genome) != len(args.prefix):
        raise SystemExit("reboundary: --prefix and --genome must pair up one to one")
    if args.posthoc:
        for p in args.prefix:
            prepare_posthoc(args.indir, p)
    run(args.indir, args.prefix, args.genome, s, args.tools_dir, args.dump_proposals,
        args.cache, args.verbose)
    if args.posthoc:
        for p, gpath in zip(args.prefix, args.genome):
            regenerate(args.indir, p, gpath, plots=not args.no_plots)
    return 0


if __name__ == "__main__":
    sys.exit(main())

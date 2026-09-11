#!/usr/bin/env python3
"""Recover strand for the LTR-RTs ltrquest.annotate's cascade left unstranded.

An opt-in post-processing stage, in two phases around the annotator:

    recover_strand --phase align   aligns, writes {prefix}_strand_recovery.tsv
    annotate                       reads that file as its fourth cascade tier
    gff3                           pools the tables, unchanged
    recover_strand --phase apply   re-orients the depth FASTAs

Splitting it is what keeps the FASTA and the tables from drifting apart. The
apply phase flips a record only to match the `strand` column ltrquest.annotate
actually wrote, so the two agree by construction rather than because two copies
of the cascade happened to reach the same answer.

Signal 1, homology orientation transfer. The unstranded element, in genome-forward
orientation, is aligned against the elements that already carry a strand, held in
coding-sense orientation. The alignment orientation IS the query's genomic strand.
Four independent views are run:
    dc-megablast on the full-length element
    dc-megablast on the internal region between the LTRs
    minimap2     on the full-length element
    minimap2     on the internal region
A view calls only when every one of its own alignments agrees on orientation. Any
disagreement between views vetoes the locus outright. The preset then sets how much
agreement a surviving locus needs.

Signal 2, polypurine tract margin, for what signal 1 cannot reach. Both strand
hypotheses are scored on the same two windows flanking the internal region, and a call
is made only when one side wins by a margin and clears an absolute floor.

Deliberately NOT used: the consensus LTR and the 5' LTR. Both produce confident,
unanimous, kilobase-scale wrong calls on exactly the elements the full-length view
declines. That holds for every aligner and every encoding of the IUPAC codes.

Leave-one-out against LTRquest's own tesorter strand, pooled over five Brassicaceae
genomes, 3616 labelled elements:
  preset         calls  wrong   error   recall
  conservative    2088      4  0.192%    0.577
  balanced        2845     20  0.703%    0.787
  sensitive       3025     32  1.058%    0.837
Every disagreement outside Chis and Ahal is zero. Disagreement also tracks how well
LTRquest could label the element: 0.37% where 5+ protein domains support the label,
2.54% where only one does.

Usage:
  python -m ltrquest.recover_strand --phase align --prefix Athal_chr2_LTRs \
      --genome Athal_chr2.fa [--preset balanced] [--ppt] [-v]
  python -m ltrquest.recover_strand --phase apply --prefix Athal_chr2_LTRs [-v]

Missing inputs degrade rather than abort: nothing left unstranded, or nothing stranded
to transfer from, logs and returns 0, because the pipeline continues past this stage.
A missing aligner is a real fault and exits non-zero, but the wrapper checks for both
aligners while parsing its own arguments, so an opted-in run fails in its first second
rather than after the last expensive stage.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

import pyfaidx

from .annotate import (
    RECOVERY_HEADER,
    UNKNOWN,
    collect_elements,
    discover_depth_tables,
    element_key,
    header_names,
    load_pass2_links,
    load_recovered_strands,
    load_tesorter_strands,
    load_unannotated,
    read_table,
    recovery_sidecar_path,
    resolve_strands,
    target_mode,
    write_table,
)
# The one true record-level flip. detect.revcomp reverses the IUPAC depth codes
# without complementing them, which is exactly what a depth FASTA needs; see
# `reorient_fasta` for why complementing them would be a silent corruption.
from .detect import revcomp as revcomp_record
from .gff3 import is_gzip
from .table import Columns, as_int

PRESETS = ("conservative", "balanced", "sensitive")
# How many views each preset is supposed to run. The published recall and
# disagreement figures assume all of them, so running short of this is worth
# saying out loud rather than leaving in a -v-only line.
NOMINAL_VIEWS = {"conservative": 2, "balanced": 4, "sensitive": 4}
PHASES = ("align", "apply")

# Tier names this stage writes into the sidecar's `source` column. They join
# 'tesorter', 'domain_order' and 'pass2' in ltrquest.annotate's one vocabulary,
# and reach the user as the GFF3's strand_source attribute.
HOMOLOGY = "homology"
PPT = "ppt"

# An internal region shorter than this is not worth aligning: below roughly a
# hundred bases dc-megablast's seeding is unreliable and the view contributes
# noise rather than evidence.
MIN_INTERNAL_BP = 100

# Full IUPAC complement, for sequence drawn from the genome. This is what the
# leave-one-out numbers above were measured with, and it must never be used on a
# depth FASTA record -- see `reorient_fasta`.
COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv",
                     "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def rc(seq: str) -> str:
    """Reverse complement of genomic sequence, IUPAC codes included."""
    return seq.translate(COMP)[::-1]


def warn(msg: str) -> None:
    print(f"[ltr_recover_strand] WARNING: {msg}", file=sys.stderr)


def log(msg: str) -> None:
    print(f"[ltr_recover_strand] {msg}")


# -----------------------------
# Genome access
# -----------------------------
def plain_genome(genome: str, workdir: str) -> str:
    """A seekable copy of `genome`, decompressing it only if it really is gzipped.

    The driver names the worker's genome '{prefix}.input_genome.fa' whatever its
    compression, which is why gff3.is_gzip sniffs magic bytes rather than
    trusting the extension. faidx cannot index a plain gzip stream, so a
    compressed genome is expanded once into this stage's own work directory and
    removed with it.
    """
    if not is_gzip(genome):
        return genome
    out = os.path.join(workdir, "genome.fa")
    log(f"genome is gzipped; expanding once into {os.path.basename(out)}")
    with gzip.open(genome, "rb") as fin, open(out, "wb") as fout:
        shutil.copyfileobj(fin, fout, 1 << 20)
    return out


class Genome:
    """Indexed random access to the genome. Never resident.

    The index goes in the work directory rather than beside the genome: the
    genome here is the user's own input, routinely on a read-only reference
    mount, and this stage has no business leaving a '.fai' next to it.
    """

    def __init__(self, path: str, workdir: str):
        plain = plain_genome(path, workdir)
        self.fa = pyfaidx.Fasta(
            plain,
            indexname=os.path.join(workdir, os.path.basename(plain) + ".fai"),
            sequence_always_upper=True, as_raw=True)

    def get(self, chrom: str, start: int, end: int) -> Optional[str]:
        """1-based inclusive slice, uppercased. None if the interval is unusable."""
        try:
            seq = self.fa[chrom][max(0, start - 1):end]
        except (KeyError, ValueError):
            return None
        return str(seq) or None

    def close(self) -> None:
        self.fa.close()


# -----------------------------
# Elements
# -----------------------------
class Element(NamedTuple):
    key: str            # 'chrom:start-end', the cascade's own element key
    name: str           # column 1 as written, classification suffix included
    chrom: str
    start: int          # 1-based inclusive, genomic
    end: int
    ltr5_end: int       # 1-based, relative to `start`; 0 when unknown
    ltr3_start: int     # 1-based, relative to `start`; 0 when unknown

    def internal_bounds(self) -> Optional[Tuple[int, int]]:
        """0-based half-open cut of the region between the LTRs, or None.

        Derived from ltr5_end/ltr3_start, the two fields ltrquest.gff3 uses to
        place the long_terminal_repeat features and explicitly not from
        ltr5_len/ltr3_len, which do not stay mutually consistent. The result is
        exactly the span between those two features, and matches the
        `ltr3_start - ltr5_end - 1` internal length plot_summary reports.

        Both fields describe the GENOME-FORWARD frame even for elements whose
        library record is stored reverse-complemented: detect.rebase_to_trimmed
        shifts them by the 5' flank and nothing else re-frames them. That is why
        this stage draws sequence from the genome and never from the depth
        FASTA, where the same offsets would be wrong for every '-' record.
        """
        elem_len = self.end - self.start + 1
        if self.ltr5_end <= 0 or self.ltr3_start <= 0:
            return None
        if self.ltr3_start > elem_len:
            return None
        lo, hi = min(self.ltr5_end, elem_len), self.ltr3_start - 1
        if hi - lo < MIN_INTERNAL_BP:
            return None
        return lo, hi


def elements_from(loaded: Sequence[Tuple[object, Optional[List[str]], List[List[str]]]],
                  verbose: bool = False) -> Dict[str, Element]:
    """Every element across the depth tables, keyed the way the cascade keys them.

    `loaded` carries rows already put through annotate.load_unannotated, so the
    field offsets match what ltrquest.annotate itself sees whether or not a
    previous run inserted the strand and family columns.
    """
    elements: Dict[str, Element] = {}
    for _table, header, rows in loaded:
        cols = Columns.of(header_names(header))
        for row in rows:
            name = row[0] if row else ""
            key = element_key(name)
            if key is None or key in elements:
                continue
            chrom, _, span = key.rpartition(":")
            lo, _, hi = span.partition("-")
            elements[key] = Element(
                key=key, name=name, chrom=chrom, start=int(lo), end=int(hi),
                ltr5_end=as_int(cols.get(row, "ltr5_end"), 0) or 0,
                ltr3_start=as_int(cols.get(row, "ltr3_start"), 0) or 0,
            )
    if verbose:
        log(f"elements: {len(elements)} across {len(loaded)} depth table(s)")
    return elements


# -----------------------------
# Signal 2: polypurine tract
# -----------------------------
def ppt_score(region: str, lmin: int, lmax: int, pen: float, slack: int) -> float:
    """Best purine-run score in `region`, whose LTR boundary is at its RIGHT end."""
    n = len(region)
    cs = [0] * (n + 1)
    for i, ch in enumerate(region):
        cs[i + 1] = cs[i] + (1 if ch in "AG" else 0)
    best = -1e9
    for end in range(max(0, n - slack), n + 1):
        for length in range(lmin, min(lmax, end) + 1):
            purines = cs[end] - cs[end - length]
            score = purines - pen * (length - purines)
            if score > best:
                best = score
    return best


def ppt_call(seq: str, l5: int, l3: int, window: int, lmin: int, lmax: int,
             pen: float, slack: int, margin: float, floor: float
             ) -> Tuple[Optional[str], float, float]:
    """Score both strand hypotheses on the two windows flanking the internal region."""
    internal = seq[l5:len(seq) - l3]
    if len(internal) < 2 * window + 20:
        return None, 0.0, 0.0
    plus = ppt_score(internal[-window:], lmin, lmax, pen, slack)
    minus = ppt_score(rc(internal[:window]), lmin, lmax, pen, slack)
    if max(plus, minus) < floor:
        return None, plus, minus
    delta = plus - minus
    if delta >= margin:
        return "+", plus, minus
    if delta <= -margin:
        return "-", plus, minus
    return None, plus, minus


# -----------------------------
# Signal 1: homology orientation transfer
# -----------------------------
def votes_to_calls(vote: Dict[str, collections.Counter], min_bp: int,
                   purity: float) -> Dict[str, str]:
    """One view's per-query call: the winning orientation, if it wins cleanly."""
    out: Dict[str, str] = {}
    for key, counter in vote.items():
        top, top_bp = counter.most_common(1)[0]
        total = sum(counter.values())
        if top_bp >= min_bp and top_bp / total >= purity:
            out[key] = top
    return out


def parse_paf(text: str) -> Dict[str, collections.Counter]:
    """Query -> matching bases per orientation, from minimap2 PAF."""
    vote: Dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 10 or f[0] == f[5]:
            continue
        vote[f[0]][f[4]] += int(f[9])
    return vote


def parse_blast(text: str) -> Dict[str, collections.Counter]:
    """Query -> matching bases per orientation, from blastn outfmt 6."""
    vote: Dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 4 or f[0] == f[1]:
            continue
        vote[f[0]]["+" if f[2] == "plus" else "-"] += int(f[3])
    return vote


def run_minimap2(binary: str, ref: str, qry: str, threads: int, verbose: bool
                 ) -> Dict[str, collections.Counter]:
    cmd = [binary, "-x", "map-ont", "-k", "11", "-w", "5", "-c", "-N", "50",
           "-p", "0.05", "--secondary=yes", "-m", "20", "-t", str(threads), ref, qry]
    if verbose:
        log("  " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("minimap2 failed: " + proc.stderr.strip().split("\n")[-1])
    return parse_paf(proc.stdout)


def run_blast(ref: str, qry: str, threads: int, evalue: str, workdir: str,
              verbose: bool) -> Dict[str, collections.Counter]:
    db = os.path.join(workdir, "refdb_" + os.path.basename(ref))
    subprocess.run(["makeblastdb", "-in", ref, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    cmd = ["blastn", "-task", "dc-megablast", "-query", qry, "-db", db,
           "-evalue", evalue, "-num_threads", str(threads),
           "-max_target_seqs", "50",
           "-outfmt", "6 qseqid sseqid sstrand nident"]
    if verbose:
        log("  " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("blastn failed: " + proc.stderr.strip().split("\n")[-1])
    return parse_blast(proc.stdout)


def views_needed(preset: str) -> int:
    """How many agreeing views a call needs.

    Only `balanced` asks for two. `conservative` runs just the two dc-megablast
    views and `sensitive` deliberately accepts a single one, since in both cases
    any second view that disagrees has already vetoed the locus.
    """
    return 2 if preset == "balanced" else 1


def combine_views(views: Sequence[Tuple[str, Dict[str, str]]], need: int
                  ) -> Tuple[Dict[str, Tuple[str, str, str]], int, int]:
    """Fold per-view calls into final ones. Returns (calls, vetoed, short).

    Views that contradict each other veto the locus outright: a locus one view
    reads forward and another reads reverse has no defensible call, and dropping
    it is what keeps the error rate where the leave-one-out table says it is.
    """
    tally: Dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for _label, calls in views:
        for key, value in calls.items():
            tally[key][value] += 1

    out: Dict[str, Tuple[str, str, str]] = {}
    vetoed = short = 0
    for key, counter in tally.items():
        if len(counter) > 1:
            vetoed += 1
            continue
        value, agreeing = counter.most_common(1)[0]
        if agreeing < need:
            short += 1
            continue
        out[key] = (value, HOMOLOGY, f"{agreeing} of {len(views)} views")
    return out, vetoed, short


# -----------------------------
# The sidecar ltrquest.annotate reads back
# -----------------------------
def write_sidecar(path: str, unstranded: Sequence[str],
                  calls: Dict[str, Tuple[str, str, str]]) -> None:
    """One row per locus this stage was asked about, called or not.

    Recording the declines too is what makes the file a report as well as an
    input: the difference between 'no call' and 'never looked' is exactly what a
    user comparing presets wants to see. It is also what lets the apply phase
    put a locus back the way it found it when a later run declines it.
    """
    rows = []
    for key in unstranded:
        if key in calls:
            strand, source, evidence = calls[key]
            rows.append([key, strand, source, evidence])
        else:
            rows.append([key, UNKNOWN, "none", ""])
    write_table(path, ["#" + RECOVERY_HEADER[0], *RECOVERY_HEADER[1:]], rows)


# -----------------------------
# Phase 2: re-orienting the depth outputs
# -----------------------------
def depth_fasta_for(tsv_path: str) -> str:
    """The depth FASTA beside a depth TSV: '..._ltr.tsv' -> '..._ltr.fa'."""
    return tsv_path[:-len(".tsv")] + ".fa"


def reorient_fasta(path: str, flip: Sequence[str], wrap: int = 60) -> List[str]:
    """Reverse-complement the named records in place. Returns the names flipped.

    Returning the names rather than a count is what lets the caller restate the
    `orientation` column for exactly the records that moved, so a name the FASTA
    turns out not to hold cannot leave the column claiming a flip that never
    happened.

    Uses detect.revcomp, the same function that oriented these records in the
    first place, and NOT this module's `rc`. A depth FASTA carries nested
    children painted as IUPAC codes whose identity IS the child's depth
    (reconcile.IUPAC_DEPTH_SEQ), so complementing them would relabel depth 1 as
    depth 3, depth 2 as depth 9, and depth 8 as 'V' -- a character reserved for
    the wrapper's far-mask and absent from that table entirely. detect.revcomp
    complements only ACGTN, which reverses the codes without renaming them.

    Rewritten through a temporary file in the same directory and renamed, so an
    interrupted run cannot leave a half-flipped FASTA behind. A file holding
    none of the named records is left exactly as it was, rather than rewritten
    at this function's own line width: the same flip list is offered to every
    depth, and only one of them holds any given element.
    """
    targets = set(flip)
    if not targets or not os.path.isfile(path):
        return []

    directory = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + ".", suffix=".tmp",
                               dir=directory)
    flipped: List[str] = []
    try:
        with os.fdopen(fd, "w") as out, open(path) as fh:
            name: Optional[str] = None
            chunks: List[str] = []

            def emit() -> None:
                if name is None:
                    return
                seq = "".join(chunks)
                if name in targets:
                    seq = revcomp_record(seq)
                    flipped.append(name)
                out.write(f">{name}\n")
                for i in range(0, len(seq), wrap):
                    out.write(seq[i:i + wrap] + "\n")

            for line in fh:
                if line.startswith(">"):
                    emit()
                    # A CRLF file would otherwise keep the '\r' here and match
                    # nothing, silently leaving the record unflipped.
                    name = line[1:].rstrip("\r\n")
                    chunks = []
                else:
                    chunks.append(line.strip())
            emit()
        if flipped:
            os.chmod(tmp, target_mode(path))
            os.replace(tmp, path)
        else:
            os.unlink(tmp)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    return flipped


def restate_orientation(path: str, header: Optional[List[str]],
                        rows: List[List[str]], desired: Dict[str, str]) -> int:
    """Rewrite the `orientation` column so it still reads what the FASTA stores.

    The twin of detect.restate_orientation_to_match_library, which does the same
    job for TEsorter2's own flips. reconcile._outer_is_revcomped treats this
    column as a direct readout of the library, so the two must never drift.
    """
    cols = Columns.of(header_names(header))
    i_orient = cols.index("orientation")
    if i_orient is None:
        warn(f"no 'orientation' column in {os.path.basename(path)}; it cannot be "
             f"kept in step with the FASTA")
        return 0

    changed = 0
    for row in rows:
        key = element_key(row[0]) if row else None
        want = desired.get(key) if key else None
        if want is None or i_orient >= len(row) or row[i_orient] == want:
            continue
        row[i_orient] = want
        changed += 1
    if changed:
        write_table(path, header, rows)
    return changed


def apply_recovery(prefix: str, indir: str = ".", verbose: bool = False) -> int:
    """Re-orient the depth FASTAs so each record matches its final `strand`.

    Driven by the strand column ltrquest.annotate wrote, not by this stage's own
    calls, so the stored orientation and the table cannot disagree. Loci the
    sidecar never mentions are left entirely alone: their orientation belongs to
    TEsorter2's own flip in detect.bounded_fasta_oriented.

    Each depth table is judged against its OWN `orientation` column and paired
    with its own FASTA. Pooling the decision across tables would mean that an
    interrupted run -- one pair rewritten, the next not yet -- read the stale
    pair's column on the retry and flipped the finished one a second time,
    leaving a record storing the opposite of what its column claims.

    Idempotent: the `orientation` column is the gate, so a record already stored
    the way its table asks for is skipped and a second run is a no-op.
    """
    tables = discover_depth_tables(prefix, indir)
    if not tables:
        print(f"[ltr_recover_strand] ERROR: no {prefix}_depth<N>[_clean]_ltr.tsv "
              f"found in {indir}", file=sys.stderr)
        return 1

    owned = load_recovered_strands(prefix, indir, verbose=verbose,
                                   warn_missing=False, called_only=False)
    if not owned:
        if os.path.isfile(recovery_sidecar_path(prefix, indir)):
            log("the strand-recovery sidecar names no loci; nothing to re-orient")
        else:
            log("no strand-recovery sidecar to apply; leaving the depth FASTAs alone")
        return 0

    n_flipped = n_columns = n_skipped = 0
    for table in tables:
        header, rows = read_table(table.path)
        cols = Columns.of(header_names(header))
        fasta = depth_fasta_for(table.path)
        name = os.path.basename(table.path)

        if "strand" not in cols:
            warn(f"{name} carries no 'strand' column; run ltrquest.annotate "
                 f"before --phase apply. Not re-oriented.")
            n_skipped += 1
            continue
        if "orientation" not in cols:
            # Flipping without being able to record it would make the next run
            # flip again, so decline the table rather than corrupt it.
            warn(f"{name} carries no 'orientation' column, so a flip could not "
                 f"be recorded. Not re-oriented.")
            n_skipped += 1
            continue
        if not os.path.isfile(fasta):
            warn(f"no {os.path.basename(fasta)} beside {name}; not re-oriented")
            n_skipped += 1
            continue

        desired: Dict[str, str] = {}
        flip_names: List[str] = []
        for row in rows:
            key = element_key(row[0]) if row else None
            if key is None or key not in owned:
                continue
            strand = cols.get(row, "strand", UNKNOWN)
            want = strand if strand in ("+", "-") else "+"
            desired[key] = want
            if cols.get(row, "orientation", "+") != want:
                flip_names.append(row[0])
        if not flip_names:
            continue

        # FASTA first, then the column that describes it: detect.py orients the
        # library and restates the column in that order for the same reason.
        flipped = reorient_fasta(fasta, flip_names)
        if not flipped:
            continue
        n_flipped += len(flipped)
        moved = {element_key(n): desired[element_key(n)] for n in flipped}
        n_columns += restate_orientation(table.path, header, rows, moved)

    if n_flipped:
        log(f"re-oriented {n_flipped} record(s) across {len(tables)} depth "
            f"table(s) to coding sense; restated {n_columns} orientation value(s)")
    elif not n_skipped:
        log(f"depth FASTAs already match their tables for all {len(owned)} "
            f"recovery-owned locus/loci")
    return 0


# -----------------------------
# Phase 1: alignment
# -----------------------------
def build_work_fastas(elements: Dict[str, Element], strand: Dict[str, str],
                      genome: Genome, workdir: str, verbose: bool = False
                      ) -> Tuple[Dict[str, str], Dict[str, str], int, int]:
    """Write the four alignment FASTAs.

    References are the already-stranded elements held in coding sense; queries
    are the unstranded ones in genome-forward orientation. That asymmetry is the
    whole method: an alignment's orientation against a coding-sense reference is
    the query's own genomic strand.

    Returns (paths, forward sequence of each query, internal reference count,
    internal query count).
    """
    paths = {name: os.path.join(workdir, name + ".fa")
             for name in ("full_ref", "full_qry", "int_ref", "int_qry")}
    handles = {name: open(path, "w") for name, path in paths.items()}
    seqs: Dict[str, str] = {}
    n_int_ref = n_int_qry = n_missing = 0

    try:
        for key, element in elements.items():
            seq = genome.get(element.chrom, element.start, element.end)
            if not seq:
                n_missing += 1
                continue
            bounds = element.internal_bounds()
            internal = None
            if bounds is not None and bounds[1] <= len(seq):
                internal = seq[bounds[0]:bounds[1]]

            if strand.get(key, UNKNOWN) in ("+", "-"):
                sense = seq if strand[key] == "+" else rc(seq)
                handles["full_ref"].write(f">{key}\n{sense}\n")
                if internal is not None:
                    oriented = internal if strand[key] == "+" else rc(internal)
                    handles["int_ref"].write(f">{key}\n{oriented}\n")
                    n_int_ref += 1
            else:
                seqs[key] = seq
                handles["full_qry"].write(f">{key}\n{seq}\n")
                if internal is not None:
                    handles["int_qry"].write(f">{key}\n{internal}\n")
                    n_int_qry += 1
    finally:
        for handle in handles.values():
            handle.close()

    if n_missing:
        warn(f"{n_missing} element(s) had no sequence at their coordinates in the "
             f"genome and were skipped")
    if verbose:
        log(f"  internal-region view: {n_int_ref} references, {n_int_qry} queries")
    return paths, seqs, n_int_ref, n_int_qry


def run_views(paths: Dict[str, str], preset: str, minimap2: Optional[str],
              n_int_qry: int, no_internal: bool, threads: int, evalue: str,
              workdir: str, verbose: bool = False
              ) -> List[Tuple[str, Dict[str, collections.Counter]]]:
    """Run each view the preset asks for, in the order the presets were tuned in.

    A view whose aligner fails is dropped with a warning rather than taking the
    run down: the remaining views still produce a defensible, if smaller, call
    set, and this stage sits after every expensive one.
    """
    want_mm = preset in ("balanced", "sensitive") and minimap2 is not None
    want_internal = not no_internal and n_int_qry > 0

    plan = []
    if want_mm:
        plan.append(("minimap2/full", "minimap2", "full_ref", "full_qry"))
    plan.append(("dc-megablast/full", "blast", "full_ref", "full_qry"))
    if want_internal:
        plan.append(("dc-megablast/internal", "blast", "int_ref", "int_qry"))
    if want_mm and want_internal:
        plan.append(("minimap2/internal", "minimap2", "int_ref", "int_qry"))

    views: List[Tuple[str, Dict[str, collections.Counter]]] = []
    for label, tool, ref, qry in plan:
        try:
            if tool == "minimap2":
                vote = run_minimap2(minimap2, paths[ref], paths[qry], threads,
                                    verbose)
            else:
                vote = run_blast(paths[ref], paths[qry], threads, evalue, workdir,
                                 verbose)
        except (RuntimeError, OSError, subprocess.CalledProcessError) as exc:
            warn(f"the {label} view failed and was dropped: {exc}")
            continue
        views.append((label, vote))

    nominal = NOMINAL_VIEWS[preset]
    if len(views) < nominal:
        ran = ", ".join(label for label, _ in views) or "none"
        warn(f"the {preset} preset runs {nominal} views, but only {len(views)} "
             f"could run here ({ran}). Its published recall and disagreement "
             f"figures assume all {nominal}; with fewer there is less to veto a "
             f"wrong call.")
    return views


def resolve_aligners(preset: str, minimap2: str) -> Optional[str]:
    """Fail fast on a missing aligner: the stage was asked for explicitly.

    Both tools are declared in environment.yml, so a miss here means a broken
    environment rather than an unusual one, and silently degrading the preset
    would hand back a different recall than the one that was requested. The
    wrapper runs the same check while parsing its arguments, so a pipeline run
    fails immediately rather than after every expensive stage has succeeded.
    """
    for tool in ("blastn", "makeblastdb"):
        if shutil.which(tool) is None:
            sys.exit(f"[ltr_recover_strand] ERROR: {tool} is not on PATH. Every "
                     f"strand-recovery preset uses dc-megablast; install blast+ "
                     f"(it is in environment.yml).")
    if preset not in ("balanced", "sensitive"):
        return None
    found = shutil.which(minimap2) or (minimap2 if os.path.exists(minimap2) else None)
    if found is None:
        sys.exit(f"[ltr_recover_strand] ERROR: {minimap2} is not on PATH. The "
                 f"'{preset}' preset needs it for two of its four views; install "
                 f"minimap2 (it is in environment.yml) or use the conservative "
                 f"preset.")
    return found


def align(prefix: str, genome: str, indir: str = ".",
          preset: str = "conservative", ppt: bool = False, threads: int = 8,
          minimap2: str = "minimap2", min_ident_bp: int = 100,
          purity: float = 1.0, evalue: str = "1e-5", no_internal: bool = False,
          ppt_margin: float = 6.0, ppt_floor: float = 12.0, ppt_window: int = 60,
          keep_tmp: bool = False, verbose: bool = False) -> int:
    """Call strand for the unstranded elements and publish the sidecar."""
    tables = discover_depth_tables(prefix, indir)
    if not tables:
        print(f"[ltr_recover_strand] ERROR: no {prefix}_depth<N>[_clean]_ltr.tsv "
              f"found in {indir}", file=sys.stderr)
        return 1
    if not os.path.isfile(genome):
        print(f"[ltr_recover_strand] ERROR: genome not found: {genome}",
              file=sys.stderr)
        return 1

    # One read of the depth tables serves both the geometry this stage needs and
    # the cascade it reproduces to know what already counts as stranded. The
    # element set is built exactly as ltrquest.annotate builds it -- every
    # variant of every depth -- so the two cascades cannot disagree.
    loaded = [(table,) + load_unannotated(table.path) for table in tables]
    elements = elements_from(loaded, verbose)
    info = collect_elements((header, rows) for _t, header, rows in loaded)

    tesorter = load_tesorter_strands(prefix, indir, verbose)
    target_of, orientation = load_pass2_links(prefix, indir, verbose)
    strand, _source = resolve_strands(info, tesorter, target_of, orientation,
                                      verbose)

    unstranded = [k for k in elements if strand.get(k, UNKNOWN) not in ("+", "-")]
    n_known = len(elements) - len(unstranded)
    log(f"elements: {len(elements)}   stranded: {n_known}   "
        f"unstranded: {len(unstranded)}")

    out = recovery_sidecar_path(prefix, indir)
    calls: Dict[str, Tuple[str, str, str]] = {}
    n_homology = vetoed = 0

    short = 0
    if not unstranded:
        log("every element already has a strand; nothing to recover")
    elif n_known == 0 and not ppt:
        warn("no stranded elements to transfer orientation from; "
             "leaving every element unstranded")
    else:
        # The polypurine tract is a property of the element's own sequence, so
        # it still has something to say when there are no donors to align
        # against; homology transfer does not.
        found_mm = resolve_aligners(preset, minimap2) if n_known else None
        workdir = os.path.join(indir, prefix + "_strand_recovery.work")
        os.makedirs(workdir, exist_ok=True)
        fa = None
        try:
            fa = Genome(genome, workdir)
            paths, seqs, _n_int_ref, n_int_qry = build_work_fastas(
                elements, strand, fa, workdir, verbose)

            if not n_known:
                warn("no stranded elements to transfer orientation from; "
                     "running the PPT fallback alone")
            else:
                views = run_views(paths, preset, found_mm, n_int_qry, no_internal,
                                  threads, evalue, workdir, verbose)
                if not views:
                    warn("no alignment view could be run; homology transfer "
                         "contributed nothing")

                per_view: List[Tuple[str, Dict[str, str]]] = []
                for label, vote in views:
                    called = votes_to_calls(vote, min_ident_bp, purity)
                    if verbose:
                        log(f"  {label:<22} aligned {len(vote):4d}, "
                            f"called {len(called):4d}")
                    per_view.append((label, called))

                calls, vetoed, short = combine_views(per_view, views_needed(preset))
                n_homology = len(calls)

            if ppt:
                for key in unstranded:
                    bounds = elements[key].internal_bounds()
                    seq = seqs.get(key)
                    # The same guard build_work_fastas applies: an element
                    # running past the contig end comes back truncated, and
                    # scoring the 3' window of a truncated record would put the
                    # tract search inside the LTR.
                    if (key in calls or bounds is None or seq is None
                            or bounds[1] > len(seq)):
                        continue
                    call, plus, minus = ppt_call(
                        seq, bounds[0], len(seq) - bounds[1], ppt_window,
                        12, 30, 3.0, 2, ppt_margin, ppt_floor)
                    if call:
                        calls[key] = (call, PPT,
                                      f"score + {plus:.0f} vs - {minus:.0f}")
        finally:
            if fa is not None:
                fa.close()
            if not keep_tmp:
                shutil.rmtree(workdir, ignore_errors=True)
            elif verbose:
                log(f"  kept the working FASTAs in {os.path.basename(workdir)}")

    write_sidecar(out, unstranded, calls)
    dropped = f"   single-view {short}" if short else ""
    log(f"recovered {len(calls)} of {len(unstranded)}   "
        f"[homology {n_homology}, ppt {len(calls) - n_homology}]   "
        f"vetoed {vetoed}{dropped}   "
        f"still unstranded {len(unstranded) - len(calls)}")
    log(f"wrote {os.path.basename(out)}")
    return 0


# -----------------------------
# Driver
# -----------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Recover strand for the LTR-RTs ltrquest.annotate's cascade "
                    "left unstranded, and re-orient the depth FASTAs to match.",
        # argparse %-formats help= strings but not the epilog, so this one
        # carries single percent signs while --ppt's help below doubles them.
        epilog="Presets, measured by leave-one-out against LTRquest's own "
               "tesorter strand over 3616 elements: conservative, the two "
               "dc-megablast views, recall 0.58 at 0.19% disagreement; "
               "balanced, all four views with two required to agree, recall "
               "0.79 at 0.70%; sensitive, all four views with one enough "
               "unless another contradicts it, recall 0.84 at 1.06%.")
    parser.add_argument("--prefix", required=True,
                        help="Wrapper output prefix, e.g. Athal_tair10_chr2_LTRs")
    parser.add_argument("--phase", choices=PHASES, default="align",
                        help="align: call strand and write the sidecar, before "
                             "ltrquest.annotate. apply: re-orient the depth "
                             "FASTAs to the strand column, after it. "
                             "(default: align)")
    parser.add_argument("--genome",
                        help="Genome FASTA the elements were called on (plain or "
                             "gzipped). Required by --phase align; indexed into "
                             "the work directory, never held in memory.")
    parser.add_argument("--indir", default=".",
                        help="Directory holding the wrapper outputs (default: .)")
    parser.add_argument("--preset", choices=PRESETS, default="conservative",
                        help="How much agreement a call needs (default: "
                             "conservative)")
    parser.add_argument("--ppt", action="store_true",
                        help="Add the polypurine-tract fallback for loci homology "
                             "cannot reach. Measured at 93%% on real genomes, so "
                             "off by default.")
    parser.add_argument("--no-internal", action="store_true",
                        help="Skip the internal-region views, leaving only the "
                             "full-length ones")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Threads for the aligners (default: 8)")
    parser.add_argument("--minimap2", default="minimap2",
                        help="minimap2 binary (default: minimap2). Any build "
                             "works; the input is pure ACGT.")
    parser.add_argument("--min-ident-bp", type=int, default=100,
                        help="Matching bases a view needs before it calls "
                             "(default: 100)")
    parser.add_argument("--purity", type=float, default=1.0,
                        help="Fraction of a view's matching bases that must agree "
                             "on orientation (default: 1.0)")
    parser.add_argument("--evalue", default="1e-5",
                        help="blastn evalue (default: 1e-5)")
    parser.add_argument("--ppt-margin", type=float, default=6.0,
                        help="PPT score margin (default: 6)")
    parser.add_argument("--ppt-floor", type=float, default=12.0,
                        help="PPT winning-score floor (default: 12)")
    parser.add_argument("--ppt-window", type=int, default=60,
                        help="bp searched at each end for the PPT (default: 60)")
    parser.add_argument("--keep-tmp", action="store_true",
                        help="Keep the working FASTAs and the BLAST database")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Per-view call counts and the commands run")
    args = parser.parse_args(argv)

    if args.phase == "apply":
        return apply_recovery(args.prefix, args.indir, args.verbose)
    if not args.genome:
        parser.error("--phase align needs --genome")
    return align(
        prefix=args.prefix, genome=args.genome, indir=args.indir,
        preset=args.preset, ppt=args.ppt, threads=args.threads,
        minimap2=args.minimap2, min_ident_bp=args.min_ident_bp,
        purity=args.purity, evalue=args.evalue, no_internal=args.no_internal,
        ppt_margin=args.ppt_margin, ppt_floor=args.ppt_floor,
        ppt_window=args.ppt_window, keep_tmp=args.keep_tmp,
        verbose=args.verbose)


if __name__ == "__main__":  # pragma: no cover - exercised via the console script
    sys.exit(main())

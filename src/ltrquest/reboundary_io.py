"""Reading and rewriting a run's clean side for ltrquest.reboundary.

Only the `_clean_` depth tables and FASTAs are ever rewritten. The raw
`_depth<N>_ltr.*` files are detection outputs that ltrquest.record hands to a
later run, so they keep the calls as detected.

Rewriting is a two-phase commit (`Commit`). Every new file is written to
`<path>.new` and fsynced before any file is swapped in, so a failure while the
files are being written -- the long phase -- leaves the run exactly as it was.
The swap phase then moves each target aside to `<path>.old` and its replacement
into place, one path at a time, and deletes the `.old` files only once every
swap has succeeded; a failure part-way through puts every original back. What a
crash during the swap phase itself can leave behind is `.new`/`.old` files
beside their targets, and the originals are always recoverable from them by
hand. A kill between the two renames of one swap takes the target away
altogether, so nothing a later run loads can reveal it: `leftover_staging` looks
at the directory itself, and the driver calls it before any work begins.
"""

from __future__ import annotations

import bisect
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import (Callable, Dict, FrozenSet, Iterable, Iterator, List, NamedTuple, Optional,
                    Sequence, Set, Tuple)

from .annotate import discover_depth_tables, element_key, header_names, read_table, target_mode
from .detect import revcomp as revcomp_record
from .kmer2ltr import COLUMNS
from .ltr_model import Member
from .table import Columns, as_float, as_int

SIDECAR_SUFFIX = "_reboundary.tsv"
SIDECAR_COLUMNS = [
    "old_seq_id", "new_seq_id", "family", "method", "templates", "decision", "reason",
    "ext5", "ext3", "end_source5", "end_source3", "id_outer5", "id_outer3", "support5",
    "support3", "k2l_status", "tsd_called", "tsd_new", "tsd_null", "k2p_called", "k2p_new",
    "unpaired5", "unpaired3", "obstacle5", "obstacle3", "merged_into",
]
# Decisions whose row names a re-bounded element ("extended": sidecars written before
# trims and merges existed).
MOVED = ("moved", "extended")
NEST_RE = re.compile(r"^(nest-outer|nest-inner):(.+)$")


@dataclass
class CleanTable:
    path: str
    fasta: str
    depth: int
    header: List[str]
    rows: List[List[str]]
    cols: Columns


def load_clean_tables(indir: str, prefix: str) -> List[CleanTable]:
    out: List[CleanTable] = []
    for t in discover_depth_tables(prefix, indir, variants=("clean",)):
        header, rows = read_table(t.path)
        if header is None:
            raise ValueError(f"{t.path}: no header line")
        cols = Columns.of(header_names(header))
        if cols.names[:len(COLUMNS)] != COLUMNS:
            raise ValueError(f"{t.path}: its first {len(COLUMNS)} columns are not Kmer2LTR's, "
                             f"and re-boundarying replaces them by position")
        missing = [c for c in ("strand", "family", "nest_status", "orientation") if c not in cols]
        if missing:
            raise ValueError(f"{t.path}: no {', '.join(missing)} column; run ltrquest.annotate "
                             f"first (re-boundarying needs the strand, family, nesting and "
                             f"orientation columns)")
        out.append(CleanTable(t.path, t.path[:-len(".tsv")] + ".fa", t.depth, header, rows, cols))
    return out


def members_from(prefix: str, tables: Sequence[CleanTable]) -> List[Member]:
    out: List[Member] = []
    for t in tables:
        for row in t.rows:
            key = element_key(row[0]) if row else None
            if key is None:
                continue
            chrom, span = key.rsplit(":", 1)
            start, end = (int(x) for x in span.split("-"))
            l5e, l3s = as_int(t.cols.get(row, "ltr5_end")), as_int(t.cols.get(row, "ltr3_start"))
            if l5e is None or l3s is None:
                continue
            out.append(Member(
                prefix=prefix, name=row[0], chrom=chrom, start=start, end=end,
                l1=start + l5e - 1, r0=start + l3s - 1,
                strand=t.cols.get(row, "strand"), orientation=t.cols.get(row, "orientation", "+"),
                family=t.cols.get(row, "family"), depth=t.depth,
                k2p=as_float(t.cols.get(row, "k2p")), tsd=t.cols.get(row, "tsd"),
                nest_status=t.cols.get(row, "nest_status"),
                tsd_offset=t.cols.get(row, "tsd_offset", "")))
    return out


def iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].rstrip("\r\n"), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def fetch_records(paths: Iterable[str], names: Set[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for path in paths:
        if names and os.path.isfile(path):
            for name, seq in iter_fasta(path):
                if name in names:
                    out[name] = seq
    return out


def forward(record: str, orientation: str) -> str:
    """A stored record turned to the genome-forward frame (depth letters untouched)."""
    return record if orientation != "-" else revcomp_record(record)


def stored(fwd: str, orientation: str) -> str:
    return fwd if orientation != "-" else revcomp_record(fwd)


def sanitize(seq: str) -> str:
    """Every non-ACGT base as N, as Kmer2LTR's own FASTA reader does."""
    return re.sub(r"[^ACGT]", "N", seq.upper())


def paint(record: str, orientation: str, rec_start: int, spans: Sequence[Tuple[int, int]],
          letter: str) -> str:
    """Overwrite genomic spans (1-based, inclusive) inside a stored record with a depth letter."""
    chars, n = list(record), len(record)
    for a, b in spans:
        lo, hi = max(0, a - rec_start), min(n, b - rec_start + 1)
        if hi <= lo:
            continue
        if orientation == "-":
            lo, hi = n - hi, n - lo
        for i in range(lo, hi):
            chars[i] = letter
    return "".join(chars)


def rekey_nest(value: str, old2new: Dict[str, str]) -> str:
    if value in ("", ".") or not old2new:
        return value
    out = []
    for tok in value.split(";"):
        m = NEST_RE.match(tok)
        if m and m.group(2) in old2new:
            tok = f"{m.group(1)}:{old2new[m.group(2)]}"
        out.append(tok)
    return ";".join(out)


def drop_nest(value: str, keys: FrozenSet[str]) -> str:
    """`value` without the tokens that name any of `keys` (calls that no longer exist).

    Dropped rather than re-keyed onto the call that absorbed them: that call's own
    token is already there, and `rekey_nest` does not deduplicate.
    """
    if value in ("", ".") or not keys:
        return value
    out = [tok for tok in value.split(";")
           if not ((m := NEST_RE.match(tok)) and m.group(2) in keys)]
    return ";".join(out) if out else "."


def hosts_of(nest_status: str) -> List[str]:
    """Keys of the elements this one is nested inside."""
    return [m.group(2) for m in (NEST_RE.match(t) for t in nest_status.split(";"))
            if m and m.group(1) == "nest-inner"]


class SpanIndex:
    """Every element's span per (prefix, chrom), for overlap queries."""

    def __init__(self, members: Iterable[Member]):
        by: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = defaultdict(list)
        for m in members:
            by[(m.prefix, m.chrom)].append((m.start, m.end, m.key))
        self._by = {k: sorted(v) for k, v in by.items()}
        self._starts = {k: [s for s, _, _ in v] for k, v in self._by.items()}
        self._longest = {k: max(e - s + 1 for s, e, _ in v) for k, v in self._by.items()}

    def overlapping(self, prefix: str, chrom: str, lo: int, hi: int) -> List[Tuple[int, int, str]]:
        k = (prefix, chrom)
        spans = self._by.get(k, [])
        out = []
        i = bisect.bisect_left(self._starts.get(k, []), lo - self._longest.get(k, 0))
        for s, e, key in spans[i:]:
            if s > hi:
                break
            if e >= lo:
                out.append((s, e, key))
        return out


def _added(m: Member, left: int, right: int) -> List[Tuple[int, int]]:
    return (([(left, m.start - 1)] if left < m.start else [])
            + ([(m.end + 1, right)] if right > m.end else []))


def conflict(index: SpanIndex, m: Member, left: int, right: int,
             ignore: FrozenSet[str] = frozenset()) -> Optional[str]:
    """Why moving `m` to left..right would clash with another element, or None.

    `ignore` holds keys of calls that are leaving (a merge partner).
    """
    added = _added(m, left, right)
    for s, e, key in index.overlapping(m.prefix, m.chrom, min(left, m.start),
                                       max(right, m.end)):
        if key == m.key or key in ignore:
            continue
        if s <= m.start and e >= m.end:           # a host of the call as it stands
            if not (s <= left and e >= right):
                return "host_exceeded"
            if s == left and e == right:          # the new key would be the host's own
                return "duplicates_host"
            continue
        if s >= m.start and e <= m.end:           # nested inside the call already
            if s < left or e > right:             # ...and a trim would leave it outside
                return "nest_broken"
            if s == left and e == right:          # ...or give the call its key
                return "duplicates_element"
            continue
        for a, b in added:
            if s >= a and e <= b:
                return "engulfs_element"
            if s <= b and e >= a:
                return "overlaps_element"
    return None


def merge_partner(index: SpanIndex, by_uid: Dict[str, Member], m: Member, left: int,
                  right: int, tol: int = 5) -> Optional[str]:
    """Key of the one call that is `m`'s own element called a second time, or None.

    A detector can call one element twice, staggered: each call has one of the
    element's true outer ends. Moving `m` to its true ends then reaches the other
    call's outer end, and the other call -- not the move -- is what is wrong. It is
    the same element when it has the same strand, depth and nesting, hosts nothing,
    carries no TSD of its own (whose ends a TSD confirms is a different insertion),
    lies inside the new span, and its outer end on the side it overlaps is the new
    outer end within `tol` bp. Family is not compared: a split call's family is
    clustered from its own, partly wrong, LTR pair. Anything short of exactly one
    such call is no merge.
    """
    added = _added(m, left, right)
    found: List[str] = []
    for s, e, key in index.overlapping(m.prefix, m.chrom, left, right):
        if key == m.key or not any(s <= b and e >= a for a, b in added):
            continue
        if s <= m.start and e >= m.end:            # a host, not a second call
            continue
        x = by_uid.get(f"{m.prefix}\t{key}")
        if x is None or not (left <= s and e <= right):
            return None
        same = (x.strand == m.strand and x.depth == m.depth
                and sorted(hosts_of(x.nest_status)) == sorted(hosts_of(m.nest_status))
                and "nest-outer:" not in (x.nest_status or "")
                and x.tsd in ("", ".", "NA"))
        ends = ((s < m.start and abs(s - left) <= tol)
                or (e > m.end and abs(e - right) <= tol))
        if not (same and ends):
            return None
        found.append(key)
    return found[0] if len(found) == 1 else None


def nested_index(members: Iterable[Member]) -> Dict[Tuple[str, str], List[Member]]:
    """(prefix, host key) -> every element nested anywhere inside that host. Built once,
    it spares `trim_restore` a scan of every member per trimmed element and host."""
    out: Dict[Tuple[str, str], List[Member]] = defaultdict(list)
    for x in members:
        for hk in hosts_of(x.nest_status):
            out[(x.prefix, hk)].append(x)
    return out


def trim_restore(m: Member, trimmed: Sequence[Tuple[int, int]], by_uid: Dict[str, Member],
                 letters: Sequence[str], fetch: Callable[[str, str, int, int], str],
                 nested: Optional[Dict[Tuple[str, str], List[Member]]] = None
                 ) -> Dict[str, Tuple[Tuple[int, str], ...]]:
    """host key -> ((genomic start, forward-frame characters), ...) for `m`'s trimmed bases.

    Every ancestor's record masked those bases with `m`'s depth letter. Once they
    are no longer `m`'s, each shows what reconcile paints there without `m`: the
    letter of the innermost other element of that host still covering it (an
    intermediate host, in a grandparent), else the genomic base.
    """
    out: Dict[str, Tuple[Tuple[int, str], ...]] = {}
    if nested is None:
        nested = nested_index(x for x in by_uid.values() if x.prefix == m.prefix)
    for hk in hosts_of(m.nest_status):
        within = [x for x in nested.get((m.prefix, hk), ()) if x.key != m.key]
        segs = []
        for a, b in trimmed:
            chars = list(fetch(m.prefix, m.chrom, a, b))
            for i, pos in enumerate(range(a, b + 1)):
                cover = [x for x in within if x.start <= pos <= x.end]
                if cover:
                    d = min(x.depth for x in cover)
                    chars[i] = letters[d] if d < len(letters) else "X"
            segs.append((a, "".join(chars)))
        out[hk] = tuple(segs)
    return out


def overwrite(record: str, orientation: str, rec_start: int,
              segments: Sequence[Tuple[int, str]]) -> str:
    """Write forward-frame `segments` (genomic start, characters) into a stored record,
    mirrored and reverse-complemented when the record is stored `-` (depth letters are
    written as they are, as `paint` and reconcile do)."""
    chars, n = list(record), len(record)
    for start, text in segments:
        if orientation == "-":
            first, text = n - (start - rec_start) - len(text), revcomp_record(text)
        else:
            first = start - rec_start
        for i, ch in enumerate(text):
            if 0 <= first + i < n:
                chars[first + i] = ch
    return "".join(chars)


def genome_order(m: Member) -> Tuple[str, str, int, int, str]:
    """Genome, chrom, start order, ties broken by end and name: whatever the order in
    which workers return, the first claim is always the same call's."""
    return m.prefix, m.chrom, m.start, m.end, m.name


def mutual_conflicts(items: Sequence[Tuple[Member, int, int]]) -> Set[str]:
    """uids whose added bases overlap an earlier candidate's added bases
    (`genome_order`)."""
    taken: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    lost: Set[str] = set()
    for m, left, right in sorted(items, key=lambda t: genome_order(t[0])):
        added = _added(m, left, right)
        k = (m.prefix, m.chrom)
        if any(a <= d and b >= c for a, b in added for c, d in taken[k]):
            lost.add(m.uid)
            continue
        taken[k].extend(added)
    return lost


def duplicate_moves(items: Sequence[Tuple[Member, int, int]]) -> Set[str]:
    """uids moving to the span an earlier candidate (`genome_order`) moves to.

    Two calls of one element, one inside the other, can reach the same ends from
    opposite sides, one extending and the other trimming. Neither move conflicts
    with the other call as it stands, but both would be written under one key.
    """
    taken: Set[Tuple[str, str, int, int]] = set()
    lost: Set[str] = set()
    for m, left, right in sorted(items, key=lambda t: genome_order(t[0])):
        span = (m.prefix, m.chrom, left, right)
        if span in taken:
            lost.add(m.uid)
        else:
            taken.add(span)
    return lost


@dataclass(frozen=True)
class Accepted:
    member: Member
    new_name: str
    fields: List[str]       # the 29 Kmer2LTR columns, rebased
    record: str             # stored-orientation FASTA record
    left: int
    right: int
    # host key -> ((genomic start, forward characters), ...) over trimmed bases
    restore: Dict[str, Tuple[Tuple[int, str], ...]] = field(default_factory=dict)


NEW_SUFFIX = ".new"
OLD_SUFFIX = ".old"
# What this stage ever stages, so an unrelated `genome.fa.old` is not mistaken for ours.
STAGED_RE = re.compile(r"(_depth\d+_clean_ltr\.(tsv|fa)|" + re.escape(SIDECAR_SUFFIX) + r")$")


class CommitError(RuntimeError):
    """A swap failed and the originals could not all be put back automatically."""


def leftover_staging(directory: str) -> List[str]:
    """`.new`/`.old` files a rewrite killed mid-swap left in `directory`, sorted.

    A kill between `os.replace(path, path.old)` and `os.replace(path.new, path)`
    leaves no `path` at all, so the next run never discovers it, never stages it
    and would annotate that depth short without a word. Only the directory
    listing shows it, which is why this runs before any work.
    """
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return []
    return [os.path.join(directory, n) for n in names
            for suffix in (NEW_SUFFIX, OLD_SUFFIX)
            if n.endswith(suffix) and STAGED_RE.search(n[:-len(suffix)])]


def leftover_message(directory: str, leftovers: Sequence[str]) -> str:
    """What is on disk and what to do with it, for a refusal that names every file."""
    return (f"reboundary: {directory} holds files left by a rewrite that was interrupted "
            f"mid-swap: " + ", ".join(os.path.basename(p) for p in leftovers)
            + f". A '{OLD_SUFFIX}' file is the original of the path it is named after and a "
            f"'{NEW_SUFFIX}' file is a rewrite that never went in. Move each '{OLD_SUFFIX}' "
            f"file back over that path, delete the '{NEW_SUFFIX}' files, then re-run. "
            f"Nothing has been changed.")


def _unlink(path: str) -> bool:
    try:
        os.unlink(path)
    except OSError:
        return False
    return True


def sync_dir(path: str) -> None:
    """Make the renames themselves durable; best effort, some filesystems refuse."""
    try:
        fd = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(fd)
    except OSError:
        pass
    finally:
        os.close(fd)


class _Staged:
    """A handle on a staged file whose bytes reach the disk before it closes."""

    def __init__(self, fh):
        self._fh = fh

    def write(self, text: str) -> int:
        return self._fh.write(text)

    def close(self) -> None:
        if self._fh.closed:
            return
        self._fh.flush()
        os.fsync(self._fh.fileno())
        self._fh.close()

    def __enter__(self) -> "_Staged":
        return self

    def __exit__(self, *exc) -> bool:
        self.close()
        return False


class Commit:
    """Write every rewritten file beside its target, then swap them all in.

    `open` hands back a handle on `<path>.new`, fsynced when it closes, so a
    failure during writing changes nothing: `abort()` just deletes the `.new`
    files. `commit()` then swaps, per path in staging order, `<path>` out to
    `<path>.old` and `<path>.new` in, deleting the `.old` files only once every
    path has been swapped. A failure mid-swap restores every `.old` in reverse
    order and deletes the `.new` files, so the run is left exactly as it was; if
    a restore itself fails, nothing is deleted and `CommitError` names every
    file left behind and what to do with it. A leftover `.old` from a swap that
    a kill cut short stops the next `commit()` rather than being overwritten.

    The staging names are derived from the target, so one path may be staged
    only once: a second swap of the same path would move the first swap's
    rewrite into the `.old` that holds the original, and the original would be
    gone with nothing left to roll back from. `open` refuses it.
    """

    def __init__(self):
        self._staged: List[Tuple[str, str]] = []
        self._paths: Set[str] = set()

    def open(self, path: str):
        if path in self._paths:
            raise ValueError(f"reboundary: {path} was staged twice in one commit; a path may "
                             f"be staged once, and its original would be lost on the second "
                             f"swap. Nothing has been written.")
        tmp = path + NEW_SUFFIX
        self._staged.append((tmp, path))
        self._paths.add(path)
        fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        return _Staged(os.fdopen(fd, "w"))

    def commit(self) -> None:
        stale = [p + OLD_SUFFIX for _, p in self._staged if os.path.lexists(p + OLD_SUFFIX)]
        if stale:                             # an earlier swap was killed and never recovered
            raise CommitError(
                "reboundary: an earlier rewrite was interrupted mid-swap and left "
                + ", ".join(stale) + ". Each one holds the original of the path it is named "
                "after: move it back over that path, delete any '" + NEW_SUFFIX
                + "' file beside it, then re-run.")
        staged, self._staged, self._paths = self._staged, [], set()
        done: List[Tuple[str, str]] = []      # (target, its .old, or "" if it was absent)
        try:
            for tmp, path in staged:
                os.chmod(tmp, target_mode(path))
            for tmp, path in staged:
                old = path + OLD_SUFFIX if os.path.lexists(path) else ""
                if old:
                    os.replace(path, old)
                done.append((path, old))
                os.replace(tmp, path)
        except BaseException as exc:
            _roll_back(done, staged, exc)
            raise
        for directory in sorted({os.path.dirname(os.path.abspath(p)) for _, p in staged}):
            sync_dir(directory)
        for _, old in done:
            if old:
                _unlink(old)

    def abort(self) -> None:
        for tmp, _ in self._staged:
            _unlink(tmp)
        self._staged, self._paths = [], set()


def _roll_back(done: Sequence[Tuple[str, str]], staged: Sequence[Tuple[str, str]],
               exc: BaseException) -> None:
    """Undo the swaps already made. Raise CommitError rather than hide a failed undo.

    The message describes only the files this commit made, by name: every
    `.old` named here was moved aside from its own path moments ago and holds
    that path's original, and every path named under `delete` had no file at all
    before this commit. Nothing is said about `.old` files in general.
    """
    move_back: List[str] = []                 # <path>.old, still holding <path>'s original
    delete: List[str] = []                    # a path this commit created and could not remove
    for path, old in reversed(list(done)):
        try:
            if old:
                os.replace(old, path)
            elif os.path.lexists(path):
                os.unlink(path)               # nothing was there before this commit
        except OSError as e:
            (move_back if old else delete).append(f"{old or path} ({e.strerror})")
    if move_back or delete:
        left = [t for t, _ in staged if os.path.lexists(t)]
        raise CommitError(
            "reboundary: a file swap failed and the originals could not all be put back. "
            "Nothing has been deleted; every file involved is named below, with what to do. "
            + ("Move back over the path it is named after, which now holds a rewrite (the '"
               + OLD_SUFFIX + "' file is that path's original): " + ", ".join(move_back) + ". "
               if move_back else "")
            + ("Delete; this run created it and there was no file there before: "
               + ", ".join(delete) + ". " if delete else "")
            + ("Rewrites still staged, to delete as well: " + ", ".join(left) + ". "
               if left else "")
            + "Then re-run.") from exc
    for tmp, _ in staged:
        _unlink(tmp)


def check_keys_stay_unique(tables: Sequence[CleanTable], old2new: Dict[str, str],
                           retired: FrozenSet[str] = frozenset()) -> None:
    """Refuse a re-key that would give two elements one key, before anything is written.

    A row and a FASTA record are found by their element key, so two rows sharing
    one would make every consumer keyed on it drop an element without a word.
    `conflict` already refuses every proposal that could produce this; the check
    is here so that no other route to it can ever land on disk.
    """
    claimed: Dict[str, str] = {}
    for t in tables:
        for row in t.rows:
            k = element_key(row[0]) if row else None
            if k is None or k in retired:
                continue
            new = old2new.get(k, k)
            if new in claimed:
                raise ValueError(f"reboundary: {row[0]} and {claimed[new]} would both be "
                                 f"written with the same key {new}; nothing was written")
            claimed[new] = row[0]


def rewrite(tables: Sequence[CleanTable], accepted: Dict[str, Accepted], letters: Sequence[str],
            commit: Commit, wrap: int = 60, retired: FrozenSet[str] = frozenset()) -> None:
    """Stage one genome's clean tables and FASTAs with `accepted` (old key -> result)
    applied and the `retired` calls (merged into another) removed."""
    old2new = {k: element_key(a.new_name) for k, a in accepted.items()}
    check_keys_stay_unique(tables, old2new, retired)
    jobs: Dict[str, List[Tuple[List[Tuple[int, int]], str]]] = defaultdict(list)
    restores: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    for a in accepted.values():
        m = a.member
        letter = letters[m.depth] if m.depth < len(letters) else "X"
        for host in hosts_of(m.nest_status):
            jobs[host].append((_added(m, a.left, a.right), letter))
        for host, segs in a.restore.items():
            restores[host].extend(segs)
    where: Dict[str, Tuple[int, str]] = {}
    for t in tables:
        for row in t.rows:
            k = element_key(row[0]) if row else None
            if k is not None:
                if k in accepted:
                    start = accepted[k].left
                else:
                    start = int(k.rsplit(":", 1)[1].split("-")[0])
                where[k] = (start, t.cols.get(row, "orientation", "+"))
    for t in tables:
        i_nest = t.cols.require("nest_status")
        rows = []
        for row in t.rows:
            k = element_key(row[0]) if row else None
            if k in retired:
                continue
            if k in accepted:
                row = list(accepted[k].fields) + row[len(COLUMNS):]
            else:
                row = list(row)
            row[i_nest] = rekey_nest(drop_nest(row[i_nest], retired), old2new)
            rows.append(row)
        if os.path.isfile(t.fasta):
            with commit.open(t.fasta) as out:
                for name, seq in iter_fasta(t.fasta):
                    k = element_key(name)
                    if k in retired:
                        continue
                    if k in accepted:
                        name, seq = accepted[k].new_name, accepted[k].record
                    for spans, letter in jobs.get(k, ()):
                        start, orientation = where[k]
                        seq = paint(seq, orientation, start, spans, letter)
                    if k in restores:
                        start, orientation = where[k]
                        seq = overwrite(seq, orientation, start, restores[k])
                    out.write(f">{name}\n")
                    for i in range(0, len(seq), wrap):
                        out.write(seq[i:i + wrap] + "\n")
        with commit.open(t.path) as out:
            out.write("\t".join(t.header) + "\n")
            for row in rows:
                out.write("\t".join(row) + "\n")


def sidecar_text(rows: Sequence[Dict[str, str]]) -> str:
    lines = ["#" + "\t".join(SIDECAR_COLUMNS)]
    lines += ["\t".join(str(r.get(c, ".")) for c in SIDECAR_COLUMNS) for r in rows]
    return "\n".join(lines) + "\n"


class Rebound(NamedTuple):
    old_name: str
    new_name: str
    ext5: int
    ext3: int


def read_map(path: str) -> Dict[str, Rebound]:
    """New element key -> Rebound, for the extended rows of a sidecar."""
    header, rows = read_table(path)
    cols = Columns.of(header_names(header))
    out: Dict[str, Rebound] = {}
    for row in rows:
        if cols.get(row, "decision") not in MOVED:
            continue
        new = cols.get(row, "new_seq_id")
        key = element_key(new)
        if key is not None:
            out[key] = Rebound(cols.get(row, "old_seq_id"), new,
                               as_int(cols.get(row, "ext5"), 0), as_int(cols.get(row, "ext3"), 0))
    return out

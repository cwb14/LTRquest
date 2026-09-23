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
beside their targets: `commit()` reports them by name if it cannot clean them up
itself, and the originals are always recoverable from them by hand.
"""

from __future__ import annotations

import bisect
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, NamedTuple, Optional, Sequence, Set, Tuple

from .annotate import discover_depth_tables, element_key, header_names, read_table, target_mode
from .detect import revcomp as revcomp_record
from .kmer2ltr import COLUMNS
from .ltr_model import Member
from .table import Columns, as_float, as_int

SIDECAR_SUFFIX = "_reboundary.tsv"
SIDECAR_COLUMNS = [
    "old_seq_id", "new_seq_id", "family", "method", "model_id", "decision", "reason",
    "ext5", "ext3", "end_source5", "end_source3", "id_outer5", "id_outer3", "credit_bits",
    "k2l_status", "tsd_called", "tsd_new", "tsd_null", "k2p_called", "k2p_new",
    "obstacle5", "obstacle3",
]
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
                             f"first (re-boundarying models families)")
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
                nest_status=t.cols.get(row, "nest_status")))
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


def conflict(index: SpanIndex, m: Member, left: int, right: int) -> Optional[str]:
    """Why extending `m` to left..right would clash with another element, or None."""
    added = _added(m, left, right)
    for s, e, key in index.overlapping(m.prefix, m.chrom, left, right):
        if key == m.key:
            continue
        if s <= m.start and e >= m.end:           # a host of the call as it stands
            if not (s <= left and e >= right):
                return "host_exceeded"
            if s == left and e == right:          # the new key would be the host's own
                return "duplicates_host"
            continue
        if s >= m.start and e <= m.end:           # nested inside the call already
            continue
        for a, b in added:
            if s >= a and e <= b:
                return "engulfs_element"
            if s <= b and e >= a:
                return "overlaps_element"
    return None


def mutual_conflicts(items: Sequence[Tuple[Member, int, int]]) -> Set[str]:
    """uids whose added bases overlap an earlier candidate's added bases
    (genome, chrom, start order)."""
    taken: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    lost: Set[str] = set()
    for m, left, right in sorted(items, key=lambda t: (t[0].prefix, t[0].chrom, t[0].start)):
        added = _added(m, left, right)
        k = (m.prefix, m.chrom)
        if any(a <= d and b >= c for a, b in added for c, d in taken[k]):
            lost.add(m.uid)
            continue
        taken[k].extend(added)
    return lost


@dataclass(frozen=True)
class Accepted:
    member: Member
    new_name: str
    fields: List[str]       # the 29 Kmer2LTR columns, rebased
    record: str             # stored-orientation FASTA record
    left: int
    right: int


NEW_SUFFIX = ".new"
OLD_SUFFIX = ".old"


class CommitError(RuntimeError):
    """A swap failed and the originals could not all be put back automatically."""


def _unlink(path: str) -> bool:
    try:
        os.unlink(path)
    except OSError:
        return False
    return True


def _sync_dir(path: str) -> None:
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
    """

    def __init__(self):
        self._staged: List[Tuple[str, str]] = []

    def open(self, path: str):
        tmp = path + NEW_SUFFIX
        self._staged.append((tmp, path))
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
        staged, self._staged = self._staged, []
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
            _sync_dir(directory)
        for _, old in done:
            if old:
                _unlink(old)

    def abort(self) -> None:
        for tmp, _ in self._staged:
            _unlink(tmp)
        self._staged = []


def _roll_back(done: Sequence[Tuple[str, str]], staged: Sequence[Tuple[str, str]],
               exc: BaseException) -> None:
    """Undo the swaps already made. Raise CommitError rather than hide a failed undo."""
    stuck: List[str] = []
    for path, old in reversed(list(done)):
        try:
            if old:
                os.replace(old, path)
            elif os.path.lexists(path):
                os.unlink(path)               # nothing was there before this commit
        except OSError as e:
            stuck.append(f"{old or path} ({e.strerror})")
    if stuck:
        left = [t for t, _ in staged if os.path.lexists(t)]
        raise CommitError(
            "reboundary: a file swap failed and the originals could not all be put back. "
            "Nothing has been deleted. Each '" + OLD_SUFFIX + "' file is an original and "
            "each '" + NEW_SUFFIX + "' file is a rewrite: move every '" + OLD_SUFFIX +
            "' file back over the path it is named after, delete the '" + NEW_SUFFIX +
            "' files, then re-run. Could not put back: " + ", ".join(stuck)
            + (". Rewrites still on disk: " + ", ".join(left) if left else "")) from exc
    for tmp, _ in staged:
        _unlink(tmp)


def atomic_write_text(path: str, text: str) -> None:
    c = Commit()
    with c.open(path) as fh:
        fh.write(text)
    c.commit()


def check_keys_stay_unique(tables: Sequence[CleanTable], old2new: Dict[str, str]) -> None:
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
            if k is None:
                continue
            new = old2new.get(k, k)
            if new in claimed:
                raise ValueError(f"reboundary: {row[0]} and {claimed[new]} would both be "
                                 f"written with the same key {new}; nothing was written")
            claimed[new] = row[0]


def rewrite(tables: Sequence[CleanTable], accepted: Dict[str, Accepted], letters: Sequence[str],
            commit: Commit, wrap: int = 60) -> None:
    """Stage one genome's clean tables and FASTAs with `accepted` (old key -> result) applied."""
    old2new = {k: element_key(a.new_name) for k, a in accepted.items()}
    check_keys_stay_unique(tables, old2new)
    jobs: Dict[str, List[Tuple[List[Tuple[int, int]], str]]] = defaultdict(list)
    for a in accepted.values():
        m = a.member
        letter = letters[m.depth] if m.depth < len(letters) else "X"
        for host in hosts_of(m.nest_status):
            jobs[host].append((_added(m, a.left, a.right), letter))
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
            if k in accepted:
                row = list(accepted[k].fields) + row[len(COLUMNS):]
            else:
                row = list(row)
            row[i_nest] = rekey_nest(row[i_nest], old2new)
            rows.append(row)
        if os.path.isfile(t.fasta):
            with commit.open(t.fasta) as out:
                for name, seq in iter_fasta(t.fasta):
                    k = element_key(name)
                    if k in accepted:
                        name, seq = accepted[k].new_name, accepted[k].record
                    for spans, letter in jobs.get(k, ()):
                        start, orientation = where[k]
                        seq = paint(seq, orientation, start, spans, letter)
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
        if cols.get(row, "decision") != "extended":
            continue
        new = cols.get(row, "new_seq_id")
        key = element_key(new)
        if key is not None:
            out[key] = Rebound(cols.get(row, "old_seq_id"), new,
                               as_int(cols.get(row, "ext5"), 0), as_int(cols.get(row, "ext3"), 0))
    return out

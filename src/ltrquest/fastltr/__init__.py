"""FastLTR, the third LTR caller: build it, and run it on one chunk.

FastLTR (Gene Myers, FASTAN; the detector behind HiRepEat) counts 8-mer matches
along every diagonal of a 32 kb window, aligns only the diagonals that have
enough, and reports the pairs of non-overlapping direct repeats it finds. With
the tandem marking below it crosses a satellite, telomere or rDNA array once
rather than once per spacing, so it finishes a 500 kb chunk in seconds where
LTRharvest and LTR_FINDER can run for tens of minutes -- the chunks
``--ltr-timeout`` cuts short, which then contribute nothing, because neither
tool writes output before it finishes.

LTRquest ships its own adaptation (``FastLTR.c`` in this directory; the changes
are listed at its top): the repeat spacing it searches matches LTRharvest's
``-mindistltr``, a tandem alignment marks the bases it covers so other spacings
do not re-align through the same array, seeds are counted over neighbouring
diagonals so indels between old LTRs do not hide them, and hits are printed as
text. It is compiled against a pinned FASTAN checkout in ``--tools-dir``.

FastLTR's sequence reader splits a contig at every N and reads any other IUPAC
code as A, so masked bases -- genes, TRF arrays, earlier rounds' elements -- are
deleted before the call and the hits lifted back. Deleting a masked inner is
also what brings its outer's two LTRs back within FastLTR's 20 kb reach, the
same collapse HiRepEat performs between its rounds.

FastLTR reports every non-overlapping direct-repeat pair, where LTRharvest keeps
only the best of overlapping predictions, so inside a degenerate tandem array
(spacers between the units) it would pair each unit with its neighbour. An
LTR-RT has two copies; an array has more, at a near-constant spacing. A pair is
therefore dropped when either copy recurs one spacing further out (`periodic`).
"""

from __future__ import annotations

import bisect
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

FASTAN_URL = "https://github.com/thegenemyers/FASTAN"
FASTAN_REF = "843d0ed28b93b05ad428995d986312f02468afcd"
SOURCE = Path(__file__).resolve().parent / "FastLTR.c"
# Printed by the adapted binary's usage, and by no other FastLTR: a stock build
# on PATH has none of -x, -d or -w and must not be picked up.
MARKER = "LTRquest adaptation of FastLTR"
# FASTAN's own sources the adapted FastLTR.c is linked against.
LINK_SOURCES = ["FastLTR.c", "alncode.c", "align.c", "ANO.c", "GDB.c", "gene_core.c", "ONElib.c"]

DEFAULT_ARGS = "-d100 -D20000 -w8"
DEFAULT_MIN_LTR = 100
_ACGT = frozenset("ACGTacgt")

# Array test (`periodic`): a copy recurs one spacing further out when a local
# alignment covers >= PERIODIC_COVER of it at >= the pair's own identity minus
# PERIODIC_SLACK, within +-(0.2*spacing + 50) of the expected place (FastLTR's
# own band is 0.2*spacing). At most PERIODIC_QUERY bp of a copy are aligned.
# Heuristics, set on A. thaliana and checked on rice: they flag 1.0% and 1.4% of
# the FastLTR pairs behind high-confidence calls, 60-82% of those behind calls
# inside tandem arrays.
PERIODIC_COVER = 0.5
PERIODIC_SLACK = 0.10
PERIODIC_QUERY = 1000


# -----------------------------
# build
# -----------------------------

def tool_usable(bin_path: Path) -> bool:
    """True when `bin_path` runs and is the adapted build."""
    bin_path = Path(bin_path)
    if not bin_path.is_file() or not os.access(str(bin_path), os.X_OK):
        return False
    try:
        r = subprocess.run([str(bin_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return False
    return MARKER in (r.stdout or "") + (r.stderr or "")


def _compilers() -> List[str]:
    out = []
    for cc in (os.environ.get("CC"), "gcc", "cc"):
        if cc and cc not in out and shutil.which(shlex.split(cc)[0]):
            out.append(cc)
    return out


def build(src_dir: Path, out_bin: Path) -> str:
    """Compile the adapted FastLTR.c inside a FASTAN checkout.

    Returns "" on success, else the last compiler error. $CFLAGS / $CPPFLAGS /
    $LDFLAGS are passed through: an activated conda compiler puts the
    environment's zlib headers and libraries there.
    """
    src_dir, out_bin = Path(src_dir).resolve(), Path(out_bin).resolve()
    shutil.copyfile(SOURCE, src_dir / "FastLTR.c")
    env_flags = [f for v in ("CFLAGS", "CPPFLAGS", "LDFLAGS")
                 for f in shlex.split(os.environ.get(v, ""))]
    err = "no C compiler found ($CC, gcc or cc)"
    for cc in _compilers():
        cmd = shlex.split(cc) + ["-O3", "-fno-strict-aliasing"] + env_flags \
            + ["-o", str(out_bin)] + LINK_SOURCES + ["-lm", "-lz", "-lpthread"]
        r = subprocess.run(cmd, cwd=str(src_dir), stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True)
        if r.returncode == 0 and tool_usable(out_bin):
            return ""
        err = (r.stderr or "").strip()[-2000:] or f"{cc} exited {r.returncode}"
    return err


def ensure(tools_dir: Path) -> str:
    """Return a usable adapted FastLTR, building it in --tools-dir if need be.

    Same order as the other helpers: a previous build in --tools-dir (the
    container pre-populates it), then PATH, then clone FASTAN at the pinned
    commit and compile against it.
    """
    # Absolute: the compiler runs inside the checkout, and the driver's default
    # --tools-dir is relative.
    tools_dir = Path(tools_dir).resolve()
    tools_dir.mkdir(parents=True, exist_ok=True)
    src_dir = tools_dir / "FASTAN"
    local = src_dir / "FastLTR"
    if tool_usable(local):
        return str(local)

    on_path = shutil.which("FastLTR")
    if on_path and tool_usable(Path(on_path)):
        return on_path

    if not (src_dir / "GDB.c").is_file():
        if shutil.which("git") is None:
            raise RuntimeError(
                "FastLTR is not built in --tools-dir and cannot be built here: git is not "
                "installed.\n    conda install -c conda-forge git make c-compiler zlib\n"
                "Or run without it: --no-fastltr")
        if src_dir.exists():
            shutil.rmtree(src_dir)
        subprocess.run(["git", "clone", "--quiet", FASTAN_URL, str(src_dir)], check=True)
        subprocess.run(["git", "-C", str(src_dir), "checkout", "--quiet", FASTAN_REF],
                       check=True)

    err = build(src_dir, local)
    if not err:
        return str(local)
    raise RuntimeError(
        f"FastLTR build failed in {src_dir} (needs a C compiler and zlib):\n{err}\n"
        "    conda install -c conda-forge c-compiler zlib\n"
        "Or run without it: --no-fastltr")


# -----------------------------
# masked bases
# -----------------------------

@dataclass
class Piece:
    """One stretch of a record FastLTR sees as a single sequence.

    `orig` and `comp` hold, for each kept run of ACGT, its start in the record
    and in the piece; runs are adjacent in the piece and apart in the record
    wherever masked bases were deleted between them.
    """
    record: int
    seq: str
    orig: List[int] = field(default_factory=list)
    comp: List[int] = field(default_factory=list)
    lens: List[int] = field(default_factory=list)

    def lift(self, x: int) -> int:
        """The record position of piece position `x` (0-based, a kept base)."""
        k = bisect.bisect_right(self.comp, x) - 1
        if k < 0 or x - self.comp[k] >= self.lens[k]:
            raise ValueError(f"position {x} is outside piece of length {len(self.seq)}")
        return self.orig[k] + (x - self.comp[k])


def compress(seq: str, record: int = 0, split_char: Optional[str] = None,
             min_piece: int = 0) -> List[Piece]:
    """Delete masked (non-ACGT) bases from `seq`; split it at runs of `split_char`.

    A run of `split_char` (the far-mask letter of rounds 2+, whose candidates are
    excluded anyway) ends one piece and starts the next, so no hit spans it.
    Pieces shorter than `min_piece` are dropped: they cannot hold a repeat pair.
    """
    split = split_char.upper() if split_char else None
    pieces: List[Piece] = []
    cur = Piece(record, "")
    parts: List[str] = []
    clen = 0
    n = len(seq)
    i = 0
    while i < n:
        if seq[i] in _ACGT:
            j = i
            while j < n and seq[j] in _ACGT:
                j += 1
            cur.orig.append(i)
            cur.comp.append(clen)
            cur.lens.append(j - i)
            parts.append(seq[i:j])
            clen += j - i
            i = j
            continue
        j = i
        while j < n and seq[j] not in _ACGT:
            j += 1
        if split is not None and split in seq[i:j].upper():
            if clen:
                cur.seq = "".join(parts).upper()
                pieces.append(cur)
            cur, parts, clen = Piece(record, ""), [], 0
        i = j
    if clen:
        cur.seq = "".join(parts).upper()
        pieces.append(cur)
    return [p for p in pieces if len(p.seq) >= max(min_piece, 8)]


# -----------------------------
# one chunk
# -----------------------------

def _records(path: str) -> Iterable[Tuple[str, str]]:
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def _aligned(query: str, target: str) -> Tuple[float, float]:
    """(identity, fraction of `query` aligned) of the best local alignment."""
    if len(query) < 20 or len(target) < 20:
        return 0.0, 0.0
    import parasail
    tb = parasail.sw_trace_striped_16(query, target, 5, 2, _matrix()).traceback
    cols = len(tb.query)
    same = sum(1 for x, y in zip(tb.query, tb.ref) if x == y and x != "-")
    return same / max(1, cols), (cols - tb.query.count("-")) / len(query)


_MATRIX = None


def _matrix():
    global _MATRIX
    if _MATRIX is None:
        import parasail
        _MATRIX = parasail.matrix_create("ACGT", 2, -3)
    return _MATRIX


def periodic(seq: str, bs: int, be: int, as_: int, ae: int, identity: float) -> bool:
    """True when a pair's copy recurs one spacing further out: a tandem array.

    `seq` holds the pair at [bs, be) and [as_, ae), 0-based; `identity` is the
    pair's own. The 5' copy is looked for one spacing upstream, the 3' copy one
    spacing downstream.
    """
    d = as_ - bs
    w = int(0.2 * d) + 50
    for cs, ce, shift in ((bs, be, -d), (as_, ae, d)):
        if ce - cs > PERIODIC_QUERY:                 # its middle stands for it
            cs += (ce - cs - PERIODIC_QUERY) // 2
            ce = cs + PERIODIC_QUERY
        lo, hi = cs + shift - w, ce + shift + w
        if hi <= 0 or lo >= len(seq):
            continue
        ident, cover = _aligned(seq[cs:ce], seq[max(0, lo):min(len(seq), hi)])
        if cover >= PERIODIC_COVER and ident >= identity - PERIODIC_SLACK:
            return True
    return False


def scn_line(piece: Piece, bs: int, be: int, as_: int, ae: int, diffs: int) -> str:
    """One LTRharvest-style SCN row (1-based, inclusive) for a hit in `piece`."""
    l5s, l5e = piece.lift(bs) + 1, piece.lift(be - 1) + 1
    l3s, l3e = piece.lift(as_) + 1, piece.lift(ae - 1) + 1
    span = max(be - bs, ae - as_)
    sim = max(0.0, min(100.0, 100.0 * (1.0 - diffs / span))) if span else 0.0
    return (f"{l5s} {l3e} {l3e - l5s + 1} "
            f"{l5s} {l5e} {l5e - l5s + 1} "
            f"{l3s} {l3e} {l3e - l3s + 1} "
            f"{sim:.2f} {piece.record}")


def hits_to_scn(text: str, pieces: List[Piece], min_ltr: int,
                stats: Optional[dict] = None) -> List[str]:
    """Turn FastLTR -x lines into SCN rows: pairs whose copies are >= min_ltr bp
    and that are not units of a tandem array (`periodic`). `stats`, if given,
    gets the number of pairs dropped as array units under "periodic"."""
    rows = []
    n_periodic = 0
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 6 or not f[0].startswith("p"):
            continue
        try:
            k = int(f[0][1:])
            bs, be, as_, ae, diffs = (int(x) for x in f[1:6])
        except ValueError:
            continue
        if not (0 <= k < len(pieces)) or be - bs < min_ltr or ae - as_ < min_ltr:
            continue
        identity = 1.0 - diffs / max(be - bs, ae - as_)
        if periodic(pieces[k].seq, bs, be, as_, ae, identity):
            n_periodic += 1
            continue
        rows.append(scn_line(pieces[k], bs, be, as_, ae, diffs))
    if stats is not None:
        stats["periodic"] = stats.get("periodic", 0) + n_periodic
    return rows


def run_chunk(fastltr_bin: str, chunk_fa: str, wdir: str, scn_path: str,
              args: List[str], min_ltr: int = DEFAULT_MIN_LTR, timeout_s: int = 0,
              split_char: Optional[str] = None, verbose: bool = False,
              stats: Optional[dict] = None) -> int:
    """Run FastLTR on one chunk FASTA and write its hits as SCN; returns the row count.

    Like the other two callers, a run past `timeout_s` is killed; FastLTR
    prints each hit as it is found, so what it found by then is kept. `stats`
    collects the pairs dropped as tandem-array units (see `hits_to_scn`).
    """
    wdir_p = Path(wdir)
    pieces: List[Piece] = []
    for rec, (_name, seq) in enumerate(_records(chunk_fa)):
        pieces += compress(seq, rec, split_char, min_piece=2 * min_ltr)
    with open(scn_path, "w") as out:
        out.write("## FastLTR (LTRquest adaptation)\n")
        out.write("## s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr\n")
    if not pieces:
        return 0

    comp_fa = wdir_p / "fastltr.in.fa"
    with open(comp_fa, "w") as fh:
        for k, p in enumerate(pieces):
            fh.write(f">p{k}\n")
            for i in range(0, len(p.seq), 80):
                fh.write(p.seq[i:i + 80] + "\n")

    cmd = [fastltr_bin, "-x", "-T1"] + list(args) + [comp_fa.name, "fastltr.1aln"]
    if verbose:
        print(f"  $ {' '.join(cmd)}  (cwd: {wdir})", file=sys.stderr)
    try:
        r = subprocess.run(cmd, cwd=str(wdir_p), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           text=True, timeout=(timeout_s if timeout_s and timeout_s > 0 else None))
        if r.returncode != 0:
            raise RuntimeError(f"FastLTR failed on {chunk_fa}:\n{(r.stderr or '').strip()}")
        text = r.stdout or ""
    except subprocess.TimeoutExpired as te:
        text = te.stdout or ""
        if isinstance(text, bytes):
            text = text.decode(errors="replace")
        print(f"[WARN] FastLTR timeout after {timeout_s}s on {chunk_fa}; kept "
              f"{text.count(chr(10))} hit(s) found before it.", file=sys.stderr)

    rows = hits_to_scn(text, pieces, min_ltr, stats)
    with open(scn_path, "a") as out:
        for row in rows:
            out.write(row + "\n")
    return len(rows)

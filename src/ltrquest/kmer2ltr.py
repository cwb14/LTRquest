"""Locating, invoking and reading Kmer2LTR.

Kmer2LTR is an installable package with a console script, but LTRquest has
always been able to run it from a bare `git clone`, and on a shared cluster
that matters: `pip install` into the active environment is often not the
caller's to make. So the clone is put on PYTHONPATH and the package is run as
a module. Its own dependencies -- parasail, pywfa, numpy -- come from
LTRquest's environment either way.

The output TSV's column order is the field order of a dataclass upstream, so
it is read by header rather than by position.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterator, Optional, Sequence

from .table import parse_header

REPO_URL = "https://github.com/cwb14/Kmer2LTR.git"

COLUMNS = [
    "seq_id", "seq_len", "status",
    "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end",
    "ltr5_len", "ltr3_len", "flank5_len", "flank3_len",
    "aln_len", "n_sites", "n_ts", "n_tv", "n_gapcols",
    "identity", "p_dist", "k2p", "k2p_se",
    "bitscore", "flank_margin_bits", "cigar", "motif", "k2p_time",
    "orientation", "tsd", "tsd_offset", "tsd_input",
]


def _module_argv(src_dir: Path) -> list[str]:
    return [f"PYTHONPATH={src_dir}", sys.executable, "-m", "kmer2ltr"]


def _importable(src_dir: Path) -> bool:
    return (src_dir / "kmer2ltr" / "cli.py").is_file()


def resolve(tools_dir: Path) -> list[str]:
    """An argv prefix that runs Kmer2LTR, cloning it if that is the only way.

    In the clone-fallback case the prefix leads with a `NAME=VALUE`
    environment assignment (see `split_env`) rather than an executable --
    the same shape `env VAR=val cmd...` takes in a shell. `run` unpacks this
    for you; a caller that hands the result straight to `subprocess` instead
    must call `split_env` on it first, or run it through a shell with `env`
    in front.
    """
    exe = shutil.which("Kmer2LTR")
    if exe:
        return [exe]

    tools_dir = Path(tools_dir)
    src_dir = tools_dir / "Kmer2LTR" / "src"
    if _importable(src_dir):
        return _module_argv(src_dir)

    tools_dir.mkdir(parents=True, exist_ok=True)
    target = tools_dir / "Kmer2LTR"
    if not target.exists():
        try:
            subprocess.run(["git", "clone", "--depth", "1", REPO_URL, str(target)],
                           check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"failed to clone {REPO_URL} into {target}: {e}. Compute nodes "
                f"on this cluster often have no outbound network access -- "
                f"clone it on a login node instead, or install Kmer2LTR "
                f"(pip install kmer2ltr) so it is on PATH. If {target} exists "
                f"from this failed attempt, delete it first: a leftover "
                f"directory there stops the next run from retrying the clone."
            ) from e
    if _importable(src_dir):
        return _module_argv(src_dir)

    raise RuntimeError(
        f"Kmer2LTR is neither on PATH nor importable from {src_dir}. "
        f"Install it (pip install kmer2ltr) or point --tools-dir at a clone."
    )


def split_env(prefix: Sequence[str]) -> tuple[list[str], dict[str, str]]:
    """Separate leading NAME=VALUE assignments from the command itself."""
    env = dict(os.environ)
    argv = list(prefix)
    while argv and "=" in argv[0] and not Path(argv[0]).exists():
        name, _, value = argv[0].partition("=")
        existing = env.get(name)
        env[name] = f"{value}{os.pathsep}{existing}" if existing else value
        argv = argv[1:]
    return argv, env


def build_argv(prefix: Sequence[str], in_fa, out_tsv, *,
               threads: int = 1,
               mutation_rate: Optional[float] = None,
               genome: Optional[str] = None,
               trim_flanks: bool = False,
               ltr_cluster: bool = False,
               internal_cluster: bool = False,
               min_seq_id: Optional[float] = None) -> list[str]:
    argv = list(prefix) + [str(in_fa), "-o", str(out_tsv), "-t", str(threads)]
    if mutation_rate is not None:
        argv += ["-u", f"{mutation_rate:g}"]
    if genome:
        argv += ["--genome", str(genome)]
    if trim_flanks:
        argv += ["--trim-flanks"]
    if ltr_cluster:
        argv += ["--ltr-cluster"]
    if internal_cluster:
        argv += ["--internal-cluster"]
    if min_seq_id is not None:
        argv += ["--min-seq-id", f"{min_seq_id:g}"]
    return argv


def run(prefix: Sequence[str], in_fa, out_tsv, *, verbose: bool = False, **kwargs) -> Path:
    argv = build_argv(prefix, in_fa, out_tsv, **kwargs)
    argv, env = split_env(argv)
    if verbose:
        argv = argv + ["-v"]
        print("+ " + " ".join(argv), flush=True)
    subprocess.run(argv, check=True, env=env)
    return Path(out_tsv)


def trimmed_path(out_tsv) -> Path:
    return Path(str(Path(out_tsv).with_suffix("")) + ".trimmed.fa")


def read_rows(tsv) -> Iterator[dict[str, str]]:
    with open(tsv) as fh:
        first = fh.readline()
        names = parse_header("#" + first if not first.startswith("#") else first)
        if not names or names[0] != "seq_id":
            raise ValueError(f"{tsv}: expected a Kmer2LTR header line, got {first!r}")
        for line in fh:
            if not line.strip():
                continue
            yield dict(zip(names, line.rstrip("\n").split("\t")))


def status_counts(tsv) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in read_rows(tsv):
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    return counts

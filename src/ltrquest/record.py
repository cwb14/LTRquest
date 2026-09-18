"""Per-genome detection records, so a later run can reuse a genome's detection.

Detection is the expensive, per-genome part of a run; the pooled stages after it
-- clustering, the false-positive call, annotation -- are cheap by comparison.
Once a genome's detection is final, the driver writes <prefix>.detect.json beside
its outputs: what went in (genome and protein checksums, every setting that
shapes detection) and what came out. A later run over the same directory,
typically with more genomes, reuses that genome rather than detecting it again
when nothing that went in has changed, and refuses when something has, so one
pool never mixes genomes detected under different settings.

Usage:
  python -m ltrquest.record check   --prefix P --genome G [--proteins F] [--setting K=V ...]
  python -m ltrquest.record write   --prefix P --genome G [--proteins F] [--setting K=V ...]
                                    [--fp-masked]
  python -m ltrquest.record outputs --prefix P
  python -m ltrquest.record clear   --prefix P

check prints 'reuse' or 'detect'. A record that no longer matches is reported
difference by difference on stderr, with exit status 1.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

from . import __version__

FORMAT = 1
SUFFIX = ".detect.json"
_CHUNK = 1 << 20

# Names detection writes for a genome, after its '<prefix>_'.
_OUTPUTS = r"r\d+(_ltr\.(tsv|fa)|\.work)|depth\d+_ltr\.(tsv|fa)|strand_recovery\.tsv"
_CLEAN = r"depth\d+_clean_ltr\.(tsv|fa)"


def record_path(indir: str, prefix: str) -> str:
    return os.path.join(indir, prefix + SUFFIX)


def fingerprint(path: Optional[str]) -> Optional[dict]:
    """Name, size and SHA-256 of a file's bytes as stored (a .gz is hashed compressed)."""
    if not path:
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            digest.update(chunk)
    return {"name": os.path.basename(path), "size": os.path.getsize(path),
            "sha256": digest.hexdigest()}


def detection_outputs(indir: str, prefix: str) -> List[str]:
    """What detection leaves for the pooled stages to read, and nothing they rewrite.

    The per-round tables and FASTAs, the round work directories (TEsorter2 calls,
    pass-2 alignments, the miniprot GFF), the reconciled depth tables, and the
    strand-recovery sidecar when there is one. The _clean_ tables, GFF3s, logs
    and plots are regenerated from these on every run.
    """
    return _named(indir, prefix, _OUTPUTS)


def _named(indir: str, prefix: str, alternatives: str) -> List[str]:
    pattern = re.compile(rf"^{re.escape(prefix)}_({alternatives})$")
    return sorted(n for n in os.listdir(indir) if pattern.match(n))


def clear(indir: str, prefix: str) -> List[str]:
    """Remove what an earlier detection of this genome left, record included.

    A genome detected again replaces its earlier outputs wholesale: a depth
    table or round from a longer earlier run would otherwise survive beside the
    new ones, and be recorded and pooled as this genome's. Only names detection
    and the FP stage write are touched, never anything else sharing the prefix.
    """
    removed = _named(indir, prefix, f"{_OUTPUTS}|{_CLEAN}")
    if os.path.isfile(record_path(indir, prefix)):
        removed.append(prefix + SUFFIX)
    for name in removed:
        path = os.path.join(indir, name)
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        else:
            os.unlink(path)
    return removed


def build(indir: str, prefix: str, genome: str, proteins: Optional[str],
          settings: Dict[str, str], fp_masked: bool = False) -> dict:
    return {
        "format": FORMAT,
        "ltrquest": __version__,
        "prefix": prefix,
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "genome": fingerprint(genome),
        "proteins": fingerprint(proteins),
        "settings": dict(sorted(settings.items())),
        # Detected on the FP-masked genome, not the one named above.
        "fp_masked": fp_masked,
        "outputs": detection_outputs(indir, prefix),
    }


def write(indir: str, prefix: str, genome: str, proteins: Optional[str],
          settings: Dict[str, str], fp_masked: bool = False) -> str:
    """Write the record through a temporary file, so a killed run leaves none."""
    rec = build(indir, prefix, genome, proteins, settings, fp_masked)
    if not rec["outputs"]:
        raise ValueError(f"no detection outputs for {prefix} in {indir}")
    path = record_path(indir, prefix)
    fd, tmp = tempfile.mkstemp(prefix=f".{prefix}.", suffix=".tmp", dir=indir)
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(rec, fh, indent=1)
            fh.write("\n")
        # mkstemp creates 0600; the record should read like any other output.
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return path


def _same(a: Optional[str], b: Optional[str]) -> bool:
    """Equal as written, or as numbers: --fp-mask-threshold 0.1 is the default 0.10."""
    if a == b:
        return True
    try:
        return a is not None and b is not None and float(a) == float(b)
    except ValueError:
        return False


def _show(value: Optional[str]) -> str:
    return "(unset)" if value is None else (value if value != "" else '""')


def _file_change(label: str, saved: Optional[dict], path: Optional[str]) -> Optional[str]:
    if saved is None and not path:
        return None
    if saved is None:
        return f"{label}: none before, {os.path.basename(path)} now"
    if not path:
        return f"{label}: {saved.get('name')} before, none now"
    # Size first: a different size settles it without reading a large genome.
    if (os.path.getsize(path) != saved.get("size")
            or fingerprint(path)["sha256"] != saved.get("sha256")):
        return f"{label}: {os.path.basename(path)} differs from the {saved.get('name')} it was detected from"
    return None


def check(indir: str, prefix: str, genome: str, proteins: Optional[str],
          settings: Dict[str, str]) -> Tuple[str, List[str], List[str]]:
    """('detect' | 'reuse' | 'refuse', problems, warnings) for one genome."""
    path = record_path(indir, prefix)
    if not os.path.isfile(path):
        return "detect", [], []
    try:
        with open(path) as fh:
            saved = json.load(fh)
    except (OSError, ValueError) as exc:
        return "refuse", [f"unreadable record {path}: {exc}"], []
    if not isinstance(saved, dict) or saved.get("format") != FORMAT:
        return "refuse", [f"{path} is not a format-{FORMAT} detection record"], []

    problems: List[str] = []
    for label, now in (("genome", genome), ("proteins", proteins)):
        change = _file_change(label, saved.get(label), now)
        if change:
            problems.append(change)

    before = saved.get("settings") or {}
    for key in sorted(set(before) | set(settings)):
        if not _same(before.get(key), settings.get(key)):
            problems.append(f"{key}: {_show(before.get(key))} before, "
                            f"{_show(settings.get(key))} now")

    outputs = saved.get("outputs") or []
    if not outputs:
        problems.append("the record lists no detection outputs")
    problems += [f"missing output: {name}" for name in outputs
                 if not os.path.exists(os.path.join(indir, name))]

    warnings = []
    if saved.get("ltrquest") != __version__:
        warnings.append(f"detected with LTRquest {saved.get('ltrquest')}; this is {__version__}")
    if saved.get("fp_masked"):
        warnings.append("detected on an FP-masked genome, so the pool no longer holds its "
                        "false positives; a genome added now may stay unmasked where a "
                        "run from scratch would mask it (see the FP fraction in "
                        "*_fpcheck.log)")
    return ("refuse" if problems else "reuse"), problems, warnings


def _settings(pairs: List[str], parser: argparse.ArgumentParser) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for pair in pairs:
        key, sep, value = pair.partition("=")
        if not sep or not key:
            parser.error(f"--setting expects KEY=VALUE, got {pair!r}")
        out[key] = value
    return out


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="ltrquest-record",
        description="Write or check the per-genome detection record that lets a later run "
                    "reuse a genome instead of detecting it again.")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, text in (("check", "print 'reuse' or 'detect'; exit 1 if the record no "
                                 "longer matches"),
                       ("write", "record a finished detection"),
                       ("outputs", "list the detection outputs a record names"),
                       ("clear", "remove an earlier detection's outputs and record")):
        cmd = sub.add_parser(name, help=text, description=text)
        cmd.add_argument("--indir", default=".", help="run directory (default: .)")
        cmd.add_argument("--prefix", required=True, help="genome output prefix, e.g. Athal_LTRs")
        if name in ("outputs", "clear"):
            continue
        cmd.add_argument("--genome", required=True, help="the genome FASTA as given to ltrquest")
        cmd.add_argument("--proteins", default=None, help="the --proteins FASTA, if any")
        cmd.add_argument("--setting", action="append", default=[], metavar="KEY=VALUE",
                         help="one detection setting; repeat for each")
        if name == "write":
            cmd.add_argument("--fp-masked", action="store_true",
                             help="detection ran on the FP-masked genome")
    args = parser.parse_args(argv)

    if args.command == "outputs":
        with open(record_path(args.indir, args.prefix)) as fh:
            print("\n".join(json.load(fh)["outputs"]))
        return 0
    if args.command == "clear":
        removed = clear(args.indir, args.prefix)
        if removed:
            print("\n".join(removed))
        return 0

    settings = _settings(args.setting, parser)
    if args.command == "write":
        try:
            write(args.indir, args.prefix, args.genome, args.proteins, settings,
                  args.fp_masked)
        except ValueError as exc:
            print(f"[ltrquest.record] ERROR: {exc}", file=sys.stderr)
            return 1
        return 0

    verdict, problems, warnings = check(args.indir, args.prefix, args.genome,
                                        args.proteins, settings)
    for message in warnings:
        print(f"[ltrquest.record] WARNING: {args.prefix}: {message}", file=sys.stderr)
    if verdict == "refuse":
        for message in problems:
            print(f"  {args.prefix}: {message}", file=sys.stderr)
        return 1
    print(verdict)
    return 0


if __name__ == "__main__":
    sys.exit(main())

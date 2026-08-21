"""Shared fixtures.

The synthetic fixtures are deliberately tiny: every unit test here runs without
GenomeTools, LTR_finder, MMseqs2, HMMER or any of the other external binaries the
pipeline drives, so the whole non-slow suite finishes in about a second and works
on a bare CI runner. The one piece of real data, `tests/data/`, is used only for
streaming and parsing checks.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pytest

import ltrquest

DATA = Path(__file__).parent / "data"
SCRIPTS = Path(ltrquest.__file__).resolve().parent / "scripts"


@pytest.fixture(scope="session")
def driver() -> Path:
    """The packaged bash driver, wherever the install put it."""
    return SCRIPTS / "ltrquest.sh"


@pytest.fixture(scope="session")
def driver_iupac_seq(driver: Path) -> list[str]:
    """The per-round mask letters the driver hands to the detector."""
    match = re.search(r"^IUPAC_SEQ=\(([^)]*)\)", driver.read_text(), re.MULTILINE)
    assert match, "IUPAC_SEQ not found in the driver; the tests below are stale"
    return match.group(1).split()


@pytest.fixture(scope="session")
def athal_genome() -> Path:
    """The shipped Arabidopsis thaliana chr2, gzipped."""
    path = DATA / "Athal_tair10_chr2.fa.gz"
    if not path.is_file():
        pytest.skip(f"test genome not present: {path}")
    return path


@pytest.fixture(scope="session")
def athal_slice(tmp_path_factory, athal_genome: Path) -> Path:
    """A 200 kb plain-FASTA slice of the real chromosome.

    Real sequence, small enough to be free. Used where a test wants genuine base
    composition and line wrapping rather than a hand-written toy contig.
    """
    out = tmp_path_factory.mktemp("athal") / "chr2_slice.fa"
    want = 200_000
    seq: list[str] = []
    total = 0
    with gzip.open(athal_genome, "rt") as handle:
        name = handle.readline().split()[0].lstrip(">")
        for line in handle:
            if line.startswith(">"):
                break
            line = line.strip()
            seq.append(line)
            total += len(line)
            if total >= want:
                break
    body = "".join(seq)[:want]
    with out.open("w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(body), 60):
            fh.write(body[i : i + 60] + "\n")
    return out


@pytest.fixture
def toy_genome(tmp_path: Path) -> Path:
    """Two short contigs with a recognisable base at every position.

    chr1 is 300 bp of A, chr2 is 120 bp of C, so any masking is obvious by eye
    and by assertion.
    """
    path = tmp_path / "toy.fa"
    path.write_text(">chr1\n" + "A" * 300 + "\n>chr2\n" + "C" * 120 + "\n")
    return path

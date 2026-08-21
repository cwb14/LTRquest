"""End-to-end run on the shipped Arabidopsis chromosome.

Marked `slow` and deselected in CI: this drives GenomeTools, LTR_finder,
MMseqs2, HMMER, miniprot and the bundled Kmer2LTR/TEsorter2 helpers, and takes
minutes rather than milliseconds. It is the check that the stages still compose,
which no amount of unit testing can replace.

Run it against the container, where the whole toolchain is present::

    docker run --rm -v "$PWD:/w" -w /w ghcr.io/cwb14/ltrquest:1.0.0 \\
        pytest -m slow tests/test_end_to_end.py

or against a conda environment built from environment.yml::

    pytest -m slow tests/test_end_to_end.py

It skips itself, rather than failing, when the toolchain is not on PATH.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess

import pytest

pytestmark = pytest.mark.slow

# Everything the detector shells out to that conda provides. Kmer2LTR and
# TEsorter2 are fetched into --tools-dir on demand, so they are not listed.
REQUIRED_TOOLS = ["gt", "ltr_finder", "mmseqs", "hmmsearch", "blastn", "minimap2", "miniprot"]

# Two rounds, so masking and cross-round reconciliation both run. One round
# would exercise detection only.
ROUNDS = 2
SLICE_BP = 2_000_000


def _threads() -> str:
    return os.environ.get("LTRQUEST_E2E_THREADS", str(min(8, os.cpu_count() or 1)))


@pytest.fixture(scope="module")
def toolchain():
    if shutil.which("ltrquest") is None:
        pytest.skip("ltrquest is not on PATH; install the package first (pip install .)")
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        pytest.skip(f"external toolchain not on PATH: {', '.join(missing)}")


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory, toolchain, athal_genome):
    """A completed LTRquest run over a slice of the real chromosome."""
    work = tmp_path_factory.mktemp("e2e")

    genome = work / "athal_slice.fa"
    written = 0
    with gzip.open(athal_genome, "rt") as src, genome.open("w") as dst:
        dst.write(">chr_test\n")
        src.readline()
        for line in src:
            if line.startswith(">") or written >= SLICE_BP:
                break
            dst.write(line)
            written += len(line.strip())

    proteins = work / "prot_slice.fa"
    with gzip.open(athal_genome.parent / "Athal.pep.gz", "rt") as src, proteins.open("w") as dst:
        for i, line in enumerate(src):
            if i >= 40_000:
                break
            dst.write(line)

    result = subprocess.run(
        ["ltrquest",
         "--genome", genome.name,
         "--proteins", proteins.name,
         "--threads", _threads(),
         "--max-rounds", str(ROUNDS),
         "--terminate_count", "1",
         "--no-plots"],
        cwd=work, capture_output=True, text=True,
    )
    if result.returncode != 0:
        pytest.fail(f"ltrquest exited {result.returncode}\n--- stderr ---\n{result.stderr[-4000:]}")
    return work, result.stdout


def test_every_round_ran(run_dir):
    _work, stdout = run_dir
    for r in range(1, ROUNDS + 1):
        assert f"Round {r}: detected" in stdout, f"round {r} did not report a count"


def test_the_genome_was_masked_between_rounds(run_dir):
    work, stdout = run_dir
    assert "Masking original genome for next round" in stdout
    masked = list(work.glob("*_r1.fa"))
    assert masked, "no masked genome was written for round 2"
    seq = "".join(
        line.strip() for line in masked[0].read_text().splitlines() if not line.startswith(">")
    )
    # Round 1 paints its hits 'N' and everything far from one 'V'.
    assert "N" in seq and "V" in seq


def test_depth_buckets_exist_and_are_populated(run_dir):
    work, _stdout = run_dir
    depth0 = work / "athal_slice_LTRs_depth0_clean_ltr.tsv"
    assert depth0.is_file(), "no depth0 table"
    rows = [ln for ln in depth0.read_text().splitlines() if ln and not ln.startswith("#")]
    assert rows, "depth0 table is empty"


def test_table_and_fasta_agree(run_dir):
    work, _stdout = run_dir
    for tsv in sorted(work.glob("athal_slice_LTRs_depth*_clean_ltr.tsv")):
        fasta = tsv.with_name(tsv.name.replace(".tsv", ".fa"))
        assert fasta.is_file(), f"no FASTA beside {tsv.name}"
        names = {
            ln.split("\t")[0] for ln in tsv.read_text().splitlines()
            if ln and not ln.startswith("#")
        }
        headers = {
            ln[1:].strip() for ln in fasta.read_text().splitlines() if ln.startswith(">")
        }
        assert names == headers, f"{tsv.name} and its FASTA disagree on membership"


def test_nested_elements_are_masked_in_their_hosts_fasta(run_dir):
    work, _stdout = run_dir
    depth1 = work / "athal_slice_LTRs_depth1_clean_ltr.fa"
    if not depth1.is_file():
        pytest.skip("this slice produced no nested elements")
    seq = "".join(
        line.strip() for line in depth1.read_text().splitlines() if not line.startswith(">")
    )
    # A depth1 record carries its depth0 insert hard-masked as 'N'.
    assert "N" in seq, "depth1 record has no masked insert"


def test_nest_status_is_reciprocal(run_dir):
    work, _stdout = run_dir
    depth1 = work / "athal_slice_LTRs_depth1_clean_ltr.tsv"
    if not depth1.is_file():
        pytest.skip("this slice produced no nested elements")
    for line in depth1.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        status = line.split("\t")[-1]
        assert status.startswith("nest-outer:"), (
            f"a depth1 element must name the element inside it, got: {status}"
        )


def test_gff3_is_written_and_well_formed(run_dir):
    work, _stdout = run_dir
    gff3 = work / "athal_slice_LTRs_all_depth_LTR_cleaned.gff3"
    assert gff3.is_file(), "no pooled GFF3"

    features = 0
    for line in gff3.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        assert len(cols) == 9, f"GFF3 line has {len(cols)} columns: {line[:80]}"
        assert cols[1] == "LTRquest", f"unexpected source column: {cols[1]}"
        assert int(cols[3]) <= int(cols[4]), f"start > end: {line[:80]}"
        assert cols[6] in "+-.", f"bad strand: {cols[6]}"
        features += 1
    assert features, "GFF3 has no features"


def test_the_family_column_was_filled_in(run_dir):
    work, _stdout = run_dir
    depth0 = work / "athal_slice_LTRs_depth0_clean_ltr.tsv"
    header = depth0.read_text().splitlines()[0].lstrip("#").split("\t")
    assert "family" in header and "strand" in header
    idx = header.index("family")
    families = {
        ln.split("\t")[idx] for ln in depth0.read_text().splitlines()
        if ln and not ln.startswith("#")
    }
    assert any(f.startswith("athal_slice_LTRs_fam") for f in families), (
        f"no family labels were assigned: {sorted(families)[:5]}"
    )

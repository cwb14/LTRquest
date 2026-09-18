"""End-to-end run on the shipped Arabidopsis chromosome.

Marked `slow` and deselected in CI: this drives GenomeTools, LTR_finder,
MMseqs2, HMMER, miniprot and the bundled Kmer2LTR/TEsorter2 helpers, and takes
minutes rather than milliseconds. It is the check that the stages still compose,
which no amount of unit testing can replace.

Run it against the container, where the whole toolchain is present::

    docker run --rm -v "$PWD:/w" -w /w ghcr.io/cwb14/ltrquest:1.0.1 \\
        pytest -m slow tests/test_end_to_end.py

or against a conda environment built from environment.yml::

    pytest -m slow tests/test_end_to_end.py

It skips itself, rather than failing, when the toolchain is not on PATH.
"""

from __future__ import annotations

import gzip
import hashlib
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


# ---------------------------------------------------------------------------
# Adding a genome to an earlier multi-genome run
# ---------------------------------------------------------------------------

# Three 1 Mb windows of chr2, each holding enough LTR-RTs to survive one round.
WINDOWS = {"a": 2_000_000, "b": 3_500_000, "c": 5_000_000}
WINDOW_BP = 1_000_000


def _write_windows(athal_genome, dirs):
    end = max(WINDOWS.values()) + WINDOW_BP
    chunks, total = [], 0
    with gzip.open(athal_genome, "rt") as src:
        src.readline()
        for line in src:
            if line.startswith(">") or total >= end:
                break
            chunks.append(line.strip())
            total += len(chunks[-1])
    seq = "".join(chunks)
    with gzip.open(athal_genome.parent / "Athal.pep.gz", "rt") as src:
        proteins = "".join(line for _, line in zip(range(40_000), src))
    for d in dirs:
        for name, start in WINDOWS.items():
            window = seq[start:start + WINDOW_BP]
            body = "\n".join(window[i:i + 60] for i in range(0, len(window), 60))
            (d / f"{name}.fa").write_text(f">chr_{name}\n{body}\n")
        (d / "prot.fa").write_text(proteins)


def _ltrquest(cwd, *genomes, extra=()):
    return subprocess.run(
        ["ltrquest", "--genome", *genomes, "--proteins", "prot.fa",
         "--threads", _threads(), "--max-rounds", "1", "--terminate_count", "1",
         "--no-plots", *extra],
        cwd=cwd, capture_output=True, text=True,
    )


def _require(result):
    if result.returncode != 0:
        pytest.fail(f"ltrquest exited {result.returncode}\n--- stderr ---\n{result.stderr[-4000:]}")
    return result


def _detection_digests(work, prefixes):
    """sha256 of every round table and work-directory file of these genomes."""
    return {
        str(path.relative_to(work)): hashlib.sha256(path.read_bytes()).hexdigest()
        for prefix in prefixes
        for top in sorted(work.glob(f"{prefix}_r*"))
        for path in ([top] if top.is_file() else sorted(top.rglob("*")))
        if path.is_file()
    }


def _families(tsv):
    lines = [ln for ln in tsv.read_text().splitlines() if ln]
    idx = lines[0].lstrip("#").split("\t").index("family")
    return sorted((ln.split("\t")[0], ln.split("\t")[idx]) for ln in lines[1:])


def _gff3_families(gff3):
    found = []
    for line in gff3.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
        if "family" in attrs:
            found.append((cols[0], cols[2], cols[3], cols[4],
                          *(attrs.get(k) for k in ("family", "family_rep", "family_size",
                                                   "family_clades"))))
    return sorted(found)


@pytest.fixture(scope="module")
def incremental(tmp_path_factory, toolchain, athal_genome):
    """Genomes a+b, then a+b+c in the same directory; and a+b+c from scratch."""
    root = tmp_path_factory.mktemp("incremental")
    inc, fresh = root / "inc", root / "fresh"
    inc.mkdir()
    fresh.mkdir()
    _write_windows(athal_genome, (inc, fresh))

    _require(_ltrquest(inc, "a.fa", "b.fa"))
    before = _detection_digests(inc, ("a_LTRs", "b_LTRs"))
    second = _require(_ltrquest(inc, "a.fa", "b.fa", "c.fa"))
    _require(_ltrquest(fresh, "a.fa", "b.fa", "c.fa"))
    return inc, fresh, second.stdout, before


def test_adding_a_genome_detects_only_that_genome(incremental):
    _inc, _fresh, stdout, _before = incremental
    assert "reusing:     a_LTRs b_LTRs" in stdout
    assert "detecting:   c_LTRs" in stdout


def test_reused_detection_files_are_untouched(incremental):
    inc, _fresh, _stdout, before = incremental
    assert before, "the first run left no detection files"
    assert _detection_digests(inc, ("a_LTRs", "b_LTRs")) == before


def test_every_genome_has_a_detection_record(incremental):
    inc, _fresh, _stdout, _before = incremental
    for prefix in ("a_LTRs", "b_LTRs", "c_LTRs"):
        assert (inc / f"{prefix}.detect.json").is_file()


def test_families_match_a_run_from_scratch(incremental):
    inc, fresh, _stdout, _before = incremental
    tables = sorted(p.name for p in fresh.glob("*_depth*_clean_ltr.tsv"))
    assert tables
    assert sorted(p.name for p in inc.glob("*_depth*_clean_ltr.tsv")) == tables
    for name in tables:
        assert _families(inc / name) == _families(fresh / name), name


def test_gff3_family_attributes_match_a_run_from_scratch(incremental):
    inc, fresh, _stdout, _before = incremental
    for prefix in ("a_LTRs", "b_LTRs", "c_LTRs"):
        name = f"{prefix}_all_depth_LTR_cleaned.gff3"
        assert _gff3_families(inc / name) == _gff3_families(fresh / name), name


def test_a_changed_detection_setting_is_refused(incremental):
    inc, _fresh, _stdout, _before = incremental
    result = _ltrquest(inc, "a.fa", "b.fa", "c.fa", extra=("--max-rounds", "2"))
    assert result.returncode != 0
    assert "max_rounds: 1 before, 2 now" in result.stderr

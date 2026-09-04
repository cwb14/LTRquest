import subprocess
import textwrap

import pytest

from ltrquest import kmer2ltr

TSV = textwrap.dedent("""\
    seq_id\tseq_len\tstatus\tltr5_start\tltr3_end\tk2p\ttsd\ttsd_input
    chr1:100-200#LTR/Gypsy\t101\tpass\t1\t101\t0.012\tAAGCT\tAAGCT
    chr1:300-500#LTR/Copia\t201\tno_pair\tNA\tNA\tNA\tNA\tNA
    chr1:600-900#LTR/Gypsy\t301\tpass\t5\t298\t0.05\t.\t.
    """)


def test_trimmed_path_strips_one_suffix():
    assert kmer2ltr.trimmed_path("a/b_kmer2ltr.tsv").name == "b_kmer2ltr.trimmed.fa"
    assert kmer2ltr.trimmed_path("a/b_all_ltr").name == "b_all_ltr.trimmed.fa"


def test_columns_are_the_documented_twentynine():
    assert len(kmer2ltr.COLUMNS) == 29
    assert kmer2ltr.COLUMNS[0] == "seq_id"
    assert kmer2ltr.COLUMNS[-1] == "tsd_input"
    assert "flank5_len" in kmer2ltr.COLUMNS
    assert "orientation" in kmer2ltr.COLUMNS


def test_read_rows_yields_dicts(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text(TSV)
    rows = list(kmer2ltr.read_rows(p))
    assert len(rows) == 3
    assert rows[0]["seq_id"] == "chr1:100-200#LTR/Gypsy"
    assert rows[1]["status"] == "no_pair"
    assert rows[2]["tsd"] == "."


def test_read_rows_accepts_a_header_that_already_has_a_hash(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text("#" + TSV)
    rows = list(kmer2ltr.read_rows(p))
    assert len(rows) == 3
    assert rows[0]["seq_id"] == "chr1:100-200#LTR/Gypsy"


def test_read_rows_rejects_a_headerless_file(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text("chr1:100-200\t101\tpass\n")
    with pytest.raises(ValueError, match="header"):
        list(kmer2ltr.read_rows(p))


def test_status_counts(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text(TSV)
    assert kmer2ltr.status_counts(p) == {"pass": 2, "no_pair": 1}


def test_build_argv_uses_positional_input_and_dash_t():
    argv = kmer2ltr.build_argv(["Kmer2LTR"], "in.fa", "out.tsv", threads=8,
                               mutation_rate=3e-8, genome="g.fa", trim_flanks=True)
    assert argv[0] == "Kmer2LTR"
    assert argv[1] == "in.fa"
    assert "-i" not in argv
    assert "-p" not in argv
    assert argv[argv.index("-t") + 1] == "8"
    assert argv[argv.index("-o") + 1] == "out.tsv"
    assert argv[argv.index("-u") + 1] == "3e-08"
    assert argv[argv.index("--genome") + 1] == "g.fa"
    assert "--trim-flanks" in argv


def test_build_argv_omits_optional_flags():
    argv = kmer2ltr.build_argv(["Kmer2LTR"], "in.fa", "out.tsv")
    assert "--genome" not in argv
    assert "--trim-flanks" not in argv
    assert "-u" not in argv
    assert "--min-seq-id" not in argv


def test_build_argv_clustering():
    argv = kmer2ltr.build_argv(["Kmer2LTR"], "in.fa", "out.tsv",
                               ltr_cluster=True, internal_cluster=True, min_seq_id=0.75)
    assert "--ltr-cluster" in argv
    assert "--internal-cluster" in argv
    assert argv[argv.index("--min-seq-id") + 1] == "0.75"


def test_resolve_prefers_path(tmp_path, monkeypatch):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    exe = fake_bin / "Kmer2LTR"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", str(fake_bin))
    assert kmer2ltr.resolve(tmp_path / "tools") == [str(exe)]


def test_resolve_falls_back_to_a_clone_layout(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    tools = tmp_path / "tools"
    src = tools / "Kmer2LTR" / "src" / "kmer2ltr"
    src.mkdir(parents=True)
    (src / "__init__.py").write_text("")
    (src / "cli.py").write_text("def main():\n    return 0\n")
    argv = kmer2ltr.resolve(tools)
    assert argv[0] == f"PYTHONPATH={tools / 'Kmer2LTR' / 'src'}"
    assert argv[-2:] == ["-m", "kmer2ltr"]


def test_run_raises_on_failure(tmp_path):
    with pytest.raises(subprocess.CalledProcessError):
        kmer2ltr.run(["false"], "in.fa", str(tmp_path / "out.tsv"))

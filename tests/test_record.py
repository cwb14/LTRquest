"""Detection records: when a genome from an earlier run can be reused."""

from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

import ltrquest
from ltrquest import record

SETTINGS = {"max_rounds": "1", "terminate_count": "100", "trf": "true"}


@pytest.fixture
def run_dir(tmp_path: Path) -> Path:
    """A finished detection for prefix 'p', plus files that are not detection outputs."""
    for name in ("p_r1_ltr.tsv", "p_r1_ltr.fa", "p_depth0_ltr.tsv", "p_depth0_ltr.fa",
                 "p_strand_recovery.tsv",
                 # regenerated downstream, or another genome's
                 "p_depth0_clean_ltr.tsv", "p_depth0_clean_ltr.fa",
                 "p_all_depth_LTR_cleaned.gff3", "p_fpcheck.log", "p2_r1_ltr.tsv"):
        (tmp_path / name).write_text(name + "\n")
    (tmp_path / "p_r1.work").mkdir()
    (tmp_path / "p_r1.work" / "x.cls.tsv").write_text("x\n")
    (tmp_path / "p_plots").mkdir()
    (tmp_path / "genome.fa").write_text(">chr1\nACGTACGT\n")
    (tmp_path / "prot.fa").write_text(">p1\nMKV\n")
    return tmp_path


def check(d: Path, genome="genome.fa", proteins="prot.fa", settings=None):
    return record.check(str(d), "p", str(d / genome),
                        str(d / proteins) if proteins else None,
                        dict(SETTINGS if settings is None else settings))


def write(d: Path, genome="genome.fa", proteins="prot.fa", settings=None):
    return record.write(str(d), "p", str(d / genome),
                        str(d / proteins) if proteins else None,
                        dict(SETTINGS if settings is None else settings))


class TestDetectionOutputs:
    def test_lists_rounds_depth_tables_workdirs_and_sidecar(self, run_dir):
        assert record.detection_outputs(str(run_dir), "p") == [
            "p_depth0_ltr.fa", "p_depth0_ltr.tsv", "p_r1.work", "p_r1_ltr.fa",
            "p_r1_ltr.tsv", "p_strand_recovery.tsv"]

    def test_prefix_is_matched_exactly(self, run_dir):
        assert record.detection_outputs(str(run_dir), "p2") == ["p2_r1_ltr.tsv"]


class TestCheck:
    def test_no_record_means_detect(self, run_dir):
        assert check(run_dir) == ("detect", [], [])

    def test_unchanged_inputs_are_reused(self, run_dir):
        write(run_dir)
        assert check(run_dir) == ("reuse", [], [])

    def test_a_renamed_but_identical_genome_is_reused(self, run_dir):
        write(run_dir)
        (run_dir / "same.fa").write_bytes((run_dir / "genome.fa").read_bytes())
        assert check(run_dir, genome="same.fa")[0] == "reuse"

    def test_a_resized_genome_is_refused(self, run_dir):
        write(run_dir)
        (run_dir / "genome.fa").write_text(">chr1\nACGTACGTAA\n")
        verdict, problems, _ = check(run_dir)
        assert verdict == "refuse"
        assert any(p.startswith("genome:") for p in problems)

    def test_same_size_different_bases_is_refused(self, run_dir):
        write(run_dir)
        (run_dir / "genome.fa").write_text(">chr1\nTTTTACGT\n")
        assert check(run_dir)[0] == "refuse"

    def test_a_changed_setting_is_named(self, run_dir):
        write(run_dir)
        verdict, problems, _ = check(run_dir, settings={**SETTINGS, "max_rounds": "2"})
        assert verdict == "refuse"
        assert problems == ["max_rounds: 1 before, 2 now"]

    def test_an_added_setting_is_refused(self, run_dir):
        write(run_dir)
        verdict, problems, _ = check(run_dir, settings={**SETTINGS, "sdust": "true"})
        assert verdict == "refuse"
        assert problems == ["sdust: (unset) before, true now"]

    def test_a_removed_setting_is_refused(self, run_dir):
        write(run_dir)
        assert check(run_dir, settings={"max_rounds": "1"})[0] == "refuse"

    def test_proteins_added_is_refused(self, run_dir):
        write(run_dir, proteins=None)
        verdict, problems, _ = check(run_dir)
        assert verdict == "refuse"
        assert any(p.startswith("proteins:") for p in problems)

    def test_proteins_removed_is_refused(self, run_dir):
        write(run_dir)
        assert check(run_dir, proteins=None)[0] == "refuse"

    def test_proteins_changed_is_refused(self, run_dir):
        write(run_dir)
        (run_dir / "prot.fa").write_text(">p1\nMKVL\n")
        assert check(run_dir)[0] == "refuse"

    def test_a_missing_output_is_named(self, run_dir):
        write(run_dir)
        (run_dir / "p_depth0_ltr.tsv").unlink()
        verdict, problems, _ = check(run_dir)
        assert verdict == "refuse"
        assert problems == ["missing output: p_depth0_ltr.tsv"]

    def test_downstream_files_are_not_required(self, run_dir):
        write(run_dir)
        (run_dir / "p_depth0_clean_ltr.tsv").unlink()
        (run_dir / "p_all_depth_LTR_cleaned.gff3").unlink()
        assert check(run_dir)[0] == "reuse"

    def test_a_different_version_only_warns(self, run_dir):
        path = write(run_dir)
        saved = json.loads(Path(path).read_text())
        saved["ltrquest"] = "0.0.1"
        Path(path).write_text(json.dumps(saved))
        verdict, problems, warnings = check(run_dir)
        assert (verdict, problems) == ("reuse", [])
        assert warnings == [f"detected with LTRquest 0.0.1; this is {ltrquest.__version__}"]

    def test_numerically_equal_settings_match(self, run_dir):
        write(run_dir, settings={**SETTINGS, "fp_mask_threshold": "0.10"})
        assert check(run_dir, settings={**SETTINGS, "fp_mask_threshold": "0.1"})[0] == "reuse"

    def test_an_fp_masked_genome_is_reused_with_a_warning(self, run_dir):
        record.write(str(run_dir), "p", str(run_dir / "genome.fa"), None, SETTINGS,
                     fp_masked=True)
        verdict, _problems, warnings = check(run_dir, proteins=None)
        assert verdict == "reuse"
        assert any("FP-masked" in w for w in warnings)

    def test_an_unreadable_record_is_refused(self, run_dir):
        (run_dir / "p.detect.json").write_text("{not json")
        assert check(run_dir)[0] == "refuse"


class TestWrite:
    def test_records_what_went_in_and_came_out(self, run_dir):
        saved = json.loads(Path(write(run_dir)).read_text())
        assert saved["prefix"] == "p"
        assert saved["genome"]["name"] == "genome.fa"
        assert saved["genome"]["size"] == (run_dir / "genome.fa").stat().st_size
        assert saved["proteins"]["name"] == "prot.fa"
        assert saved["settings"] == SETTINGS
        assert "p_r1.work" in saved["outputs"]

    def test_is_readable_like_any_other_output(self, run_dir):
        umask = os.umask(0)
        os.umask(umask)
        mode = stat.S_IMODE(os.stat(write(run_dir)).st_mode)
        assert mode == 0o666 & ~umask

    def test_refuses_a_genome_with_no_outputs(self, tmp_path):
        (tmp_path / "genome.fa").write_text(">chr1\nACGT\n")
        with pytest.raises(ValueError):
            record.write(str(tmp_path), "p", str(tmp_path / "genome.fa"), None, {})


class TestClear:
    def test_removes_what_detection_wrote_and_nothing_else(self, run_dir):
        write(run_dir)
        (run_dir / "p_depth3_ltr.tsv").write_text("stale\n")
        (run_dir / "p_depth0_clean_ltr.filtered.tsv").write_text("mine\n")
        (run_dir / "p_r1_my_blast").mkdir()
        removed = record.clear(str(run_dir), "p")
        assert removed == [
            "p_depth0_clean_ltr.fa", "p_depth0_clean_ltr.tsv", "p_depth0_ltr.fa",
            "p_depth0_ltr.tsv", "p_depth3_ltr.tsv", "p_r1.work", "p_r1_ltr.fa",
            "p_r1_ltr.tsv", "p_strand_recovery.tsv", "p.detect.json"]
        assert sorted(n.name for n in run_dir.iterdir()) == [
            "genome.fa", "p2_r1_ltr.tsv", "p_all_depth_LTR_cleaned.gff3",
            "p_depth0_clean_ltr.filtered.tsv", "p_fpcheck.log", "p_plots",
            "p_r1_my_blast", "prot.fa"]


class TestCli:
    def cli(self, *args):
        return subprocess.run([sys.executable, "-m", "ltrquest.record", *args],
                              capture_output=True, text=True)

    def common(self, d: Path, max_rounds="1"):
        return ["--indir", str(d), "--prefix", "p", "--genome", str(d / "genome.fa"),
                "--setting", f"max_rounds={max_rounds}"]

    def test_check_write_check(self, run_dir):
        first = self.cli("check", *self.common(run_dir))
        assert (first.returncode, first.stdout) == (0, "detect\n")
        assert self.cli("write", *self.common(run_dir)).returncode == 0
        second = self.cli("check", *self.common(run_dir))
        assert (second.returncode, second.stdout) == (0, "reuse\n")

    def test_refusal_exits_nonzero_and_says_why(self, run_dir):
        self.cli("write", *self.common(run_dir))
        result = self.cli("check", *self.common(run_dir, max_rounds="2"))
        assert result.returncode == 1
        assert result.stdout == ""
        assert "max_rounds: 1 before, 2 now" in result.stderr

    def test_outputs_prints_the_recorded_names(self, run_dir):
        self.cli("write", *self.common(run_dir))
        result = self.cli("outputs", "--indir", str(run_dir), "--prefix", "p")
        assert result.stdout.split() == record.detection_outputs(str(run_dir), "p")

    def test_write_records_fp_masking(self, run_dir):
        self.cli("write", *self.common(run_dir), "--fp-masked")
        saved = json.loads((run_dir / "p.detect.json").read_text())
        assert saved["fp_masked"] is True

    def test_clear_prints_what_it_removed(self, run_dir):
        result = self.cli("clear", "--indir", str(run_dir), "--prefix", "p")
        assert result.returncode == 0
        assert "p_r1.work" in result.stdout.split()
        assert not (run_dir / "p_r1.work").exists()

    def test_a_setting_without_equals_is_a_usage_error(self, run_dir):
        result = self.cli("check", "--indir", str(run_dir), "--prefix", "p",
                          "--genome", str(run_dir / "genome.fa"), "--setting", "oops")
        assert result.returncode == 2

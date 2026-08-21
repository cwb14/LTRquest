"""The installed package: entry points, packaged scripts, and their wiring.

These are the checks that catch a broken `pip install .` before a user does.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import ltrquest
from ltrquest import cli

STAGES = ["detect", "mask", "reconcile", "annotate", "gff3",
          "plot_struct", "plot_summary", "tegv"]
CONSOLE_SCRIPTS = ["ltrquest", "ltrquest-plots", "ltrquest-detect", "ltrquest-mask",
                   "ltrquest-reconcile", "ltrquest-annotate", "ltrquest-gff3"]
REPO_ROOT = Path(__file__).resolve().parents[1]


def run(*args, **kwargs):
    return subprocess.run(args, capture_output=True, text=True, **kwargs)


class TestPackagedScripts:
    def test_both_shell_scripts_ship_with_the_package(self):
        scripts = Path(ltrquest.__file__).resolve().parent / "scripts"
        assert (scripts / "ltrquest.sh").is_file()
        assert (scripts / "plots.sh").is_file()

    @pytest.mark.parametrize("name", ["ltrquest.sh", "plots.sh"])
    def test_scripts_are_executable(self, name):
        scripts = Path(ltrquest.__file__).resolve().parent / "scripts"
        assert (scripts / name).stat().st_mode & 0o111

    @pytest.mark.parametrize("name", ["ltrquest.sh", "plots.sh"])
    def test_scripts_parse(self, name):
        scripts = Path(ltrquest.__file__).resolve().parent / "scripts"
        result = run("bash", "-n", str(scripts / name))
        assert result.returncode == 0, result.stderr

    def test_cli_resolves_the_driver(self, driver):
        assert cli._script("ltrquest.sh") == str(driver)

    def test_cli_reports_a_missing_script_rather_than_crashing(self):
        with pytest.raises(SystemExit) as excinfo:
            cli._script("not_a_script.sh")
        assert "missing" in str(excinfo.value)


class TestStagesAreImportable:
    @pytest.mark.parametrize("stage", STAGES)
    def test_stage_imports(self, stage):
        __import__(f"ltrquest.{stage}")

    @pytest.mark.parametrize("stage", STAGES)
    def test_stage_runs_as_a_module(self, stage):
        result = run(sys.executable, "-m", f"ltrquest.{stage}", "--help")
        assert result.returncode == 0, result.stderr
        assert "usage" in (result.stdout + result.stderr).lower()

    def test_cross_module_imports_resolve_inside_the_package(self):
        # reconcile borrows containment helpers from detect, and gff3 borrows
        # the table readers from annotate. Both were sys.path hacks before.
        from ltrquest.gff3 import collect_elements  # noqa: F401
        from ltrquest.reconcile import iter_fasta  # noqa: F401


class TestConsoleScripts:
    @pytest.mark.parametrize("name", CONSOLE_SCRIPTS)
    def test_console_script_is_on_path(self, name):
        assert shutil.which(name), f"{name} was not installed"

    @pytest.mark.parametrize("name", CONSOLE_SCRIPTS)
    def test_console_script_answers_help(self, name):
        result = run(name, "--help")
        assert result.returncode == 0, result.stderr

    def test_driver_help_documents_the_required_flag(self):
        result = run("ltrquest", "--help")
        assert "--genome" in result.stdout + result.stderr

    def test_driver_help_does_not_mention_the_removed_script_path_flag(self):
        text = run("ltrquest", "--help").stdout + run("ltrquest", "--help").stderr
        assert "--script_path" not in text

    def test_driver_echoes_the_command_it_was_invoked_as(self):
        # The first line is the reproducibility record; it must read as the
        # command a user could paste back, not as a site-packages path.
        result = run("ltrquest", "--genome", "definitely_missing.fa")
        assert result.stdout.splitlines()[0] == "Command: ltrquest --genome definitely_missing.fa"


class TestDriverGuards:
    def test_no_arguments_prints_usage_and_fails(self):
        result = run("ltrquest")
        assert result.returncode != 0
        assert "Usage:" in result.stderr

    def test_a_missing_genome_is_rejected_by_name(self, tmp_path):
        result = run("ltrquest", "--genome", str(tmp_path / "nope.fa"))
        assert result.returncode != 0
        assert "Genome not found" in result.stderr

    def test_an_unknown_flag_is_rejected(self, toy_genome):
        result = run("ltrquest", "--genome", str(toy_genome), "--not-a-flag")
        assert result.returncode != 0
        assert "Unknown argument" in result.stderr

    def test_an_out_of_range_threshold_is_rejected(self, toy_genome):
        result = run("ltrquest", "--genome", str(toy_genome), "--fp-mask-threshold", "7")
        assert result.returncode != 0
        assert "fp-mask-threshold" in result.stderr

    def test_the_driver_refuses_an_interpreter_without_the_package(self, toy_genome, driver):
        # The driver dispatches every Python stage through $LTRQUEST_PYTHON; if
        # that interpreter cannot import ltrquest it must say so up front rather
        # than failing several minutes into round 1.
        result = run("bash", str(driver), "--genome", str(toy_genome),
                     env={"PATH": "/usr/bin:/bin", "LTRQUEST_PYTHON": "/nonexistent/python"})
        assert result.returncode != 0
        assert "not importable" in result.stderr


class TestVersion:
    def test_package_exposes_a_version(self):
        assert re.fullmatch(r"\d+\.\d+\.\d+", ltrquest.__version__)

    @pytest.mark.skipif(not (REPO_ROOT / "pyproject.toml").is_file(),
                        reason="not running from a source checkout")
    def test_version_matches_pyproject(self):
        text = (REPO_ROOT / "pyproject.toml").read_text()
        declared = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
        assert declared, "no version in pyproject.toml"
        assert declared.group(1) == ltrquest.__version__

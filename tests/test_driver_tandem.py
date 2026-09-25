"""The driver's --no-tandem-filter wiring: flag, FP-stage pass-through, record."""

from __future__ import annotations

import subprocess


def run(*args, **kw):
    return subprocess.run(args, capture_output=True, text=True, **kw)


def block(text: str, start: str, end: str) -> str:
    return text[text.index(start):text.index(end, text.index(start))]


def test_help_documents_the_flag(driver):
    result = run("bash", str(driver), "--help")
    assert "--no-tandem-filter" in result.stdout + result.stderr


def test_the_flag_parses(driver):
    assert run("bash", str(driver), "--no-tandem-filter", "--help").returncode == 0


def test_the_pooled_fp_call_judges_every_genome_unless_the_flag_is_given(driver):
    stage = block(driver.read_text(), "run_fp_stage() {", "write_clean_fastas() {")
    assert 'if [[ "$RUN_TANDEM" == true ]]; then\n    tandem_opts=( --tandem-genome )' in stage
    assert 'for i in "${!OUT_PREFIXES[@]}"; do' in stage
    # the first flag_fp call (the one writing _clean_ tables), not the masking ones
    first_call = block(stage, '--domains-tsv "${dtsvs[@]}"', "--threads")
    assert '"${tandem_opts[@]}"' in first_call


def test_it_is_not_a_detection_setting(driver):
    # Detection never sees the filter, so a reused genome needs no re-detection
    # to gain or lose it: the pooled stage re-judges every genome on each run.
    settings = block(driver.read_text(), "detection_settings() {", "plan_reuse() {")
    assert "TANDEM" not in settings

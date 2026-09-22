"""The driver's --reboundary wiring: flag, early check, stage order, file hygiene."""

from __future__ import annotations

import subprocess


def run(*args, **kw):
    return subprocess.run(args, capture_output=True, text=True, **kw)


def test_help_documents_the_flag(driver):
    result = run("bash", str(driver), "--help")
    assert "--reboundary" in result.stdout + result.stderr


def test_the_flag_parses(driver):
    assert run("bash", str(driver), "--reboundary", "--help").returncode == 0


def test_a_missing_mafft_fails_before_any_work(driver, toy_genome):
    result = run("bash", str(driver), "--genome", str(toy_genome), "--reboundary",
                 env={"PATH": "/usr/bin:/bin", "LTRQUEST_PYTHON": "/nonexistent/python"})
    assert result.returncode != 0 and "mafft" in result.stderr


def test_it_is_not_a_detection_setting(driver):
    text = driver.read_text()
    block = text[text.index("detection_settings() {"):text.index("plan_reuse() {")]
    assert "REBOUNDARY" not in block


def test_carry_forward_and_promotion_drop_its_files(driver):
    text = driver.read_text()
    assert "*_reboundary.*|*_pre_reboundary) continue;;" in text
    assert 'rm -f "${dst}/${p}_reboundary.tsv"' in text
    assert 'rm -rf "${dst}/${p}_pre_reboundary"' in text


def test_stage_order(driver):
    text = driver.read_text()
    stage = text[text.index("run_annotation_stage() {"):text.index("run_reboundary_stage() {")]
    off_gff3 = stage.index('"${GFF3[@]}" --prefix "$p" --indir . \\\n      --genome "${p}.input_genome.fa" "${fam_opts[@]}" "${rec_opts[@]}"\n')
    on_skip = stage.index('if [[ "$REBOUNDARY" == true ]]; then')
    assert on_skip < off_gff3                            # with the flag, loop 1 stops before gff3
    assert stage.index("run_reboundary_stage\n") > off_gff3   # then pooled, then gff3 with the map
    assert '--reboundary-map "${p}_reboundary.tsv"' in stage

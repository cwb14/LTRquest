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


def test_a_missing_blastn_fails_before_any_work(driver, toy_genome, tmp_path):
    result = run("bash", str(driver), "--genome", str(toy_genome), "--reboundary",
                 env={"PATH": "/usr/bin:/bin", "LTRQUEST_PYTHON": "/nonexistent/python"},
                 cwd=tmp_path)
    assert result.returncode != 0 and "blastn" in result.stderr


def test_an_incompatible_kmer2ltr_stops_the_run_before_detection(driver, toy_genome, tmp_path):
    """A stand-in interpreter imports ltrquest fine but reports a Kmer2LTR that
    force-pairs credited flanks: the run must stop with exit 3, not detect for hours
    and then keep the calls as detected."""
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    for tool in ("blastn", "makeblastdb"):
        stub = bin_ / tool
        stub.write_text("#!/bin/sh\nexit 0\n")
        stub.chmod(0o755)
    py = bin_ / "python-stub"
    py.write_text('#!/bin/sh\ncase "$2" in *IncompatibleKmer2LTR*) '
                  'echo "ERROR: --reboundary: force-pairs credited flanks" >&2; exit 3;; esac\n'
                  'exit 0\n')
    py.chmod(0o755)
    result = run("bash", str(driver), "--genome", str(toy_genome), "--reboundary",
                 env={"PATH": f"{bin_}:/usr/bin:/bin", "LTRQUEST_PYTHON": str(py)}, cwd=tmp_path)
    assert result.returncode == 3 and "force-pairs credited flanks" in result.stderr


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

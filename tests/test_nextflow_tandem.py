"""The Nextflow pipeline's tandem-array wiring: every genome reaches flag_fp.

Without a genome, ltrquest.flag_fp skips the tandem-array purge silently, so a
pipeline that forgot to stage them would still run and still exit 0.
"""

import json
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE = REPO_ROOT / "modules" / "local" / "ltrquest" / "flagfp" / "main.nf"
WORKFLOW = REPO_ROOT / "workflows" / "ltrquest.nf"

pytestmark = pytest.mark.skipif(not MODULE.is_file(), reason="not running from a source checkout")


def section(text: str, start: str, end: str) -> str:
    i = text.index(start)
    return text[i:text.index(end, i)]


def test_the_module_stages_the_genomes_apart_from_the_masking_genome():
    inputs = section(MODULE.read_text(), "input:", "output:")
    # numbered directories: two samples' genome.fa must not collide, nor with `genome`
    assert "path(genome), path(tandem_genomes, stageAs: 'genomes/g??/*')" in inputs


def test_the_module_hands_them_to_flag_fp_and_publishes_the_scores():
    text = MODULE.read_text()
    script = section(text, "script:", "stub:")
    assert """def tandem = tandem_genomes ? "--tandem-genome ${tandem_genomes}" : ''""" in script
    call = section(script, "python -m ltrquest.flag_fp", "2>&1 | tee")
    assert "${tandem}" in call
    assert 'path("${prefix}_fpcheck.tandem.tsv"), emit: tandem, optional: true' in text


def test_the_workflow_feeds_every_sample_unless_skipped():
    text = WORKFLOW.read_text()
    feed = section(text, "ch_tandem_genomes = ", "LTRQUEST_FLAGFP(")
    assert "params.skip_tandem_filter" in feed and "Channel.value([])" in feed
    assert "ch_samples.map { _meta, genome, _proteins -> genome }" in feed
    assert ".toSortedList" in feed                         # stable -resume hash
    call = section(text, "LTRQUEST_FLAGFP(", "ch_versions")
    assert ".combine(ch_tandem_genomes.map { files -> [files] })" in call


def test_the_parameter_is_declared_once_documented_and_published():
    config = (REPO_ROOT / "nextflow.config").read_text()
    assert re.search(r"^\s*skip_tandem_filter\s*=\s*false\s*$", config, re.MULTILINE)
    schema = json.loads((REPO_ROOT / "nextflow_schema.json").read_text())
    props = {k: v for d in schema["$defs"].values() for k, v in d.get("properties", {}).items()}
    assert props["skip_tandem_filter"]["type"] == "boolean"
    assert "--skip_tandem_filter" in (REPO_ROOT / "main.nf").read_text()
    modules_config = (REPO_ROOT / "conf" / "modules.config").read_text()
    assert "pattern: '*_fpcheck.{log,tandem.tsv}'" in modules_config

"""The Nextflow reconcile step with a genome that completed a single round.

`--max_rounds 1`, or a first round that finds fewer than `--terminate_count`
elements, leaves one table per genome. Nextflow then hands the process a single
Path rather than a list, and collect() on a Path walks its name components: the
reconciler was called as `--tsv round_1 x_r1_ltr.tsv` and died on the directory.
The three-round stub test never reaches this case.
"""

from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[1] / "modules" / "local" / "ltrquest" / "reconcile" / "main.nf"


@pytest.mark.skipif(not MODULE.is_file(), reason="not running from a source checkout")
def test_a_single_table_is_wrapped_before_it_is_iterated():
    script = MODULE.read_text()
    for name in ("tsvs", "fastas"):
        assert f"({name}" in script and f"instanceof List ? {name}" in script
    assert "tsvs.collect" not in script and "fastas.collect" not in script

"""The benchmark harness runs end to end on the synthetic run."""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import pytest

from reboundary_fixtures import build

BENCH = Path(__file__).resolve().parents[1] / "benchmarks" / "reboundary" / "bench.py"


def load():
    spec = importlib.util.spec_from_file_location("bench", BENCH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bench"] = mod          # multiprocessing pickles bench._b1_family by name
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.skipif(not BENCH.is_file(), reason="benchmarks are not in this checkout")
def test_b1_b2_and_collect_run_on_the_synthetic_run(tmp_path, k2l_api, mafft, monkeypatch,
                                                    capsys):
    bench = load()
    fx = build(tmp_path / "syn")
    common = ["--run", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
              "--tools-dir", os.environ.get("LTRQUEST_TOOLS_DIR", "."), "--threads", "2",
              "--mafft", mafft]
    monkeypatch.setattr(sys, "argv", ["bench.py", "b1", *common, "--out",
                                      str(tmp_path / "t" / "b1_x"), "--n", "3"])
    bench.main()
    s1 = json.loads((tmp_path / "t" / "b1_x" / "summary.json").read_text())
    assert s1["n_truth"] == 3 and s1["exact_or_1"] is not None and s1["false_change"] == 0.0
    monkeypatch.setattr(sys, "argv", ["bench.py", "b2", *common, "--out",
                                      str(tmp_path / "t" / "b2_x")])
    bench.main()
    s2 = json.loads((tmp_path / "t" / "b2_x" / "summary.json").read_text())
    assert s2["extended"] == 4 and s2["tsd_gain"] == 1.0
    monkeypatch.setattr(sys, "argv", ["bench.py", "collect", "--out", str(tmp_path / "t")])
    bench.main()
    out = capsys.readouterr().out
    assert "| x |" in out and "B1 planted obstacles" in out and "B2 real run" in out

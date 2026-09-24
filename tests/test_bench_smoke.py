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
    sys.modules["bench"] = mod          # multiprocessing pickles bench._b1_place by name
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.skipif(not BENCH.is_file(), reason="benchmarks are not in this checkout")
def test_b1_b2_and_collect_run_on_the_synthetic_run(tmp_path, k2l_api, blastn, monkeypatch,
                                                    capsys):
    bench = load()
    fx = build(tmp_path / "syn")
    common = ["--run", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
              "--tools-dir", os.environ.get("LTRQUEST_TOOLS_DIR", "."), "--threads", "2"]
    monkeypatch.setattr(sys, "argv", ["bench.py", "b1", *common, "--out",
                                      str(tmp_path / "t" / "b1_x"), "--n", "3"])
    bench.main()
    s1 = json.loads((tmp_path / "t" / "b1_x" / "summary.json").read_text())
    assert s1["n_truth"] == 3 and s1["exact_or_1"] is not None and s1["false_change"] == 0.0
    monkeypatch.setattr(sys, "argv", ["bench.py", "b2", *common, "--out",
                                      str(tmp_path / "t" / "b2_x")])
    bench.main()
    s2 = json.loads((tmp_path / "t" / "b2_x" / "summary.json").read_text())
    assert s2["moved"] == 9 and s2["merged"] == 1 and s2["tsd_gain"] == 1.0
    assert s2["spike_at_0"] > 0.5 and s2["spike_background"] < 0.05
    monkeypatch.setattr(sys, "argv", ["bench.py", "collect", "--out", str(tmp_path / "t")])
    bench.main()
    out = capsys.readouterr().out
    assert "| x |" in out and "B1 planted obstacles" in out and "B2 real run" in out


@pytest.mark.skipif(not BENCH.is_file(), reason="benchmarks are not in this checkout")
def test_b1_refuses_an_incompatible_kmer2ltr_before_starting_workers(tmp_path, blastn,
                                                                      monkeypatch):
    """Pool workers whose initializer raises are restarted forever, so a Kmer2LTR the
    workers would refuse must stop b1 before its pool starts."""
    bench = load()
    fx = build(tmp_path / "syn")

    def too_old(*args, **kwargs):
        raise bench.rb.k2l.IncompatibleKmer2LTR("this Kmer2LTR force-pairs credited flanks")

    def no_pool(*args, **kwargs):
        raise AssertionError("a worker pool started before the Kmer2LTR check")

    monkeypatch.setattr(bench.rb.k2l, "api", too_old)
    monkeypatch.setattr(bench, "Pool", no_pool)
    monkeypatch.setattr(sys, "argv", ["bench.py", "b1", "--run", str(fx.indir), "--prefix",
                                      fx.prefix, "--genome", str(fx.genome), "--tools-dir",
                                      ".", "--threads", "2", "--out", str(tmp_path / "b1"),
                                      "--n", "3"])
    with pytest.raises(SystemExit, match="force-pairs"):
        bench.main()

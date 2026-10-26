"""FastLTR, the third LTR caller: masked-base removal, coordinate lifting, SCN
rows, build resolution, the per-chunk mode switch, and the driver wiring.

Everything but the `slow` class runs without a compiler or FASTAN checkout.
"""

from __future__ import annotations

import random
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from ltrquest import detect
from ltrquest import fastltr as fl


def fake_binary(path: Path, usage: str) -> Path:
    """An executable that prints `usage` on stderr and exits 1, like FastLTR with no args."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!/bin/sh\necho '{usage}' >&2\nexit 1\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


class TestCompress:
    def test_masked_runs_are_deleted_and_every_base_lifts_back(self):
        seq = "ACGTNNNNacgtRRRRRRggcc"
        (piece,) = fl.compress(seq, min_piece=0)
        assert piece.seq == "ACGTACGTGGCC"
        lifted = [piece.lift(i) for i in range(len(piece.seq))]
        assert [seq[i].upper() for i in lifted] == list(piece.seq)
        assert lifted == [0, 1, 2, 3, 8, 9, 10, 11, 18, 19, 20, 21]

    def test_a_far_mask_run_splits_instead_of_joining(self):
        left, right = "ACGT" * 5, "TTGCA" * 4
        pieces = fl.compress(left + "VVVVVVVVVV" + right, split_char="V", min_piece=0)
        assert [p.seq for p in pieces] == [left, right]
        assert pieces[1].lift(0) == len(left) + 10

    def test_without_a_split_char_the_far_mask_is_deleted_like_any_mask(self):
        pieces = fl.compress("ACGT" * 5 + "VVVV" + "TTGCA" * 4, min_piece=0)
        assert len(pieces) == 1

    def test_pieces_too_short_for_a_pair_are_dropped(self):
        pieces = fl.compress("ACGTACGTAC" + "N" * 50 + "A" * 300, min_piece=200)
        assert [len(p.seq) for p in pieces] == [310]

    def test_lift_refuses_a_position_outside_the_piece(self):
        (piece,) = fl.compress("ACGTACGTACGT", min_piece=0)
        with pytest.raises(ValueError):
            piece.lift(len(piece.seq))


class TestScnRows:
    def test_row_is_one_based_inclusive_in_chunk_coordinates(self):
        (piece,) = fl.compress("A" * 1000, min_piece=0)
        row = fl.scn_line(piece, 100, 300, 700, 900, diffs=10).split()
        assert row[:9] == ["101", "900", "800", "101", "300", "200", "701", "900", "200"]
        assert float(row[9]) == pytest.approx(95.0)
        assert row[10] == "0"

    def test_a_copy_spanning_a_deleted_mask_is_reported_over_the_mask(self):
        # 100 bp masked in the middle of the 5' copy: the row spans it.
        (piece,) = fl.compress("C" * 150 + "N" * 100 + "C" * 850, min_piece=0)
        row = fl.scn_line(piece, 100, 300, 700, 900, diffs=0).split()
        assert row[3:6] == ["101", "400", "300"]
        assert row[6:8] == ["801", "1000"]

    def test_short_copies_and_junk_lines_are_skipped(self):
        rng = random.Random(9)
        ltr = "".join(rng.choice("ACGT") for _ in range(140))
        filler = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        # unique sequence with one pair at [10, 150) and [1010, 1150)
        pieces = fl.compress(filler(10) + ltr + filler(860) + ltr + filler(850), min_piece=0)
        text = "\n".join([
            "p0\t10\t150\t1010\t1150\t3\t1000",     # kept
            "p0\t10\t60\t1010\t1150\t3\t1000",      # 5' copy 50 bp < 100
            "p7\t10\t150\t1010\t1150\t3\t1000",     # no such piece
            "garbage",
            "p0\tx\t150\t1010\t1150\t3\t1000",
        ])
        rows = fl.hits_to_scn(text, pieces, min_ltr=100)
        assert len(rows) == 1 and rows[0].startswith("11 1150 ")


class TestArrayGuard:
    """An LTR-RT has two copies; a tandem array has more at a near-constant spacing."""

    @staticmethod
    def rnd(rng, n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def test_a_pair_of_neighbouring_array_units_is_periodic(self):
        rng = random.Random(1)
        unit = self.rnd(rng, 300)
        array = "".join(unit + self.rnd(rng, 150) for _ in range(6))
        assert fl.periodic(array, 450, 750, 900, 1200, identity=1.0)

    def test_the_end_units_of_an_array_still_count(self):
        # the 5' copy of the first pair has nothing upstream; the 3' copy has a
        # third unit downstream
        rng = random.Random(2)
        unit = self.rnd(rng, 300)
        array = "".join(unit + self.rnd(rng, 150) for _ in range(3))
        assert fl.periodic(array, 0, 300, 450, 750, identity=1.0)

    def test_an_ltr_pair_in_unique_sequence_is_not(self):
        rng = random.Random(3)
        ltr = self.rnd(rng, 300)
        seq = self.rnd(rng, 2000) + ltr + self.rnd(rng, 3000) + ltr + self.rnd(rng, 2000)
        assert not fl.periodic(seq, 2000, 2300, 5300, 5600, identity=1.0)

    def test_a_diverged_third_copy_does_not_make_an_array(self):
        # identity far below the pair's own: chance similarity, not a unit
        rng = random.Random(4)
        ltr = self.rnd(rng, 300)
        far = TestRealFastLTR.mutate(ltr, 0.45, rng)
        seq = self.rnd(rng, 600) + far + self.rnd(rng, 600) + ltr + self.rnd(rng, 600) + ltr + self.rnd(rng, 600)
        bs = 600 + 300 + 600
        assert not fl.periodic(seq, bs, bs + 300, bs + 900, bs + 1200, identity=1.0)

    def test_hits_to_scn_drops_array_units_and_counts_them(self):
        rng = random.Random(5)
        unit = self.rnd(rng, 300)
        (piece,) = fl.compress("".join(unit + self.rnd(rng, 150) for _ in range(6)), min_piece=0)
        stats = {}
        rows = fl.hits_to_scn("p0\t450\t750\t900\t1200\t0\t450", [piece], 100, stats)
        assert rows == [] and stats == {"periodic": 1}


class TestBuildResolution:
    def test_only_the_adapted_build_is_usable(self, tmp_path):
        ours = fake_binary(tmp_path / "a" / "FastLTR", f"Usage ... ({fl.MARKER}, FASTAN x)")
        stock = fake_binary(tmp_path / "b" / "FastLTR", "Usage: FastLTR [-vM] [-T(8)]")
        assert fl.tool_usable(ours)
        assert not fl.tool_usable(stock)
        assert not fl.tool_usable(tmp_path / "missing")

    def test_a_previous_build_in_tools_dir_wins(self, tmp_path, monkeypatch):
        built = fake_binary(tmp_path / "FASTAN" / "FastLTR", fl.MARKER)
        monkeypatch.setattr("shutil.which", lambda n: None)
        assert fl.ensure(tmp_path) == str(built.resolve())

    def test_a_stock_fastltr_on_path_is_not_taken(self, tmp_path, monkeypatch):
        stock = fake_binary(tmp_path / "bin" / "FastLTR", "Usage: FastLTR [-vM]")
        monkeypatch.setattr("shutil.which", lambda n: str(stock) if n == "FastLTR" else None)
        with pytest.raises(RuntimeError) as excinfo:
            fl.ensure(tmp_path / "tools")
        assert "git is not installed" in str(excinfo.value)
        assert "--no-fastltr" in str(excinfo.value)

    def test_the_source_ships_with_the_package(self):
        text = fl.SOURCE.read_text()
        assert fl.MARKER in text
        assert "FASTAN commit 843d0ed" in text
        assert fl.FASTAN_REF.startswith("843d0ed")


class TestWhenFastLTRRuns:
    @pytest.mark.parametrize("mode,cap,active", [
        ("all", 0, True), ("all", 30, True),
        ("rescue", 30, True), ("rescue", 0, False),
        ("off", 30, False), ("off", 0, False)])
    def test_rescue_without_a_cap_never_runs_or_builds_it(self, mode, cap, active):
        assert detect.fastltr_active(mode, cap) is active


class TestTimeoutText:
    def test_partial_stdout_arrives_as_bytes_and_is_decoded(self):
        try:
            subprocess.run(["sh", "-c", "echo found-one; sleep 5"], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True, timeout=1)
        except subprocess.TimeoutExpired as te:
            assert detect._timeout_text(te) == "found-one\n"
        else:  # pragma: no cover
            pytest.fail("the call was expected to time out")


class TestModeSwitch:
    """process_one_chunk with the three callers stubbed out."""

    @pytest.fixture
    def chunk(self, tmp_path):
        fa = tmp_path / "c.fa"
        fa.write_text(">chr1__chunk1__1-1000\n" + "ACGT" * 250 + "\n")
        return detect.ChunkInfo("chr1.chunk1.1-1000", "chr1", 0, 1000, str(fa))

    def run(self, monkeypatch, tmp_path, chunk, mode, capped):
        calls = []
        monkeypatch.setattr(detect, "build_suffixerator_index", lambda *a, **k: None)
        monkeypatch.setattr(detect, "run_ltrharvest_scn_and_gff3",
                            lambda *a, **k: (Path(a[2]).write_text(""), capped)[1])
        monkeypatch.setattr(detect, "run_ltrfinder_raw",
                            lambda *a, **k: (Path(a[2]).write_text(""), False)[1])
        monkeypatch.setattr(detect.fl, "run_chunk", lambda *a, **k: calls.append((a, k)) or 0)
        out = detect.process_one_chunk(
            "gt", "ltr_finder", chunk, str(tmp_path / "work"), [], ["-w", "2"], "both",
            ltr_timeout_s=30, fastltr_bin="FastLTR", fastltr_args=[], fastltr_mode=mode)
        return out, calls

    def test_all_runs_fastltr_on_every_chunk(self, monkeypatch, tmp_path, chunk):
        out, calls = self.run(monkeypatch, tmp_path, chunk, "all", capped=False)
        assert len(calls) == 1 and out[4] is not None and out[5] == set() and out[6] == {}

    def test_rescue_skips_a_chunk_the_other_two_finished(self, monkeypatch, tmp_path, chunk):
        out, calls = self.run(monkeypatch, tmp_path, chunk, "rescue", capped=False)
        assert calls == [] and out[4] is None

    def test_rescue_searches_a_chunk_ltrharvest_was_cut_short_on(self, monkeypatch, tmp_path, chunk):
        out, calls = self.run(monkeypatch, tmp_path, chunk, "rescue", capped=True)
        assert len(calls) == 1 and out[5] == {"LTRharvest"}

    def test_fastltr_gets_the_same_cap_as_the_other_two(self, monkeypatch, tmp_path, chunk):
        _, calls = self.run(monkeypatch, tmp_path, chunk, "all", capped=False)
        assert calls[0][1]["timeout_s"] == 30


def test_run_chunk_passes_the_cap_through(monkeypatch, tmp_path):
    fa = tmp_path / "c.fa"
    fa.write_text(">x\n" + "ACGT" * 300 + "\n")
    seen = {}

    def fake_run(cmd, **kw):
        seen.update(kw)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(fl.subprocess, "run", fake_run)
    fl.run_chunk("FastLTR", str(fa), str(tmp_path), str(tmp_path / "o.scn"), [], timeout_s=30)
    assert seen["timeout"] == 30


class TestDriverWiring:
    def test_help_documents_the_flags(self, driver):
        out = subprocess.run(["bash", str(driver), "--help"], capture_output=True, text=True)
        text = out.stdout + out.stderr
        for flag in ("--fastltr MODE", "--no-fastltr", "--fastltr-args"):
            assert flag in text

    def test_a_bad_mode_is_refused_before_any_work(self, driver, tmp_path):
        genome = tmp_path / "g.fa"
        genome.write_text(">c\nACGT\n")
        out = subprocess.run(["bash", str(driver), "--genome", str(genome), "--fastltr", "some"],
                             capture_output=True, text=True, cwd=tmp_path)
        assert out.returncode != 0
        assert "--fastltr must be all, rescue or off" in out.stdout + out.stderr

    def test_the_mode_is_a_detection_setting_only_when_it_can_act(self, driver):
        # off, or rescue without a cap, detects exactly what a pre-FastLTR run did,
        # so a record made before FastLTR existed must still match
        text = driver.read_text()
        settings = text[text.index("detection_settings() {"):text.index("plan_reuse() {")]
        assert ('if [[ "$FASTLTR_MODE" == all || ( "$FASTLTR_MODE" == rescue && '
                '"$directives" == *ltr-timeout* ) ]]; then') in settings
        assert '--setting "fastltr=${FASTLTR_MODE}"' in settings

    def test_every_round_gets_the_mode(self, driver):
        text = driver.read_text()
        call = text[text.index('"${DETECT[@]}" \\'):text.index('|| ltrharvest_exit=$?')]
        assert '"${FASTLTR_OPTS[@]}"' in call


@pytest.fixture(scope="module")
def binary(tmp_path_factory):
    """The adapted FastLTR, built from a fresh FASTAN clone."""
    if shutil.which("git") is None or not any(shutil.which(c) for c in ("gcc", "cc")):
        pytest.skip("needs git and a C compiler")
    return fl.ensure(tmp_path_factory.mktemp("tools"))


@pytest.mark.slow
class TestRealFastLTR:
    """Build the adapted FastLTR (clones FASTAN) and find a planted LTR pair."""

    @staticmethod
    def mutate(seq, rate, rng):
        return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b for b in seq)

    def test_a_planted_pair_is_found_through_a_masked_inner(self, binary, tmp_path):
        rng = random.Random(7)
        rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        ltr = rand(400)
        left, inner, right = rand(3000), rand(4000), rand(3000)
        # The internal region carries a 6 kb masked inner (an earlier round's
        # element); FastLTR only sees the pair once it is deleted.
        seq = left + ltr + inner[:2000] + "R" * 6000 + inner[2000:] + self.mutate(ltr, 0.03, rng) + right
        fa = tmp_path / "chunk.fa"
        fa.write_text(">chr1__chunk1__1-%d\n%s\n" % (len(seq), seq))
        scn = tmp_path / "chunk.fastltr.scn"
        n = fl.run_chunk(binary, str(fa), str(tmp_path), str(scn), fl.DEFAULT_ARGS.split())
        rows = [r.split() for r in scn.read_text().splitlines() if not r.startswith("#")]
        assert n == len(rows) >= 1
        l5s, l3e = 3001, len(left) + 400 + 4000 + 6000 + 400
        best = min(rows, key=lambda r: abs(int(r[0]) - l5s) + abs(int(r[1]) - l3e))
        assert abs(int(best[0]) - l5s) <= 20 and abs(int(best[1]) - l3e) <= 20
        assert int(best[2]) > 10000  # spans the masked inner in chunk coordinates

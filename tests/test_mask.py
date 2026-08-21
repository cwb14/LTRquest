"""Masking: the step that turns round N's hits into round N+1's genome."""

from __future__ import annotations

import gzip
import subprocess
import sys

import pytest

from ltrquest.mask import (
    clip_intervals,
    fasta_iter,
    mask_contig,
    merge_intervals,
    parse_features_fasta_headers_only,
    wrap_fasta,
)


class TestMergeIntervals:
    def test_empty(self):
        assert merge_intervals([]) == []

    def test_disjoint_intervals_are_kept_apart(self):
        assert merge_intervals([(0, 10), (20, 30)]) == [(0, 10), (20, 30)]

    def test_overlapping_intervals_are_merged(self):
        assert merge_intervals([(0, 10), (5, 20)]) == [(0, 20)]

    def test_adjacent_intervals_are_merged(self):
        # Half-open, so [0,10) and [10,20) touch without overlapping.
        assert merge_intervals([(0, 10), (10, 20)]) == [(0, 20)]

    def test_unsorted_input_is_handled(self):
        assert merge_intervals([(20, 30), (0, 10), (8, 22)]) == [(0, 30)]

    def test_nested_interval_does_not_shrink_its_container(self):
        assert merge_intervals([(0, 100), (10, 20)]) == [(0, 100)]


class TestClipIntervals:
    def test_interval_past_the_end_is_truncated(self):
        assert clip_intervals([(90, 200)], 100) == [(90, 100)]

    def test_negative_start_is_clamped(self):
        assert clip_intervals([(-5, 10)], 100) == [(0, 10)]

    def test_interval_entirely_outside_is_dropped(self):
        assert clip_intervals([(200, 300)], 100) == []

    def test_empty_interval_is_dropped(self):
        assert clip_intervals([(50, 50)], 100) == []


class TestWrapFasta:
    def test_wraps_at_the_requested_width(self):
        assert wrap_fasta("ACGTACGTAC", 4) == "ACGT\nACGT\nAC\n"

    def test_exact_multiple_has_no_trailing_blank_line(self):
        assert wrap_fasta("ACGTACGT", 4) == "ACGT\nACGT\n"

    def test_non_positive_width_emits_a_single_line(self):
        assert wrap_fasta("ACGT", 0) == "ACGT\n"


class TestParseFeatureHeaders:
    def test_coordinates_become_zero_based_half_open(self, tmp_path):
        # A 1-based inclusive header 11-20 covers 10 bases -> [10, 20).
        fa = tmp_path / "feats.fa"
        fa.write_text(">chr1:11-20#LTR/Copia/Ale\nACGT\n")
        assert parse_features_fasta_headers_only(str(fa)) == {"chr1": [(10, 20)]}

    def test_reversed_coordinates_are_normalised(self, tmp_path):
        fa = tmp_path / "feats.fa"
        fa.write_text(">chr1:20-11#LTR/Copia/Ale\nACGT\n")
        assert parse_features_fasta_headers_only(str(fa)) == {"chr1": [(10, 20)]}

    def test_headers_without_coordinates_are_skipped(self, tmp_path):
        fa = tmp_path / "feats.fa"
        fa.write_text(">not_a_coordinate\nACGT\n>chr1:1-5#LTR\nACGT\n")
        assert parse_features_fasta_headers_only(str(fa)) == {"chr1": [(0, 5)]}

    def test_several_features_on_one_contig(self, tmp_path):
        fa = tmp_path / "feats.fa"
        fa.write_text(">chr1:1-5#LTR\nA\n>chr1:11-15#LTR\nA\n>chr2:1-3#LTR\nA\n")
        parsed = parse_features_fasta_headers_only(str(fa))
        assert parsed == {"chr1": [(0, 5), (10, 15)], "chr2": [(0, 3)]}


class TestMaskContig:
    """feature > keep > far, in that order of priority."""

    SEQ = "A" * 100

    def test_feature_region_takes_the_feature_char(self):
        out = mask_contig(self.SEQ, [(10, 20)], "N", "V", dist=30)
        assert out[10:20] == "N" * 10

    def test_bases_within_dist_of_a_feature_survive_unmasked(self):
        out = mask_contig(self.SEQ, [(40, 50)], "N", "V", dist=20)
        assert out[20:40] == "A" * 20
        assert out[50:70] == "A" * 20

    def test_bases_beyond_dist_take_the_far_char(self):
        out = mask_contig(self.SEQ, [(40, 50)], "N", "V", dist=20)
        assert out[0:20] == "V" * 20
        assert out[70:100] == "V" * 30

    def test_length_is_never_changed(self):
        out = mask_contig(self.SEQ, [(10, 20), (60, 70)], "N", "V", dist=5)
        assert len(out) == len(self.SEQ)

    def test_no_features_masks_everything_far(self):
        out = mask_contig(self.SEQ, [], "N", "V", dist=10)
        assert out == "V" * 100

    def test_empty_sequence_is_returned_unchanged(self):
        assert mask_contig("", [(0, 10)], "N", "V", dist=10) == ""

    def test_each_round_can_use_a_different_iupac_char(self):
        # This is what lets the reconciler tell round-1 nests from round-2 nests.
        for char in "NRDYSWKMBH":
            out = mask_contig(self.SEQ, [(10, 20)], char, "V", dist=30)
            assert out[10:20] == char * 10


class TestMaskCli:
    def _run(self, *args):
        return subprocess.run(
            [sys.executable, "-m", "ltrquest.mask", *args],
            capture_output=True,
            text=True,
            check=True,
        )

    def test_masks_a_toy_genome_end_to_end(self, tmp_path, toy_genome):
        feats = tmp_path / "feats.fa"
        feats.write_text(">chr1:101-150#LTR/Copia/Ale\nACGT\n")

        out = self._run(
            "--features-fasta", str(feats),
            "--genome", str(toy_genome),
            "--feature-character", "N",
            "--far-character", "V",
            "--distance", "25",
        ).stdout

        seqs = dict(fasta_iter(out.splitlines(keepends=True)))
        assert set(seqs) == {"chr1", "chr2"}
        assert len(seqs["chr1"]) == 300
        assert seqs["chr1"][100:150] == "N" * 50          # the feature itself
        assert seqs["chr1"][75:100] == "A" * 25           # kept 25 bp flank
        assert seqs["chr1"][150:175] == "A" * 25
        assert seqs["chr1"][:75] == "V" * 75              # beyond --distance
        assert seqs["chr2"] == "V" * 120                  # no features at all

    def test_streams_a_real_chromosome_slice(self, tmp_path, athal_slice):
        original = dict(fasta_iter(athal_slice.read_text().splitlines(keepends=True)))
        (seqid,) = original

        feats = tmp_path / "feats.fa"
        feats.write_text(f">{seqid}:1001-2000#LTR/Gypsy/Tekay\nACGT\n")

        out = self._run(
            "--features-fasta", str(feats),
            "--genome", str(athal_slice),
            "--feature-character", "N",
            "--far-character", "V",
            "--distance", "15000",
        ).stdout

        ((name, seq),) = fasta_iter(out.splitlines(keepends=True))
        assert name == seqid
        assert len(seq) == len(original[name])
        assert seq[1000:2000] == "N" * 1000
        # The flank outside the feature keeps its real bases.
        assert seq[2000:2100] == original[name][2000:2100]

    def test_gzipped_genome_is_accepted(self, tmp_path, toy_genome):
        gz = tmp_path / "toy.fa.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write(toy_genome.read_text())
        feats = tmp_path / "feats.fa"
        feats.write_text(">chr1:1-10#LTR\nACGT\n")

        out = self._run(
            "--features-fasta", str(feats),
            "--genome", str(gz),
            "--feature-character", "N",
            "--far-character", "V",
            "--distance", "5",
        ).stdout
        assert out.startswith(">chr1")
        assert "N" * 10 in out

    def test_missing_genome_fails_loudly(self, tmp_path):
        feats = tmp_path / "feats.fa"
        feats.write_text(">chr1:1-10#LTR\nACGT\n")
        with pytest.raises(subprocess.CalledProcessError):
            self._run(
                "--features-fasta", str(feats),
                "--genome", str(tmp_path / "nope.fa"),
                "--feature-character", "N",
                "--far-character", "V",
                "--distance", "5",
            )

"""Strand recovery: the views, the vote, and the FASTA re-orientation."""

from __future__ import annotations

import collections
import os
import textwrap

import pytest

from ltrquest import recover_strand as rs
from ltrquest.annotate import (
    ElementInfo,
    load_recovered_strands,
    recovery_sidecar_path,
    resolve_strands,
)
from ltrquest.detect import revcomp as revcomp_record
from ltrquest.reconcile import IUPAC_DEPTH_SEQ
from ltrquest.recover_strand import (
    Element,
    Genome,
    align,
    apply_recovery,
    build_work_fastas,
    combine_views,
    depth_fasta_for,
    elements_from,
    parse_blast,
    parse_paf,
    ppt_call,
    ppt_score,
    rc,
    reorient_fasta,
    resolve_aligners,
    run_views,
    views_needed,
    votes_to_calls,
    write_sidecar,
)

# The depth-table schema, trimmed to the columns this stage reads and writes.
# nest_status stays last, which ltrquest.reconcile treats as a hard invariant.
HEADER = ("#name\tltr5_start\tltr5_end\tltr3_start\tltr3_end\t"
          "orientation\ttsd\tdomains\tnest_status")
ANNOTATED_HEADER = ("#name\tltr5_start\tltr5_end\tltr3_start\tltr3_end\t"
                    "orientation\ttsd\tstrand\tfamily\tdomains\tnest_status")


def write_depth_table(path, rows, header=HEADER):
    path.write_text(header + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")


class TestTheTwoReverseComplements:
    """`rc` is for genomic sequence; detect.revcomp is for a stored record."""

    def test_rc_complements_iupac_codes(self):
        assert rc("ACGT") == "ACGT"
        assert rc("RYKM") == "KMRY"

    def test_the_record_flip_leaves_the_depth_alphabet_alone(self):
        # Each character in IUPAC_DEPTH_SEQ names a nesting depth. Complementing
        # it would relabel the depth: R (1) would read as Y (3), and B (8) would
        # become V, which reconcile reserves for the wrapper's far-mask and does
        # not list at all.
        alphabet = "".join(IUPAC_DEPTH_SEQ)
        flipped = revcomp_record(alphabet)
        assert sorted(flipped) == sorted(alphabet)
        assert "V" not in flipped

    def test_rc_would_have_corrupted_it(self):
        # The guard this test exists for: proof the two functions really differ
        # on exactly the characters a depth FASTA carries.
        assert "V" in rc("".join(IUPAC_DEPTH_SEQ))

    def test_the_record_flip_is_its_own_inverse(self):
        record = "ACGTNRDYSWKMBH"
        assert revcomp_record(revcomp_record(record)) == record


class TestInternalBounds:
    """The cut between the LTRs, from ltr5_end/ltr3_start."""

    def element(self, start=1, end=2000, ltr5_end=300, ltr3_start=1701):
        return Element(key=f"chr1:{start}-{end}", name="chr1", chrom="chr1",
                       start=start, end=end, ltr5_end=ltr5_end,
                       ltr3_start=ltr3_start)

    def test_cuts_strictly_between_the_two_ltr_features(self):
        # 1-based positions 301..1700, i.e. ltr3_start - ltr5_end - 1 bases,
        # which is the internal length plot_summary reports from the same pair.
        assert self.element().internal_bounds() == (300, 1700)

    def test_length_matches_plot_summarys_formula(self):
        lo, hi = self.element().internal_bounds()
        assert hi - lo == 1701 - 300 - 1

    def test_offsets_are_element_relative_not_genomic(self):
        assert self.element(start=5000, end=6999).internal_bounds() == (300, 1700)

    @pytest.mark.parametrize("ltr5_end,ltr3_start", [(0, 1701), (300, 0)])
    def test_a_missing_field_declines_the_view(self, ltr5_end, ltr3_start):
        assert self.element(ltr5_end=ltr5_end, ltr3_start=ltr3_start).internal_bounds() is None

    def test_an_ltr3_start_past_the_element_declines(self):
        # gff3 omits the right LTR feature in this case, so there is no span.
        assert self.element(end=1000, ltr3_start=1701).internal_bounds() is None

    def test_too_short_to_align_declines(self):
        assert self.element(ltr5_end=300, ltr3_start=350).internal_bounds() is None


class TestVoteParsing:
    def test_blast_counts_matching_bases_per_orientation(self):
        text = "q1\tr1\tplus\t400\nq1\tr2\tplus\t150\nq2\tr1\tminus\t900\n"
        vote = parse_blast(text)
        assert vote["q1"]["+"] == 550
        assert vote["q2"]["-"] == 900

    def test_blast_ignores_a_self_hit(self):
        assert parse_blast("q1\tq1\tplus\t400\n") == {}

    def test_paf_reads_orientation_and_residue_matches(self):
        row = "q1\t1000\t0\t900\t-\tr1\t1000\t0\t900\t850\t900\t60\n"
        assert parse_paf(row)["q1"]["-"] == 850

    def test_paf_ignores_a_self_hit(self):
        row = "q1\t1000\t0\t900\t+\tq1\t1000\t0\t900\t850\t900\t60\n"
        assert parse_paf(row) == {}


class TestVotesToCalls:
    def counter(self, **kw):
        return {"q1": collections.Counter(kw)}

    def test_a_unanimous_view_calls(self):
        assert votes_to_calls(self.counter(**{"+": 400}), 100, 1.0) == {"q1": "+"}

    def test_too_few_matching_bases_declines(self):
        assert votes_to_calls(self.counter(**{"+": 50}), 100, 1.0) == {}

    def test_any_disagreement_declines_at_full_purity(self):
        assert votes_to_calls(self.counter(**{"+": 400, "-": 1}), 100, 1.0) == {}

    def test_a_relaxed_purity_tolerates_a_minority(self):
        vote = self.counter(**{"+": 400, "-": 1})
        assert votes_to_calls(vote, 100, 0.9) == {"q1": "+"}


class TestCombineViews:
    def views(self, *calls):
        return [(f"view{i}", c) for i, c in enumerate(calls)]

    def test_contradicting_views_veto_the_locus(self):
        calls, vetoed, _short = combine_views(
            self.views({"q1": "+"}, {"q1": "-"}), need=1)
        assert calls == {}
        assert vetoed == 1

    def test_one_view_is_enough_when_need_is_one(self):
        calls, _v, short = combine_views(self.views({"q1": "+"}, {}), need=1)
        assert calls["q1"][0] == "+"
        assert short == 0

    def test_one_view_is_not_enough_when_need_is_two(self):
        calls, _vetoed, short = combine_views(self.views({"q1": "+"}, {}), need=2)
        assert calls == {}
        assert short == 1

    def test_two_agreeing_views_satisfy_need_two(self):
        calls, _v, _s = combine_views(self.views({"q1": "-"}, {"q1": "-"}), need=2)
        assert calls["q1"][0] == "-"

    def test_only_balanced_asks_for_two_agreeing_views(self):
        # The accuracy contract: balanced buys its lower disagreement rate by
        # requiring a second view, which is what separates it from sensitive.
        assert views_needed("balanced") == 2
        assert views_needed("conservative") == 1
        assert views_needed("sensitive") == 1

    def test_the_evidence_names_how_many_views_agreed(self):
        calls, _v, _s = combine_views(self.views({"q1": "-"}, {"q1": "-"}), need=2)
        assert calls["q1"][1] == "homology"
        assert calls["q1"][2] == "2 of 2 views"


class TestPolypurineTract:
    def test_a_purine_run_scores_above_a_mixed_one(self):
        assert ppt_score("A" * 20, 12, 30, 3.0, 2) > ppt_score("ACAC" * 5, 12, 30, 3.0, 2)

    def test_a_plus_strand_tract_calls_plus(self):
        # The tract sits immediately 5' of the 3' LTR, so a purine run at the
        # right end of the internal region is the plus-strand signature. The
        # 5' end is purine-rich, which makes its reverse complement -- the
        # minus hypothesis -- purine-poor.
        internal = "A" * 60 + "CTCT" * 20 + "A" * 25
        seq = "C" * 10 + internal + "C" * 10
        call, plus, minus = ppt_call(seq, 10, 10, 60, 12, 30, 3.0, 2, 6.0, 12.0)
        assert call == "+"
        assert plus > minus

    def test_a_tie_declines(self):
        # Build the 5' window as the reverse complement of the 3' one, so the
        # two hypotheses are scored on identical sequence and neither can win.
        tail = "CTCT" * 9 + "A" * 24
        internal = rc(tail) + "ACGT" * 20 + tail
        seq = "C" * 10 + internal + "C" * 10
        call, plus, minus = ppt_call(seq, 10, 10, 60, 12, 30, 3.0, 2, 6.0, 12.0)
        assert plus == minus
        assert call is None

    def test_a_region_too_short_for_two_windows_declines(self):
        seq = "L" * 10 + "ACGT" * 10 + "R" * 10
        assert ppt_call(seq, 10, 10, 60, 12, 30, 3.0, 2, 6.0, 12.0)[0] is None

    def test_nothing_clearing_the_floor_declines(self):
        internal = "TCTC" * 60
        seq = "L" * 10 + internal + "R" * 10
        assert ppt_call(seq, 10, 10, 60, 12, 30, 3.0, 2, 6.0, 99.0)[0] is None


class TestGenome:
    def test_slices_one_based_inclusive(self, tmp_path):
        path = tmp_path / "g.fa"
        path.write_text(">chr1\n" + "\n".join(textwrap.wrap("ACGT" * 25, 10)) + "\n")
        genome = Genome(str(path), str(tmp_path))
        try:
            assert genome.get("chr1", 1, 4) == "ACGT"
            assert genome.get("chr1", 5, 8) == "ACGT"
        finally:
            genome.close()

    def test_the_index_goes_in_the_work_directory(self, tmp_path):
        # The genome is the user's own input and may sit on a read-only mount.
        genome_dir = tmp_path / "ref"
        genome_dir.mkdir()
        work = tmp_path / "work"
        work.mkdir()
        path = genome_dir / "g.fa"
        path.write_text(">chr1\nACGTACGTAC\n")
        Genome(str(path), str(work)).close()
        assert not (genome_dir / "g.fa.fai").exists()
        assert (work / "g.fa.fai").exists()

    def test_an_unknown_contig_returns_nothing(self, tmp_path):
        path = tmp_path / "g.fa"
        path.write_text(">chr1\nACGTACGTAC\n")
        genome = Genome(str(path), str(tmp_path))
        try:
            assert genome.get("chrZ", 1, 4) is None
        finally:
            genome.close()


class TestSidecarRoundTrip:
    def test_calls_and_declines_both_survive(self, tmp_path):
        path = recovery_sidecar_path("run", str(tmp_path))
        write_sidecar(path, ["chr1:1-100", "chr1:200-300"],
                      {"chr1:1-100": ("-", "homology", "2 of 4 views")})

        called = load_recovered_strands("run", str(tmp_path))
        assert called == {"chr1:1-100": ("-", "homology")}

        everything = load_recovered_strands("run", str(tmp_path), called_only=False)
        assert everything["chr1:200-300"] == (".", "none")

    def test_a_missing_sidecar_is_the_normal_case(self, tmp_path):
        assert load_recovered_strands("run", str(tmp_path)) == {}

    def test_the_path_can_be_pinned_explicitly(self, tmp_path):
        # A sidecar left behind by an earlier recovered run must not be picked
        # up by a later run that asked for none, so nothing globs for it.
        pinned = tmp_path / "elsewhere.tsv"
        write_sidecar(str(pinned), ["chr1:1-100"],
                      {"chr1:1-100": ("+", "ppt", "score + 30 vs - 10")})
        assert load_recovered_strands("run", str(tmp_path),
                                      path=str(pinned)) == {"chr1:1-100": ("+", "ppt")}


class TestCascadeTierOrder:
    """Recovery fills what the first three tiers leave; it never overrules them."""

    def elements(self, *keys):
        return {k: ElementInfo(name=k, superfamily="Copia", domains=[]) for k in keys}

    def test_recovery_fills_an_unstranded_element(self):
        strand, source = resolve_strands(
            self.elements("chr1:1-100"), {}, {}, {},
            recovered={"chr1:1-100": ("-", "homology")})
        assert strand["chr1:1-100"] == "-"
        assert source["chr1:1-100"] == "homology"

    def test_recovery_does_not_override_tesorter(self):
        strand, source = resolve_strands(
            self.elements("chr1:1-100"), {"chr1:1-100": "+"}, {}, {},
            recovered={"chr1:1-100": ("-", "homology")})
        assert strand["chr1:1-100"] == "+"
        assert source["chr1:1-100"] == "tesorter"

    def test_recovery_does_not_override_pass2(self):
        strand, source = resolve_strands(
            self.elements("chr1:1-100", "chr2:1-100"),
            {"chr2:1-100": "+"},
            {"chr1:1-100": "chr2:1-100"},
            {("chr1:1-100", "chr2:1-100"): "+"},
            recovered={"chr1:1-100": ("-", "homology")})
        assert strand["chr1:1-100"] == "+"
        assert source["chr1:1-100"] == "pass2"

    def test_a_declined_locus_contributes_nothing(self):
        strand, _source = resolve_strands(
            self.elements("chr1:1-100"), {}, {}, {},
            recovered={"chr1:1-100": (".", "none")})
        assert "chr1:1-100" not in strand


class TestReorientFasta:
    # Deliberately NOT a repeat of ACGT: that is its own reverse complement, so
    # asserting on it would pass whether or not the record was flipped.
    LEFT = "AAAACCCCGGGGTTTTACAC"
    RIGHT = "AAAACCCCGG"

    def fasta(self, tmp_path, name="run_depth0_ltr.fa"):
        # Wrapped the way ltrquest.reconcile writes a depth FASTA, so a rewrite
        # that flips nothing is byte-identical.
        path = tmp_path / name
        path.write_text(f">chr1:1-20#LTR/Copia\n{self.LEFT}\n"
                        f">chr2:1-10#LTR/Gypsy\n{self.RIGHT}\n")
        return path

    def test_only_the_named_record_flips(self, tmp_path):
        path = self.fasta(tmp_path)
        assert reorient_fasta(str(path), ["chr1:1-20#LTR/Copia"]) == [
            "chr1:1-20#LTR/Copia"]
        text = path.read_text()
        assert revcomp_record(self.LEFT) in text
        assert self.LEFT not in text
        assert self.RIGHT in text

    def test_flipping_twice_returns_the_original(self, tmp_path):
        path = self.fasta(tmp_path)
        before = path.read_text()
        reorient_fasta(str(path), ["chr2:1-10#LTR/Gypsy"])
        assert path.read_text() != before
        reorient_fasta(str(path), ["chr2:1-10#LTR/Gypsy"])
        assert path.read_text() == before

    def test_naming_nothing_touches_nothing(self, tmp_path):
        path = self.fasta(tmp_path)
        before = path.read_text()
        assert reorient_fasta(str(path), []) == []
        assert path.read_text() == before

    def test_a_name_the_file_does_not_hold_is_not_reported(self, tmp_path):
        path = self.fasta(tmp_path)
        before = path.read_text()
        assert reorient_fasta(str(path), ["chr9:1-5#LTR/Copia"]) == []
        assert path.read_text() == before

    def test_a_crlf_header_still_matches(self, tmp_path):
        path = tmp_path / "run_depth0_ltr.fa"
        path.write_text(f">chr1:1-20#LTR/Copia\r\n{self.LEFT}\r\n")
        assert reorient_fasta(str(path), ["chr1:1-20#LTR/Copia"]) == [
            "chr1:1-20#LTR/Copia"]
        assert revcomp_record(self.LEFT) in path.read_text()

    def test_the_flip_preserves_the_nesting_depth_alphabet(self, tmp_path):
        # The guard against reorient_fasta reaching for this module's `rc`:
        # on pure ACGT the two agree, so only a record carrying depth codes can
        # tell them apart. Complementing would turn R (depth 1) into Y (depth 3)
        # and B (depth 8) into V, which reconcile reserves for the far-mask.
        record = "ACGT" + "".join(IUPAC_DEPTH_SEQ) + "TTTT"
        path = tmp_path / "run_depth0_ltr.fa"
        path.write_text(f">chr1:1-24#LTR/Copia\n{record}\n")
        reorient_fasta(str(path), ["chr1:1-24#LTR/Copia"])
        after = "".join(l.strip() for l in open(path) if not l.startswith(">"))
        assert after == revcomp_record(record)
        assert sorted(c for c in after if c in set(IUPAC_DEPTH_SEQ)) == \
            sorted(IUPAC_DEPTH_SEQ)
        assert "V" not in after

    def test_depth_fasta_for_pairs_a_table_with_its_records(self):
        assert depth_fasta_for("run_depth0_clean_ltr.tsv") == "run_depth0_clean_ltr.fa"


class TestApplyRecovery:
    """Phase 2 follows the strand column, so the two can never disagree."""

    def setup_run(self, tmp_path, strand="-", orientation="+"):
        rows = [["chr1:1-2000#LTR/Copia", "1", "300", "1701", "2000",
                 orientation, "TATA", strand, "fam1", ".", "."]]
        tsv = tmp_path / "run_depth0_clean_ltr.tsv"
        write_depth_table(tsv, rows, ANNOTATED_HEADER)
        fasta = tmp_path / "run_depth0_clean_ltr.fa"
        fasta.write_text(">chr1:1-2000#LTR/Copia\nACGTACGTAC\n")
        write_sidecar(recovery_sidecar_path("run", str(tmp_path)),
                      ["chr1:1-2000"],
                      {"chr1:1-2000": (strand, "homology", "2 of 4 views")})
        return tsv, fasta

    def test_a_recovered_minus_record_is_stored_in_coding_sense(self, tmp_path):
        tsv, fasta = self.setup_run(tmp_path)
        assert apply_recovery("run", str(tmp_path)) == 0
        assert fasta.read_text().splitlines()[1] == "GTACGTACGT"
        assert tsv.read_text().splitlines()[1].split("\t")[5] == "-"

    def test_running_it_twice_changes_nothing(self, tmp_path):
        tsv, fasta = self.setup_run(tmp_path)
        apply_recovery("run", str(tmp_path))
        after_once = (tsv.read_text(), fasta.read_text())
        apply_recovery("run", str(tmp_path))
        assert (tsv.read_text(), fasta.read_text()) == after_once

    def test_a_recovered_plus_record_is_left_forward(self, tmp_path):
        tsv, fasta = self.setup_run(tmp_path, strand="+")
        before = fasta.read_text()
        apply_recovery("run", str(tmp_path))
        assert fasta.read_text() == before
        assert tsv.read_text().splitlines()[1].split("\t")[5] == "+"

    def test_a_locus_the_sidecar_never_mentions_is_untouched(self, tmp_path):
        rows = [["chr9:1-50#LTR/Copia", "1", "10", "41", "50", "+", "TATA",
                 "-", "fam1", ".", "."]]
        tsv = tmp_path / "run_depth0_clean_ltr.tsv"
        write_depth_table(tsv, rows, ANNOTATED_HEADER)
        fasta = tmp_path / "run_depth0_clean_ltr.fa"
        fasta.write_text(">chr9:1-50#LTR/Copia\nACGTACGTAC\n")
        write_sidecar(recovery_sidecar_path("run", str(tmp_path)), [], {})
        before = (tsv.read_text(), fasta.read_text())
        apply_recovery("run", str(tmp_path))
        assert (tsv.read_text(), fasta.read_text()) == before

    def test_a_locus_no_longer_called_goes_back_to_forward(self, tmp_path):
        # A run that switches from sensitive to conservative declines loci an
        # earlier run called. The record has to go back the way it was found.
        rows = [["chr1:1-2000#LTR/Copia", "1", "300", "1701", "2000", "-",
                 "TATA", ".", "fam1", ".", "."]]
        tsv = tmp_path / "run_depth0_clean_ltr.tsv"
        write_depth_table(tsv, rows, ANNOTATED_HEADER)
        fasta = tmp_path / "run_depth0_clean_ltr.fa"
        fasta.write_text(">chr1:1-2000#LTR/Copia\nGTACGTACGT\n")
        write_sidecar(recovery_sidecar_path("run", str(tmp_path)),
                      ["chr1:1-2000"], {})
        apply_recovery("run", str(tmp_path))
        assert fasta.read_text().splitlines()[1] == "ACGTACGTAC"
        assert tsv.read_text().splitlines()[1].split("\t")[5] == "+"

    def test_it_refuses_to_run_before_the_annotator(self, tmp_path):
        rows = [["chr1:1-2000#LTR/Copia", "1", "300", "1701", "2000", "+",
                 "TATA", ".", "."]]
        tsv = tmp_path / "run_depth0_clean_ltr.tsv"
        write_depth_table(tsv, rows)
        fasta = tmp_path / "run_depth0_clean_ltr.fa"
        fasta.write_text(">chr1:1-2000#LTR/Copia\nACGTACGTAC\n")
        write_sidecar(recovery_sidecar_path("run", str(tmp_path)),
                      ["chr1:1-2000"],
                      {"chr1:1-2000": ("-", "homology", "2 of 4 views")})
        before = fasta.read_text()
        assert apply_recovery("run", str(tmp_path)) == 0
        assert fasta.read_text() == before

    def test_it_touches_only_the_depth_fastas(self, tmp_path):
        """The round library and the pooled clustering input are not its business.

        Both describe the state at the moment they were written, before
        annotation ran, and the reconciler documents the per-round files as
        left untouched. A widened glob here would silently rewrite them.
        """
        rows = [["chr1:1-2000#LTR/Copia", "1", "300", "1701", "2000", "+",
                 "TATA", "-", "fam1", ".", "."]]
        write_depth_table(tmp_path / "run_depth0_clean_ltr.tsv", rows,
                          ANNOTATED_HEADER)
        (tmp_path / "run_depth0_clean_ltr.fa").write_text(
            ">chr1:1-2000#LTR/Copia\nAAAACCCCGGGGTTTTACAC\n")
        bystanders = {}
        for name in ("run_r1_ltr.fa", "run_all_ltr.fa", "run_all_ltr.consensus.fa",
                     "run_fpcheck.fp_LTRs.fa", "run_r2_ltr.fa"):
            path = tmp_path / name
            path.write_text(">chr1:1-2000#LTR/Copia\nAAAACCCCGGGGTTTTACAC\n")
            bystanders[path] = path.read_text()
        write_sidecar(recovery_sidecar_path("run", str(tmp_path)),
                      ["chr1:1-2000"],
                      {"chr1:1-2000": ("-", "homology", "4 of 4 views")})

        apply_recovery("run", str(tmp_path))

        assert (tmp_path / "run_depth0_clean_ltr.fa").read_text() != \
            ">chr1:1-2000#LTR/Copia\nAAAACCCCGGGGTTTTACAC\n", "the depth FASTA should have flipped"
        for path, before in bystanders.items():
            assert path.read_text() == before, f"{path.name} was rewritten"

    def test_no_sidecar_means_nothing_to_do(self, tmp_path):
        rows = [["chr1:1-2000#LTR/Copia", "1", "300", "1701", "2000", "+",
                 "TATA", "-", "fam1", ".", "."]]
        write_depth_table(tmp_path / "run_depth0_clean_ltr.tsv", rows,
                          ANNOTATED_HEADER)
        fasta = tmp_path / "run_depth0_clean_ltr.fa"
        fasta.write_text(">chr1:1-2000#LTR/Copia\nACGTACGTAC\n")
        before = fasta.read_text()
        assert apply_recovery("run", str(tmp_path)) == 0
        assert fasta.read_text() == before


class TestAlignDegradesRatherThanAborting:
    """The stage sits after every expensive one, so a benign outcome exits 0.

    The wrapper runs under `set -euo pipefail`: a non-zero exit here would throw
    away a completed detection run. Neither of these paths reaches an aligner,
    so they can be exercised without one.
    """

    def fixture(self, tmp_path, stranded):
        rows = [["chr1:1-2000#LTR/Copia/Ale", "1", "300", "1701", "2000", "+",
                 "TATA", ".", "."],
                ["chr1:5001-7000#LTR/Copia/Ale", "1", "300", "1701", "2000", "+",
                 "TATA", ".", "."]]
        write_depth_table(tmp_path / "run_depth0_clean_ltr.tsv", rows)
        (tmp_path / "genome.fa").write_text(">chr1\n" + "ACGT" * 2000 + "\n")
        if stranded:
            work = tmp_path / "run_r1.work"
            work.mkdir()
            (work / "run.cls.tsv").write_text(
                "#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\n"
                + "".join(f"{k}\tLTR\tCopia\tAle\tyes\t+\tGAG\n"
                          for k in ("chr1:1-2000", "chr1:5001-7000")))
        return str(tmp_path / "genome.fa")

    def test_everything_already_stranded_exits_zero(self, tmp_path):
        genome = self.fixture(tmp_path, stranded=True)
        assert align("run", genome, str(tmp_path)) == 0
        sidecar = tmp_path / "run_strand_recovery.tsv"
        assert sidecar.exists()
        assert load_recovered_strands("run", str(tmp_path)) == {}

    def test_no_donors_to_transfer_from_exits_zero(self, tmp_path):
        genome = self.fixture(tmp_path, stranded=False)
        assert align("run", genome, str(tmp_path)) == 0
        # Every element is recorded as looked-at-and-declined, so a later apply
        # knows which loci it owns.
        everything = load_recovered_strands("run", str(tmp_path), called_only=False)
        assert set(everything) == {"chr1:1-2000", "chr1:5001-7000"}
        assert all(v[0] == "." for v in everything.values())

    def test_a_missing_genome_is_a_real_fault(self, tmp_path):
        self.fixture(tmp_path, stranded=True)
        assert align("run", str(tmp_path / "absent.fa"), str(tmp_path)) == 1

    def test_no_depth_tables_is_a_real_fault(self, tmp_path):
        (tmp_path / "genome.fa").write_text(">chr1\nACGT\n")
        assert align("run", str(tmp_path / "genome.fa"), str(tmp_path)) == 1


class TestBuildWorkFastas:
    """The asymmetry the whole method rests on.

    References are the already-stranded elements held in coding sense; queries
    are the unstranded ones in genome-forward orientation. Invert either and
    every call comes out backwards while nothing else in the suite notices.
    """

    def setup(self, tmp_path):
        # chr1 carries three 400 bp elements at known offsets. The body must NOT
        # be its own reverse complement, or every orientation assertion below
        # passes whether or not the code oriented anything: A^n C^n G^n T^n is
        # self-complementary, so the run lengths are deliberately unequal.
        body = ("A" * 120 + "C" * 80 + "G" * 120 + "T" * 80)
        assert revcomp_record(body) != body
        assert revcomp_record(body[50:350]) != body[50:350]
        seq = "".join(["N" * 100, body, "N" * 100, body, "N" * 100, body])
        (tmp_path / "g.fa").write_text(">chr1\n" + seq + "\n")
        starts = [101, 601, 1101]
        elements = {}
        for start in starts:
            key = f"chr1:{start}-{start + 399}"
            elements[key] = Element(key=key, name=key + "#LTR/Copia", chrom="chr1",
                                    start=start, end=start + 399,
                                    ltr5_end=50, ltr3_start=351)
        return elements, Genome(str(tmp_path / "g.fa"), str(tmp_path)), body

    def read(self, path):
        out, name = {}, None
        for line in open(path):
            if line.startswith(">"):
                name = line[1:].strip()
                out[name] = []
            else:
                out[name].append(line.strip())
        return {k: "".join(v) for k, v in out.items()}

    def test_a_plus_reference_is_stored_as_read_off_the_genome(self, tmp_path):
        elements, genome, body = self.setup(tmp_path)
        keys = list(elements)
        strand = {keys[0]: "+"}
        try:
            paths, seqs, n_ref, n_qry = build_work_fastas(
                elements, strand, genome, str(tmp_path))
        finally:
            genome.close()
        assert self.read(paths["full_ref"])[keys[0]] == body

    def test_a_minus_reference_is_stored_in_coding_sense(self, tmp_path):
        elements, genome, body = self.setup(tmp_path)
        keys = list(elements)
        strand = {keys[0]: "-"}
        try:
            paths, _seqs, _n_ref, _n_qry = build_work_fastas(
                elements, strand, genome, str(tmp_path))
        finally:
            genome.close()
        assert self.read(paths["full_ref"])[keys[0]] == rc(body)

    def test_queries_stay_in_genome_forward_orientation(self, tmp_path):
        elements, genome, body = self.setup(tmp_path)
        keys = list(elements)
        try:
            paths, seqs, _n_ref, n_qry = build_work_fastas(
                elements, {keys[0]: "+"}, genome, str(tmp_path))
        finally:
            genome.close()
        queries = self.read(paths["full_qry"])
        assert set(queries) == {keys[1], keys[2]}
        assert queries[keys[1]] == body
        assert seqs[keys[1]] == body
        assert n_qry == 2

    def test_the_stranded_and_unstranded_sets_are_disjoint(self, tmp_path):
        elements, genome, _body = self.setup(tmp_path)
        keys = list(elements)
        try:
            paths, _seqs, n_ref, n_qry = build_work_fastas(
                elements, {keys[0]: "+", keys[1]: "-"}, genome, str(tmp_path))
        finally:
            genome.close()
        assert set(self.read(paths["full_ref"])) == {keys[0], keys[1]}
        assert set(self.read(paths["full_qry"])) == {keys[2]}
        assert (n_ref, n_qry) == (2, 1)

    def test_the_internal_reference_is_oriented_the_same_way(self, tmp_path):
        elements, genome, body = self.setup(tmp_path)
        keys = list(elements)
        try:
            paths, _seqs, _n_ref, _n_qry = build_work_fastas(
                elements, {keys[0]: "-"}, genome, str(tmp_path))
        finally:
            genome.close()
        internal = body[50:350]
        assert self.read(paths["int_ref"])[keys[0]] == rc(internal)
        assert self.read(paths["int_qry"])[keys[1]] == internal

    def test_an_element_off_the_end_of_the_contig_is_skipped(self, tmp_path):
        elements, genome, _body = self.setup(tmp_path)
        elements["chrZ:1-400"] = Element(key="chrZ:1-400", name="chrZ:1-400",
                                         chrom="chrZ", start=1, end=400,
                                         ltr5_end=50, ltr3_start=351)
        keys = [k for k in elements if k != "chrZ:1-400"]
        try:
            paths, seqs, _n_ref, _n_qry = build_work_fastas(
                elements, {keys[0]: "+"}, genome, str(tmp_path))
        finally:
            genome.close()
        assert "chrZ:1-400" not in seqs
        assert "chrZ:1-400" not in self.read(paths["full_qry"])


class TestViewPlan:
    """Which views each preset runs is the accuracy contract."""

    def plan(self, monkeypatch, preset, n_int_qry=5, no_internal=False,
             minimap2="minimap2"):
        seen = []

        def fake_mm(binary, ref, qry, threads, verbose):
            seen.append(("minimap2", os.path.basename(ref)))
            return {}

        def fake_bl(ref, qry, threads, evalue, workdir, verbose):
            seen.append(("blast", os.path.basename(ref)))
            return {}

        monkeypatch.setattr(rs, "run_minimap2", fake_mm)
        monkeypatch.setattr(rs, "run_blast", fake_bl)
        paths = {n: f"/tmp/{n}.fa"
                 for n in ("full_ref", "full_qry", "int_ref", "int_qry")}
        run_views(paths, preset, minimap2, n_int_qry, no_internal, 1, "1e-5", "/tmp")
        return seen

    def test_conservative_runs_the_two_dc_megablast_views_only(self, monkeypatch):
        assert self.plan(monkeypatch, "conservative") == [
            ("blast", "full_ref.fa"), ("blast", "int_ref.fa")]

    def test_balanced_runs_all_four(self, monkeypatch):
        assert len(self.plan(monkeypatch, "balanced")) == 4

    def test_sensitive_runs_all_four(self, monkeypatch):
        assert len(self.plan(monkeypatch, "sensitive")) == 4

    def test_no_internal_drops_the_internal_views(self, monkeypatch):
        seen = self.plan(monkeypatch, "balanced", no_internal=True)
        assert [r for _t, r in seen] == ["full_ref.fa", "full_ref.fa"]

    def test_no_internal_queries_drops_them_too(self, monkeypatch):
        seen = self.plan(monkeypatch, "balanced", n_int_qry=0)
        assert [r for _t, r in seen] == ["full_ref.fa", "full_ref.fa"]

    def test_a_short_run_is_reported(self, monkeypatch, capsys):
        self.plan(monkeypatch, "conservative", n_int_qry=0)
        assert "only 1 could run" in capsys.readouterr().err

    def test_a_full_run_is_not_reported(self, monkeypatch, capsys):
        self.plan(monkeypatch, "balanced")
        assert "could run here" not in capsys.readouterr().err


class TestResolveAligners:
    def test_a_missing_blastn_is_fatal(self, monkeypatch):
        monkeypatch.setattr(rs.shutil, "which", lambda name: None)
        with pytest.raises(SystemExit):
            resolve_aligners("conservative", "minimap2")

    def test_conservative_does_not_need_minimap2(self, monkeypatch):
        monkeypatch.setattr(rs.shutil, "which",
                            lambda name: "/bin/" + name if "blast" in name else None)
        assert resolve_aligners("conservative", "minimap2") is None

    def test_balanced_does_need_minimap2(self, monkeypatch):
        monkeypatch.setattr(rs.shutil, "which",
                            lambda name: "/bin/" + name if "blast" in name else None)
        with pytest.raises(SystemExit):
            resolve_aligners("balanced", "minimap2")


class TestElementsFrom:
    def load(self, tmp_path, rows, header=HEADER):
        path = tmp_path / "run_depth0_ltr.tsv"
        write_depth_table(path, rows, header)
        from ltrquest.annotate import DepthTable, load_unannotated
        table = DepthTable(str(path), 0, "raw")
        return elements_from([(table,) + load_unannotated(str(path))])

    def test_the_key_splits_on_the_last_colon(self, tmp_path):
        rows = [["scaffold:7:100-200#LTR/Copia", "1", "20", "181", "200", "+",
                 "TATA", ".", "."]]
        el = self.load(tmp_path, rows)["scaffold:7:100-200"]
        assert (el.chrom, el.start, el.end) == ("scaffold:7", 100, 200)

    def test_a_missing_ltr_field_becomes_zero_not_a_crash(self, tmp_path):
        rows = [["chr1:1-2000#LTR/Copia", "1", "NA", ".", "2000", "+",
                 "TATA", ".", "."]]
        el = self.load(tmp_path, rows)["chr1:1-2000"]
        assert (el.ltr5_end, el.ltr3_start) == (0, 0)
        assert el.internal_bounds() is None

    def test_the_full_name_is_kept_for_the_fasta_header(self, tmp_path):
        rows = [["chr1:1-2000#LTR/Copia/Ale", "1", "300", "1701", "2000", "+",
                 "TATA", ".", "."]]
        assert self.load(tmp_path, rows)["chr1:1-2000"].name == "chr1:1-2000#LTR/Copia/Ale"


class TestApplyAcrossSeveralTables:
    """The raw and clean copies of a depth are judged and rewritten separately.

    Pooling the decision would let an interrupted run -- one pair rewritten, the
    next not -- flip the finished pair a second time on the retry.
    """

    SEQ = "AAAACCCCGGGGTTTTACAC"

    def setup_pair(self, tmp_path, clean_orientation, raw_orientation):
        for variant, orientation in (("clean_", clean_orientation),
                                     ("", raw_orientation)):
            rows = [["chr1:1-2000#LTR/Copia", "1", "300", "1701", "2000",
                     orientation, "TATA", "-", "fam1", ".", "."]]
            write_depth_table(tmp_path / f"run_depth0_{variant}ltr.tsv", rows,
                              ANNOTATED_HEADER)
            seq = self.SEQ if orientation == "+" else revcomp_record(self.SEQ)
            (tmp_path / f"run_depth0_{variant}ltr.fa").write_text(
                f">chr1:1-2000#LTR/Copia\n{seq}\n")
        write_sidecar(recovery_sidecar_path("run", str(tmp_path)),
                      ["chr1:1-2000"],
                      {"chr1:1-2000": ("-", "homology", "4 of 4 views")})

    def read_seq(self, path):
        return "".join(l.strip() for l in open(path) if not l.startswith(">"))

    def test_both_copies_are_flipped(self, tmp_path):
        self.setup_pair(tmp_path, "+", "+")
        apply_recovery("run", str(tmp_path))
        for variant in ("clean_", ""):
            assert self.read_seq(tmp_path / f"run_depth0_{variant}ltr.fa") == \
                revcomp_record(self.SEQ)

    def test_a_half_finished_run_completes_rather_than_double_flipping(self, tmp_path):
        # Exactly the state an apply killed between the two pairs leaves behind.
        self.setup_pair(tmp_path, clean_orientation="-", raw_orientation="+")
        apply_recovery("run", str(tmp_path))
        clean = self.read_seq(tmp_path / "run_depth0_clean_ltr.fa")
        raw = self.read_seq(tmp_path / "run_depth0_ltr.fa")
        assert clean == revcomp_record(self.SEQ), "the finished copy was flipped again"
        assert raw == revcomp_record(self.SEQ), "the unfinished copy was not flipped"

    def test_a_table_without_an_orientation_column_is_declined(self, tmp_path):
        header = ("#name\tltr5_start\tltr5_end\tltr3_start\tltr3_end\t"
                  "tsd\tstrand\tfamily\tdomains\tnest_status")
        rows = [["chr1:1-2000#LTR/Copia", "1", "300", "1701", "2000", "TATA",
                 "-", "fam1", ".", "."]]
        write_depth_table(tmp_path / "run_depth0_clean_ltr.tsv", rows, header)
        fasta = tmp_path / "run_depth0_clean_ltr.fa"
        fasta.write_text(f">chr1:1-2000#LTR/Copia\n{self.SEQ}\n")
        write_sidecar(recovery_sidecar_path("run", str(tmp_path)),
                      ["chr1:1-2000"],
                      {"chr1:1-2000": ("-", "homology", "4 of 4 views")})
        before = fasta.read_text()
        assert apply_recovery("run", str(tmp_path)) == 0
        # Flipping without being able to record it would flip again next run.
        assert fasta.read_text() == before

    def test_a_table_whose_fasta_is_missing_is_declined(self, tmp_path):
        rows = [["chr1:1-2000#LTR/Copia", "1", "300", "1701", "2000", "+",
                 "TATA", "-", "fam1", ".", "."]]
        tsv = tmp_path / "run_depth0_clean_ltr.tsv"
        write_depth_table(tsv, rows, ANNOTATED_HEADER)
        write_sidecar(recovery_sidecar_path("run", str(tmp_path)),
                      ["chr1:1-2000"],
                      {"chr1:1-2000": ("-", "homology", "4 of 4 views")})
        assert apply_recovery("run", str(tmp_path)) == 0
        # The column must not claim a flip that never happened.
        assert tsv.read_text().splitlines()[1].split("\t")[5] == "+"

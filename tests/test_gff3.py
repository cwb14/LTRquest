"""GFF3 emission: escaping, attribute rendering, and record ordering."""

from __future__ import annotations

import gzip

import pytest

from ltrquest.gff3 import (
    SOURCE,
    Raw,
    SeqidRanker,
    build_element_blocks,
    escape,
    format_clade_composition,
    gff_line,
    is_gzip,
    order_block_lines,
    read_seq_lengths,
    render_attributes,
)


class TestEscape:
    @pytest.mark.parametrize(
        "raw,encoded",
        [
            (";", "%3B"),
            ("=", "%3D"),
            ("%", "%25"),
            ("&", "%26"),
            (",", "%2C"),
            ("\t", "%09"),
            ("\n", "%0A"),
            ("\r", "%0D"),
        ],
    )
    def test_column_nine_metacharacters_are_percent_encoded(self, raw, encoded):
        assert escape(raw) == encoded

    def test_ordinary_text_is_untouched(self):
        assert escape("LTR/Gypsy/Tekay") == "LTR/Gypsy/Tekay"

    def test_percent_is_escaped_so_the_encoding_stays_reversible(self):
        # Without this, '%3B' in a source value would decode back to ';'.
        assert escape("a%3Bb") == "a%253Bb"

    def test_numbers_are_stringified(self):
        assert escape(438418) == "438418"

    def test_raw_values_bypass_escaping(self):
        # Used where ',' is GFF3's own multi-value separator rather than data.
        assert escape(Raw("50%_SIRE,30%_Tekay")) == "50%_SIRE,30%_Tekay"


class TestRenderAttributes:
    def test_pairs_are_joined_with_semicolons(self):
        assert render_attributes([("ID", "x"), ("Name", "y")]) == "ID=x;Name=y"

    def test_unknown_values_are_omitted_entirely(self):
        # Absence is GFF3's encoding for 'not known'; 'tsd=.' says nothing.
        assert render_attributes([("ID", "x"), ("tsd", "."), ("family", "")]) == "ID=x"

    def test_none_is_omitted(self):
        assert render_attributes([("ID", "x"), ("strand_source", None)]) == "ID=x"

    def test_values_are_escaped(self):
        assert render_attributes([("Name", "a;b")]) == "Name=a%3Bb"

    def test_keys_are_escaped_too(self):
        assert render_attributes([("a=b", "c")]) == "a%3Db=c"

    def test_everything_absent_gives_an_empty_string(self):
        assert render_attributes([("tsd", "."), ("family", "")]) == ""


class TestGffLine:
    def test_the_source_column_names_this_tool(self):
        line = gff_line("chr1", "LTR_retrotransposon", 100, 200, "+", "ID=x")
        assert line.split("\t")[1] == SOURCE == "LTRquest"

    def test_all_nine_columns_are_present(self):
        line = gff_line("chr1", "LTR_retrotransposon", 100, 200, "+", "ID=x")
        cols = line.split("\t")
        assert len(cols) == 9
        assert cols[0] == "chr1"
        assert cols[3:5] == ["100", "200"]
        assert cols[5] == "."   # score
        assert cols[6] == "+"   # strand
        assert cols[7] == "."   # phase


class TestFormatCladeComposition:
    def test_percentages_are_of_family_size(self):
        clades = [("SIRE", 5), ("Tekay", 3), ("CRM", 2)]
        assert format_clade_composition(clades, 10) == "50%_SIRE,30%_Tekay,20%_CRM"

    def test_rounding_is_half_up(self):
        # 1/3 -> 33.33 -> 33; 2/3 -> 66.67 -> 67.
        assert format_clade_composition([("A", 2), ("B", 1)], 3) == "67%_A,33%_B"

    def test_a_clade_rounding_to_zero_is_dropped(self):
        out = format_clade_composition([("Big", 299), ("Rare", 1)], 300)
        assert out == "100%_Big"

    def test_output_is_sorted_by_percentage_then_name(self):
        clades = [("Zeta", 1), ("Alpha", 1), ("Big", 8)]
        assert format_clade_composition(clades, 10) == "80%_Big,10%_Alpha,10%_Zeta"

    def test_an_empty_family_yields_nothing(self):
        assert format_clade_composition([], 0) == ""
        assert format_clade_composition([("A", 1)], 0) == ""


class TestSeqidRanker:
    def test_known_sequences_keep_genome_order(self):
        ranker = SeqidRanker(["chr1", "chr2", "chr3"])
        assert [ranker.rank(c) for c in ("chr3", "chr1", "chr2")] == [2, 0, 1]

    def test_unknown_sequences_are_appended_in_first_seen_order(self):
        ranker = SeqidRanker(["chr1"])
        assert ranker.rank("scaffold_9") == 1
        assert ranker.rank("scaffold_4") == 2
        assert ranker.rank("scaffold_9") == 1   # stable on re-query

    def test_with_no_genome_order_everything_is_first_seen(self):
        ranker = SeqidRanker()
        assert ranker.rank("b") == 0
        assert ranker.rank("a") == 1


class TestOrderBlockLines:
    def _line(self, start, end):
        return gff_line("chr1", "long_terminal_repeat", start, end, "+", "ID=x")

    def test_the_parent_line_stays_first(self):
        parent = self._line(100, 900)
        child = self._line(100, 200)
        assert order_block_lines([parent, child])[0] == parent

    def test_children_are_sorted_by_position(self):
        parent = self._line(100, 900)
        lines = order_block_lines([parent, self._line(800, 900), self._line(100, 200)])
        starts = [line.split("\t")[3] for line in lines[1:]]
        assert starts == ["100", "800"]

    def test_a_longer_feature_precedes_a_shorter_one_at_the_same_start(self):
        parent = self._line(1, 1000)
        lines = order_block_lines([parent, self._line(100, 200), self._line(100, 500)])
        ends = [line.split("\t")[4] for line in lines[1:]]
        assert ends == ["500", "200"]

    def test_a_lone_line_is_returned_as_is(self):
        parent = self._line(100, 900)
        assert order_block_lines([parent]) == [parent]

    def test_no_lines_is_fine(self):
        assert order_block_lines([]) == []


class TestReadSeqLengths:
    def test_lengths_from_a_plain_fasta(self, toy_genome):
        assert read_seq_lengths(str(toy_genome)) == {"chr1": 300, "chr2": 120}

    def test_lengths_from_a_gzipped_fasta(self, tmp_path, toy_genome):
        gz = tmp_path / "toy.fa.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write(toy_genome.read_text())
        assert read_seq_lengths(str(gz)) == {"chr1": 300, "chr2": 120}

    def test_gzip_is_detected_by_magic_bytes_not_extension(self, tmp_path, toy_genome):
        # The FP-masked genome is always named '.fa' whether or not it is
        # compressed, so the extension cannot be trusted.
        lying = tmp_path / "actually_gzipped.fa"
        with gzip.open(lying, "wt") as fh:
            fh.write(toy_genome.read_text())
        assert is_gzip(str(lying)) is True
        assert read_seq_lengths(str(lying)) == {"chr1": 300, "chr2": 120}

    def test_a_plain_file_is_not_mistaken_for_gzip(self, toy_genome):
        assert is_gzip(str(toy_genome)) is False

    def test_a_missing_file_is_not_gzip(self, tmp_path):
        assert is_gzip(str(tmp_path / "nope.fa")) is False

    def test_description_after_the_name_is_ignored(self, tmp_path):
        fa = tmp_path / "desc.fa"
        fa.write_text(">NC_003071.7 Arabidopsis thaliana chromosome 2\nACGTACGTAC\n")
        assert read_seq_lengths(str(fa)) == {"NC_003071.7": 10}

    def test_wrapped_sequence_lines_are_summed(self, tmp_path):
        fa = tmp_path / "wrapped.fa"
        fa.write_text(">chr1\n" + "ACGT\n" * 25)
        assert read_seq_lengths(str(fa)) == {"chr1": 100}

    def test_a_real_chromosome_slice(self, athal_slice):
        lengths = read_seq_lengths(str(athal_slice))
        assert list(lengths.values()) == [200_000]


class TestBuildElementBlocksColumnLookups:
    """Pins every renamed element-table column lookup in build_element_blocks,
    so a typo in a column name breaks this test rather than silently emitting
    '.' in the GFF3."""

    def test_renamed_columns_reach_the_emitted_attributes(self, tmp_path):
        from ltrquest.annotate import DepthTable
        from ltrquest.detect import DETECT_COLUMNS

        row = ["."] * len(DETECT_COLUMNS)
        row[DETECT_COLUMNS.index("seq_id")] = "chr1:100-600#LTR/Copia/Ale"
        row[DETECT_COLUMNS.index("ltr5_len")] = "340"
        row[DETECT_COLUMNS.index("ltr3_len")] = "338"
        row[DETECT_COLUMNS.index("n_ts")] = "5"
        row[DETECT_COLUMNS.index("n_tv")] = "2"
        row[DETECT_COLUMNS.index("k2p")] = "0.0521"
        row[DETECT_COLUMNS.index("k2p_time")] = "868333"
        row[DETECT_COLUMNS.index("tsd")] = "TGCAA"

        path = tmp_path / "depth0_ltr.tsv"
        path.write_text("#" + "\t".join(DETECT_COLUMNS) + "\n" + "\t".join(row) + "\n")
        table = DepthTable(str(path), depth=0, variant="raw")

        blocks, skipped = build_element_blocks(
            "Athal", [table], SeqidRanker(), {}, ({}, {}))

        assert skipped == 0
        assert len(blocks) == 1
        attrs = blocks[0].payload[0]
        assert "ltr5_len=340" in attrs
        assert "ltr3_len=338" in attrs
        assert "subs=7" in attrs  # n_ts(5) + n_tv(2)
        assert "ti=5" in attrs
        assert "tv=2" in attrs
        assert "K2P_d=0.0521" in attrs
        assert "K2P_T=868333" in attrs
        assert "tsd=TGCAA" in attrs

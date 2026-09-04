"""Annotation: the `strand` and `family` columns added to every depth table."""

from __future__ import annotations

import pytest

from ltrquest.annotate import (
    UNCLASSIFIED_CLADE,
    UNKNOWN,
    DepthTable,
    Family,
    clade_composition,
    clade_of,
    element_key,
    infer_strand_from_domains,
    parse_domains_field,
    read_table,
    superfamily_of,
    write_annotated_table,
    write_table,
)


class TestElementKey:
    def test_strips_the_classification_suffix(self):
        assert element_key("chr1:100-200#LTR/Copia/Ale") == "chr1:100-200"

    def test_works_without_a_suffix(self):
        assert element_key("chr1:100-200") == "chr1:100-200"

    def test_tolerates_colons_in_the_sequence_name(self):
        assert element_key("scaffold:1:100-200#LTR") == "scaffold:1:100-200"

    def test_rejects_something_that_is_not_a_coordinate(self):
        assert element_key("just_a_name") is None

    def test_rejects_a_non_numeric_range(self):
        assert element_key("chr1:start-end") is None


class TestClassificationFields:
    @pytest.mark.parametrize(
        "name,expected",
        [
            ("chr1:1-2#LTR/Copia/Ale", "Copia"),
            ("chr1:1-2#LTR/Gypsy/Tekay", "Gypsy"),
            ("chr1:1-2#LTR", ""),
            ("chr1:1-2", ""),
        ],
    )
    def test_superfamily(self, name, expected):
        assert superfamily_of(name) == expected

    @pytest.mark.parametrize(
        "name,expected",
        [
            ("chr1:1-2#LTR/Copia/Ale", "Ale"),
            ("chr1:1-2#LTR/Gypsy/Tekay", "Tekay"),
            ("chr1:1-2#LTR/Copia", UNCLASSIFIED_CLADE),
            ("chr1:1-2#LTR/Copia/", UNCLASSIFIED_CLADE),
            ("chr1:1-2", UNCLASSIFIED_CLADE),
        ],
    )
    def test_clade_falls_back_rather_than_vanishing(self, name, expected):
        # An unclassified member must still count toward its family's
        # denominator, or clade percentages come out inflated.
        assert clade_of(name) == expected


class TestParseDomainsField:
    def test_parses_a_single_domain(self):
        assert parse_domains_field("RT|Bianca@910695-911483") == [("RT", 910695, 911483)]

    def test_parses_several_domains(self):
        field = "RH|Bianca@910002-910364;RT|Bianca@910695-911483;INT|Bianca@912348-912953"
        assert parse_domains_field(field) == [
            ("RH", 910002, 910364),
            ("RT", 910695, 911483),
            ("INT", 912348, 912953),
        ]

    def test_a_clade_less_domain_still_parses(self):
        assert parse_domains_field("RT|@100-200") == [("RT", 100, 200)]

    @pytest.mark.parametrize("field", [UNKNOWN, "", "some other column", "RT|Ale@notnum"])
    def test_a_field_that_is_not_domains_yields_nothing(self, field):
        # Callers may hand over a mis-indexed column; returning [] beats
        # returning half-parsed garbage.
        assert parse_domains_field(field) == []

    def test_one_bad_token_rejects_the_whole_field(self):
        assert parse_domains_field("RT|Ale@100-200;garbage") == []


class TestInferStrandFromDomains:
    """Pairwise concordance vote against the superfamily's canonical order."""

    COPIA_FORWARD = [("GAG", 100, 150), ("PROT", 200, 250), ("INT", 300, 350),
                     ("RT", 400, 450), ("RH", 500, 550)]

    def test_canonical_order_is_plus(self):
        assert infer_strand_from_domains(self.COPIA_FORWARD, "Copia") == "+"

    def test_reversed_order_is_minus(self):
        reversed_positions = [
            (gene, 1000 - start, 1000 - start + 50)
            for gene, start, _ in self.COPIA_FORWARD
        ]
        assert infer_strand_from_domains(reversed_positions, "Gypsy") == "-"

    def test_superfamily_changes_the_answer(self):
        # Copia puts INT before RT; Gypsy puts it last. Same layout, opposite call.
        domains = [("INT", 100, 150), ("RT", 200, 250)]
        assert infer_strand_from_domains(domains, "Copia") == "+"
        assert infer_strand_from_domains(domains, "Gypsy") == "-"

    def test_case_of_the_superfamily_does_not_matter(self):
        domains = [("INT", 100, 150), ("RT", 200, 250)]
        assert infer_strand_from_domains(domains, "gypsy") == "-"
        assert infer_strand_from_domains(domains, "GYPSY") == "-"

    def test_one_domain_is_not_enough_evidence(self):
        assert infer_strand_from_domains([("RT", 100, 200)], "Copia") == UNKNOWN

    def test_no_domains_is_unknown(self):
        assert infer_strand_from_domains([], "Copia") == UNKNOWN

    def test_unranked_domains_are_ignored(self):
        assert infer_strand_from_domains([("AP", 100, 150), ("aRH", 200, 250)],
                                         "Copia") == UNKNOWN

    def test_identical_positions_cast_no_vote(self):
        assert infer_strand_from_domains([("GAG", 100, 150), ("RT", 100, 150)],
                                         "Copia") == UNKNOWN

    def test_an_unknown_superfamily_uses_the_core_table(self):
        # INT is absent from the core ranking, leaving RT alone and unrankable.
        assert infer_strand_from_domains([("INT", 100, 150), ("RT", 200, 250)],
                                         "") == UNKNOWN
        assert infer_strand_from_domains([("GAG", 100, 150), ("RT", 200, 250)],
                                         "") == "+"

    def test_the_five_prime_most_copy_of_a_domain_wins(self):
        domains = [("RT", 900, 950), ("RT", 100, 150), ("GAG", 500, 550)]
        # RT is taken at 100, so GAG at 500 comes after it -> anti-canonical.
        assert infer_strand_from_domains(domains, "Copia") == "-"


class TestCladeComposition:
    def test_counts_are_sorted_most_abundant_first(self):
        members = [
            "chr1:1-2#LTR/Gypsy/Tekay",
            "chr1:3-4#LTR/Gypsy/Tekay",
            "chr1:5-6#LTR/Gypsy/Tekay",
            "chr1:7-8#LTR/Copia/Ale",
        ]
        assert clade_composition(members) == (("Tekay", 3), ("Ale", 1))

    def test_ties_break_alphabetically(self):
        members = ["chr1:1-2#LTR/Copia/Sire", "chr1:3-4#LTR/Copia/Ale"]
        assert clade_composition(members) == (("Ale", 1), ("Sire", 1))

    def test_unclassified_members_are_counted(self):
        members = ["chr1:1-2#LTR/Copia/Ale", "chr1:3-4"]
        assert dict(clade_composition(members)) == {"Ale": 1, UNCLASSIFIED_CLADE: 1}

    def test_no_members_gives_no_composition(self):
        assert clade_composition([]) == ()


class TestTableRoundTrip:
    def test_a_table_survives_read_then_write(self, tmp_path):
        path = tmp_path / "depth0_ltr.tsv"
        path.write_text(
            "#name\tLTR_len\ttsd\n"
            "chr1:1-100#LTR/Copia/Ale\t340\tTGCAA\n"
            "chr1:200-300#LTR/Gypsy/Tekay\t512\t.\n"
        )
        header, rows = read_table(str(path))
        assert header == ["#name", "LTR_len", "tsd"]
        assert len(rows) == 2
        assert rows[0][0] == "chr1:1-100#LTR/Copia/Ale"

        out = tmp_path / "again.tsv"
        write_table(str(out), header, rows)
        assert read_table(str(out)) == (header, rows)

    def test_a_headerless_table_is_still_readable(self, tmp_path):
        path = tmp_path / "depth0_ltr.tsv"
        path.write_text("chr1:1-100#LTR/Copia/Ale\t340\tTGCAA\n")
        header, rows = read_table(str(path))
        assert header is None
        assert rows == [["chr1:1-100#LTR/Copia/Ale", "340", "TGCAA"]]

    def test_an_empty_table_yields_no_rows(self, tmp_path):
        path = tmp_path / "depth0_ltr.tsv"
        path.write_text("")
        _header, rows = read_table(str(path))
        assert rows == []


class TestWriteAnnotatedTable:
    def test_headerless_rows_get_strand_and_family_before_the_trailing_pair(self, tmp_path):
        # No header means no schema to find `domains` in, but nest_status must
        # still end up last -- the same invariant the with-header path keeps.
        path = tmp_path / "depth0_ltr.tsv"
        path.write_text("")
        table = DepthTable(str(path), depth=0, variant="raw")
        rows = [["chr1:1-100#LTR/Copia/Ale", "RT|Ale@10-20", "."]]
        family = Family("fam00001", "chr1:1-100#LTR/Copia/Ale", 1)

        write_annotated_table(
            table, None, rows,
            strand={"chr1:1-100": "+"},
            family_by_name={"chr1:1-100#LTR/Copia/Ale": family},
            family_by_coord={},
        )

        header, written_rows = read_table(str(path))
        assert header is None
        assert written_rows == [
            ["chr1:1-100#LTR/Copia/Ale", "+", "fam00001", "RT|Ale@10-20", "."]
        ]


def test_strand_and_family_land_before_domains():
    from ltrquest.annotate import FAMILY_COL, STRAND_COL, annotation_insert_index
    from ltrquest.detect import DETECT_COLUMNS, ELEMENT_COLUMNS

    i = annotation_insert_index(DETECT_COLUMNS)
    assert DETECT_COLUMNS[i] == "domains"
    assert DETECT_COLUMNS[i - 1] == "tsd_input"

    widened = list(DETECT_COLUMNS)
    widened[i:i] = [STRAND_COL, FAMILY_COL]
    assert widened == ELEMENT_COLUMNS


def test_annotation_insert_index_falls_back_to_the_end():
    from ltrquest.annotate import annotation_insert_index

    assert annotation_insert_index(["name", "tsd"]) == 2

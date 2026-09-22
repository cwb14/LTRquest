"""The synthetic run is only useful if its truth is what it claims."""

from __future__ import annotations

import pytest

from ltrquest.detect import revcomp as revcomp_record

from reboundary_fixtures import CHROM, FAMILY, build


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build(tmp_path_factory.mktemp("syn"))


def contig(fx):
    lines = fx.genome.read_text().splitlines()
    assert lines[0] == f">{CHROM}"
    return "".join(lines[1:])


def test_every_family_copy_but_the_decayed_one_has_a_tsd_at_its_true_ends(fx):
    seq = contig(fx)
    for e in fx.family_members():
        left = seq[e.true_start - 6:e.true_start - 1]
        right = seq[e.true_end:e.true_end + 5]
        assert (left == right) == (e.kind != "decayed"), e.kind


def test_calls_sit_inside_their_true_spans_and_only_truncated_ones_differ(fx):
    for e in fx.family_members():
        assert e.true_start <= e.start and e.end <= e.true_end
        truncated = e.kind not in ("normal", "decayed")
        assert ((e.start, e.end) != (e.true_start, e.true_end)) == truncated, e.kind


def test_the_intact_ltrs_start_tg_and_end_ca(fx):
    seq = contig(fx)
    e = [x for x in fx.family_members()
         if x.kind == "normal" and x.strand == "+"][0]
    assert seq[e.start - 1:e.start + 1] == "TG"
    assert seq[e.end - 2:e.end] == "CA"
    assert e.l1 - e.start + 1 == 400 and e.end - e.r0 + 1 == 400


def test_fasta_records_are_the_called_spans_stored_as_the_tables_say(fx):
    seq = contig(fx)
    records = {}
    for path in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.fa"):
        name = None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:]
                records[name] = ""
            else:
                records[name] += line
    for e in fx.elements:
        fwd = seq[e.start - 1:e.end]
        stored = records[e.name]
        if e.kind == "host":
            inner = fx.kind("nested_del_left")
            assert stored[inner.start - e.start:inner.end - e.start + 1] == (
                "N" * (inner.end - inner.start + 1))
            continue
        assert stored == (revcomp_record(fwd) if e.strand == "-" else fwd), e.kind


def test_the_tables_carry_the_33_columns_and_the_nesting(fx):
    rows = (fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv").read_text().splitlines()
    assert rows[0].startswith("#seq_id\tseq_len\tstatus")
    assert rows[0].endswith("nest_status")
    assert len(rows[0].split("\t")) == 33
    host = fx.kind("host")
    inner = fx.kind("nested_del_left")
    assert host.depth == 1
    assert inner.nest_status == f"nest-inner:{host.key}"
    assert all(e.family == FAMILY for e in fx.family_members())

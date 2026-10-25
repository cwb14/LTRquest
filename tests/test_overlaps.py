"""Overlapping calls: nested ones are both kept, any other overlap keeps the stronger."""

from __future__ import annotations

import pytest

from ltrquest import ends
from ltrquest import overlaps as ov


def call(start, end, ltr=300, rank=(0, 0.0, 0), chrom="c"):
    return ov.Call(f"{chrom}:{start}-{end}", chrom, start, end,
                   (start, start + ltr - 1, end - ltr + 1, end), rank)


def test_hallmarks_count_only_a_tsd_at_the_called_ends():
    assert ends.hallmarks("AATGC", "0,0", "tg...ca") == (True, True)
    assert ends.hallmarks("AATGC", "1,0", "tg...ca") == (False, True)    # found 1 bp off
    assert ends.hallmarks(".", "NA", "tg...ct") == (False, False)
    assert ends.hallmarks("NA", "NA", "TG...CA") == (False, True)


def test_evidence_ranks_tsd_above_tg_ca_and_both_above_either():
    both = ov.evidence("AATGC", "0,0", "tg...ca", "10", 0)
    tsd = ov.evidence("AATGC", "0,0", "ta...ca", "900", 5)
    tgca = ov.evidence(".", "NA", "tg...ca", "900", 5)
    none = ov.evidence(".", "NA", "ga...ct", "900", 5)
    assert both > tsd > tgca > none


def test_n_domains_counts_distinct_genes():
    assert ov.n_domains("RT|Ale@1-9;RH|Ale@20-30;RT|Ale@40-50") == 2
    assert ov.n_domains(".") == 0


@pytest.mark.parametrize("a, b, want", [
    ((100, 900), (1000, 1900), None),                  # disjoint
    ((100, 2000), (700, 1500), "nested"),              # inside, own LTRs
    ((100, 2000), (300, 1500), "conflict"),            # inside, but shares the left LTR
    ((100, 1500), (900, 2400), "conflict"),            # crossing
    ((100, 1500), (100, 1500), "conflict"),            # the same span twice
    ((100, 1500), (100, 1400), "conflict"),            # one end shared: the LTRs overlap
])
def test_relation(a, b, want):
    assert ov.relation(call(*a), call(*b)) == want
    assert ov.relation(call(*b), call(*a)) == want


def test_unknown_ltrs_are_taken_as_distinct():
    a = ov.Call("c:100-2000", "c", 100, 2000, None, (0,))
    b = ov.Call("c:300-1500", "c", 300, 1500, None, (0,))
    assert ov.relation(a, b) == "nested"


def test_a_crossing_pair_keeps_the_call_with_the_hallmarks():
    weak = call(100, 1500, rank=ov.evidence(".", "NA", "tg...ca", "2000", 4))
    strong = call(900, 2400, rank=ov.evidence("AATGC", "0,0", "tg...ca", "300", 0))
    kept, dropped = ov.resolve([weak, strong])
    assert kept == [strong] and dropped == [(weak, strong)]


def test_without_hallmarks_the_better_ltr_pair_wins():
    a = call(100, 1500, rank=ov.evidence(".", "NA", "ga...ct", "400", 0))
    b = call(900, 2400, rank=ov.evidence(".", "NA", "ga...ct", "550", 0))
    assert ov.resolve([a, b])[0] == [b]


def test_domains_break_a_tie_and_genomic_order_breaks_the_rest():
    a = call(100, 1500, rank=(0, 500.0, 1))
    b = call(900, 2400, rank=(0, 500.0, 3))
    assert ov.resolve([a, b])[0] == [b]
    c = call(900, 2400, rank=(0, 500.0, 1))
    assert ov.resolve([c, a])[0] == [a]


def test_a_chain_is_resolved_strongest_first():
    # a-b and b-c cross, a and c are disjoint: b is strongest, so it alone survives;
    # were b weakest, a and c would both stay.
    a, c = call(100, 1500, rank=(1,)), call(2000, 3400, rank=(1,))
    assert ov.resolve([a, call(1200, 2300, rank=(3,)), c])[0] == [call(1200, 2300, rank=(3,))]
    assert ov.resolve([a, call(1200, 2300, rank=(0,)), c])[0] == [a, c]


def test_nesting_survives_and_a_duplicate_inner_does_not():
    outer = call(100, 9000, rank=(3,))
    inner = call(3000, 5000, rank=(1,))
    twin = call(3002, 5000, rank=(0,))        # the same inner, called 2 bp short
    kept, dropped = ov.resolve([outer, inner, twin])
    assert kept == [outer, inner] and dropped == [(twin, inner)]
    assert ov.nest_pairs(kept) == [(inner, outer)]


def test_identical_spans_keep_one_copy():
    a, b = call(100, 1500, rank=(1,)), call(100, 1500, rank=(1,))
    kept, dropped = ov.resolve([a, b])
    assert len(kept) == 1 and len(dropped) == 1


def test_calls_on_other_chromosomes_never_conflict():
    a, b = call(100, 1500, chrom="c1"), call(100, 1500, chrom="c2")
    assert ov.resolve([a, b])[0] == [a, b]


def test_long_calls_spanning_index_bins_still_meet():
    a = call(50_000, 260_000, rank=(1,))
    b = call(250_000, 300_000, rank=(0,))
    assert ov.resolve([a, b])[0] == [a]

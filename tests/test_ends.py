"""Canonical ends: TG..CA termini with an exact TSD, within 2 bp of the call."""

from __future__ import annotations

import pytest

from ltrquest import ends

BODY = "TGTTAGCCGATCCATGGACTTAGCCA"      # an "element": TG ... CA, 26 bp


def element(tsd="AATGC", body=BODY, left="CCGTTAGGACTC", right="GGATCCTTAGAC"):
    """Genome `left + tsd + body + tsd + right`; returns (seq, start, end), 1-based."""
    seq = left + tsd + body + tsd + right
    start = len(left) + len(tsd) + 1
    return seq, start, start + len(body) - 1


def snap_call(seq, start, end):
    ctx, b5, b3 = ends.end_windows(seq, start, end)
    return ends.snap(ctx, b5, b3)


def test_canonical_ends_stay_put():
    seq, s, e = element()
    assert snap_call(seq, s, e) == ends.Ends(0, 0, "AATGC")


@pytest.mark.parametrize("d5, d3", [(1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (-1, 0), (0, -2),
                                    (-1, 1), (2, -2)])
def test_a_call_up_to_2_bp_off_at_either_end_snaps_back(d5, d3):
    # A call d5 bp outside the true left end has d5 = -(that) to move: the shift
    # reported is how far the called boundary must move, positive into the element.
    seq, s, e = element()
    got = snap_call(seq, s - d5, e + d3)
    assert got == ends.Ends(d5, d3, "AATGC")


def test_three_bp_off_is_not_reached():
    seq, s, e = element()
    assert snap_call(seq, s - 3, e) is None


def test_the_4_bp_redescription_of_a_wide_call_loses_to_the_canonical_ends():
    # P. supina Ps1B:542705-550914: called 1 bp wide on the left, the called ends carry
    # an exact 4-bp TSD (AATG|C...CA|AATG) but start C, not TG. The canonical ends one
    # base in carry the 5-bp TSD.
    seq, s, e = element()
    ctx, b5, b3 = ends.end_windows(seq, s - 1, e)
    assert ends.tsd_at(ctx, b5, b3) == "AATG"
    assert ends.snap(ctx, b5, b3) == ends.Ends(1, 0, "AATGC")


def test_a_real_4_bp_tsd_is_kept_over_a_chance_5_mer_at_a_shifted_end():
    # TAAA|TG...CA|TAAA, with the base after the right copy a T: moving the left end
    # in by one gives the 5-mer TAAAT on both sides, at non-canonical ends (G...CA).
    seq, s, e = element(tsd="TAAA", right="TCCTTAGACGGA")
    ctx, b5, b3 = ends.end_windows(seq, s, e)
    assert ends.tsd_at(ctx, b5 + 1, b3) == "TAAAT"
    assert ends.snap(ctx, b5, b3) == ends.Ends(0, 0, "TAAA")


def test_tg_ca_without_a_tsd_is_not_canonical():
    seq, s, e = element()
    seq = seq[:e] + "T" + seq[e + 1:]           # break the right copy of the TSD
    assert snap_call(seq, s, e) is None


@pytest.mark.parametrize("tsd", ["AAAAA", "ANNGC"])
def test_homopolymer_and_gap_duplications_are_refused(tsd):
    seq, s, e = element(tsd=tsd)
    assert snap_call(seq, s, e) is None


def test_case_is_ignored_and_the_tsd_reported_uppercase():
    seq, s, e = element()
    assert snap_call(seq.lower(), s, e) == ends.Ends(0, 0, "AATGC")


def test_an_element_at_the_contig_start_is_read_without_wrapping():
    seq, s, e = element(left="")
    assert s == 6
    assert snap_call(seq, s, e) == ends.Ends(0, 0, "AATGC")
    assert snap_call(seq, s + 1, e) == ends.Ends(-1, 0, "AATGC")


def test_the_nearest_of_two_canonical_ends_wins():
    # TGTG at the left end: TG two bases in is also a candidate, but its TSD
    # (5 bases ending ..TG) does not match, so only the called ends qualify; with a
    # matching duplication there instead, the smaller move still wins.
    seq, s, e = element(body="TGTG" + BODY[2:])
    assert snap_call(seq, s, e) == ends.Ends(0, 0, "AATGC")


@pytest.mark.parametrize("tsd, offset, motif, want", [
    ("AATGC", "0,0", "tg...ca", "high"),
    ("AATGC", "0,0", "ga...ca", "medium"),
    (".", "NA", "tg...ca", "medium"),
    ("AATGC", "1,0", "tg...ca", "medium"),      # its TSD sits at other ends
    ("NA", "NA", "ag...ct", "low"),
])
def test_confidence_counts_the_hallmarks_at_the_reported_ends(tsd, offset, motif, want):
    assert ends.confidence(tsd, offset, motif) == want


def test_tsd_between_reads_canonical_ends_through_a_fetcher():
    seq, s, e = element()

    def fetch(lo, hi):
        a, b = max(1, lo), min(len(seq), hi)
        return seq[a - 1:b], a

    assert ends.tsd_between(fetch, s, e) == "AATGC"
    assert ends.tsd_between(fetch, s - 1, e) is None          # not TG at the left end
    short, s2, e2 = element(left="")
    assert ends.tsd_between(lambda lo, hi: (short[max(1, lo) - 1:min(len(short), hi)],
                                            max(1, lo)), s2, e2) == "AATGC"

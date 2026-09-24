"""Placing trusted templates on one element: per-side placement, anchors, votes, moves."""

from __future__ import annotations

import dataclasses
import random

import pytest

from ltrquest import ltr_model as lm
from ltrquest import ltr_place as lp
from ltrquest import templates as tp
from reboundary_fixtures import FAMILY, build, members


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build(tmp_path_factory.mktemp("syn_place"))


@pytest.fixture(scope="module")
def g(fx):
    return lm.Genomes({fx.prefix: str(fx.genome)})


@pytest.fixture(scope="module")
def mem(fx):
    return members(fx)


def member(fx, mem, kind):
    e = fx.kind(kind)
    return next(m for m in mem if m.name == e.name), e


def models(fx, mem, g, target, frame="+"):
    """Every trusted copy of the family, oriented by the true strands (no BLAST)."""
    out = []
    for t in tp.trusted(mem, g):
        if t.family != FAMILY:
            continue
        rel = "+" if t.strand == frame else "-"
        out.append((t.name, *tp.models_for(t, rel, g)))
    return out


def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _genome(tmp_path, seq, name="chrT"):
    fa = tmp_path / "g.fa"
    fa.write_text(f">{name}\n{seq}\n")
    return lm.Genomes({"p": str(fa)})


def _m(start, end, l1, r0, strand="+", chrom="chrT"):
    return lm.Member(prefix="p", name=f"{chrom}:{start}-{end}#LTR/x", chrom=chrom,
                     start=start, end=end, l1=l1, r0=r0, strand=strand, orientation="+",
                     family="f", depth=0, k2p=0.01, tsd=".")


# ---------------------------------------------------------------- pieces

def test_local_reports_query_and_window_coordinates():
    q = "ACGTACGTTTGACCAGTTAG"
    w = "GGGGG" + q[2:] + "CCCC"
    h = lp.local(q, w)
    assert (h.q_start, h.q_end, h.w_start, h.w_end) == (2, len(q) - 1, 5, 5 + len(q) - 3)
    assert h.identity == 1.0


@pytest.mark.parametrize("positions, want", [
    ([], (None, 0)),
    ([100], (None, 1)),
    ([100, 101], (100, 2)),
    ([100, 101, 150], (101, 2)),
    ([100, 200, 300], (None, 1)),
    ([300, 100, 102, 99], (100, 3)),
])
def test_the_vote_is_the_low_median_and_needs_agreeing_support(positions, want):
    assert lp.vote(positions, agree=2, need=2) == want


# ---------------------------------------------------------------- the fixture's kinds

@pytest.mark.parametrize("kind, side, shift", [
    ("del_left", "left", 50), ("nested_del_left", "left", 50),
    ("unstranded_del_left", "left", 50), ("short3_left", "left", 3),
    ("patch_right", "right", 70), ("del_left_minus", "right", 50),
])
def test_a_short_call_is_moved_to_its_true_end(fx, g, mem, kind, side, shift):
    m, e = member(fx, mem, kind)
    frame = "-" if e.strand == "-" else "+"
    p = lp.propose(m, models(fx, mem, g, m, frame), g, lp.Params())
    assert (p.left, p.right) == ((e.true_start, e.end) if side == "left" else (e.start, e.true_end))
    assert (p.ext_left, p.ext_right) == ((shift, 0) if side == "left" else (0, shift))
    assert p.sup_left + p.sup_right >= 2


def test_a_minus_strand_call_reports_its_move_at_its_biological_five_prime_end(fx, g, mem):
    m, e = member(fx, mem, "del_left_minus")
    p = lp.propose(m, models(fx, mem, g, m, "-"), g, lp.Params())
    assert p.orient == "-" and (p.ext5, p.ext3) == (50, 0)


def test_a_small_over_call_is_trimmed(fx, g, mem):
    m, e = member(fx, mem, "long3_right")
    p = lp.propose(m, models(fx, mem, g, m), g, lp.Params())
    assert (p.left, p.right) == (e.start, e.true_end) and p.ext_right == -3


def test_trimming_can_be_switched_off(fx, g, mem):
    m, _ = member(fx, mem, "long3_right")
    assert lp.propose(m, models(fx, mem, g, m), g, lp.Params(max_trim=0)) is None


def test_a_large_insertion_is_crossed_by_the_terminal_anchor(fx, g, mem):
    m, e = member(fx, mem, "ins_left")
    p = lp.propose(m, models(fx, mem, g, m), g, lp.Params())
    assert p.left == e.true_start and p.left_src == "anchor"


def test_an_intact_call_is_not_a_candidate(fx, g, mem):
    m, _ = member(fx, mem, "decayed")
    assert lp.propose(m, models(fx, mem, g, m), g, lp.Params()) is None


def test_one_template_is_not_enough_support(fx, g, mem):
    m, e = member(fx, mem, "del_left")
    p = lp.propose(m, models(fx, mem, g, m)[:1], g, lp.Params())
    assert p.left == m.start and p.reason_left == "no_support" and p.raw_left == e.true_start


def test_the_identity_gate_blocks_a_weak_end(fx, g, mem):
    m, e = member(fx, mem, "del_left")
    p = lp.propose(m, models(fx, mem, g, m), g, lp.Params(min_identity=1.01))
    assert p.left == m.start and p.reason_left == "gate_identity" and p.raw_left == e.true_start


def test_the_vote_does_not_depend_on_template_order(fx, g, mem):
    m, _ = member(fx, mem, "patch_right")
    ms = models(fx, mem, g, m)
    a = lp.propose(m, ms, g, lp.Params())
    b = lp.propose(m, list(reversed(ms)), g, lp.Params())
    assert (a.left, a.right, a.sup_left, a.sup_right) == (b.left, b.right, b.sup_left, b.sup_right)


def test_a_placement_that_misses_the_call_on_its_side_is_refused(fx, g, mem):
    """The window is widened by either LTR's deficit, so it can reach a neighbouring
    copy. A placement that does not overlap the element's own call on that side is
    not a placement of that element."""
    n = next(x for x in mem if x.family == FAMILY and x.len_left == 400 and x.strand == "+")
    ghost = dataclasses.replace(n, start=n.end + 60, l1=n.end + 99, r0=n.end + 900,
                                end=n.end + 1299)
    left, right = tp.models_for(n, "+", g)
    assert lp.place_left(ghost, left, g) is None and lp.place_right(ghost, right, g) is None


def test_obstacle_names_the_deletion(fx, g, mem):
    m, _ = member(fx, mem, "del_left")
    p = lp.propose(m, models(fx, mem, g, m), g, lp.Params())
    left, right = lp.obstacle(m, p, g)
    assert left.startswith("gap:") and int(left.split(":")[1]) >= 50 and right == "."


def test_tsd_probe_reads_the_planted_duplication_and_not_the_null(fx, g, mem, k2l_api):
    m, e = member(fx, mem, "del_left")
    found = lp.tsd_at(k2l_api, g, m, e.true_start, e.true_end)
    assert found not in (".", "NA") and len(found) >= 5
    assert lp.tsd_at(k2l_api, g, m, m.start, m.end) == "."


# ---------------------------------------------------------------- synthetic edge cases

def _element(ltr, internal, outer_mut=None):
    left = (outer_mut + ltr[len(outer_mut):]) if outer_mut else ltr
    return left + internal + ltr


def test_an_anchor_never_jumps_onto_a_complete_neighbouring_copy(tmp_path):
    """The element's outer 30 bp are too mutated for the core; a full solo LTR sits
    1 kb upstream. Its first 30 bp match the model perfectly -- but a complete copy
    starting there is a neighbour, not this element's terminus."""
    ltr, internal = "TG" + _rnd(396, 1) + "CA", _rnd(1500, 2)
    mutated = "".join("A" if c != "A" else "C" for c in ltr[:30])
    head = _rnd(2000, 3)
    solo_start = len(head) + 1
    left_pad = head + ltr + _rnd(1000, 4)
    true_start = len(left_pad) + 1
    seq = left_pad + _element(ltr, internal, mutated) + _rnd(2000, 5)
    g = _genome(tmp_path, seq)
    true_end = true_start + 400 + 1500 + 400 - 1
    m = _m(true_start + 60, true_end, true_start + 399, true_end - 399 + 60)
    tpl = [(f"t{i}", ltr, ltr) for i in range(3)]
    p = lp.propose(m, tpl, g, lp.Params())
    assert p is None or p.left in (m.start, true_start)
    assert p is None or p.left != solo_start


def test_a_call_at_a_contig_edge_moves_no_further_than_the_contig(tmp_path):
    ltr, internal = "TG" + _rnd(396, 11) + "CA", _rnd(1200, 12)
    seq = _rnd(20, 13) + ltr + internal + ltr + _rnd(500, 14)
    g = _genome(tmp_path, seq)
    true_start, true_end = 21, 20 + 400 + 1200 + 400
    m = _m(true_start + 40, true_end, true_start + 399, true_end - 399 + 40)
    p = lp.propose(m, [(f"t{i}", ltr, ltr) for i in range(3)], g, lp.Params())
    assert p.left == true_start and p.left >= 1


def test_an_anchor_never_places_an_end_off_the_contig(tmp_path):
    """The contig starts 2 bp inside the element's 5' LTR, and a 1 kb insertion near
    that LTR's outer end sends the left side to the anchor: the model's outer 30 bp
    match the contig's first 28, so extrapolating the model's terminal base would put
    the end at position -1."""
    ltr = "TG" + _rnd(396, 31) + "CA"
    seq = ltr[2:100] + _rnd(1000, 32) + ltr[100:] + _rnd(1500, 33) + ltr + _rnd(800, 34)
    g = _genome(tmp_path, seq)
    m = _m(1099, 3298, 1398, 2999)
    p = lp.propose(m, [(f"t{i}", ltr, ltr) for i in range(3)], g, lp.Params())
    assert p is None or (p.left >= 1 and p.right <= len(seq))


def test_an_anchor_never_places_an_end_past_the_contig_end(tmp_path):
    """The mirror image: the contig ends 2 bp before the element's 3' end."""
    ltr = "TG" + _rnd(396, 35) + "CA"
    seq = _rnd(800, 36) + ltr + _rnd(1500, 37) + ltr[:300] + _rnd(1000, 38) + ltr[300:398]
    g = _genome(tmp_path, seq)
    m = _m(801, 3000, 1100, 2701)
    p = lp.propose(m, [(f"t{i}", ltr, ltr) for i in range(3)], g, lp.Params())
    assert p is None or (p.left >= 1 and p.right <= len(seq))


def test_bases_an_extension_would_add_across_an_assembly_gap_are_refused(tmp_path):
    ltr, internal = "TG" + _rnd(396, 21) + "CA", _rnd(1200, 22)
    gapped = ltr[:10] + "N" * 5 + ltr[15:]
    seq = _rnd(500, 23) + gapped + internal + ltr + _rnd(500, 24)
    g = _genome(tmp_path, seq)
    true_start, true_end = 501, 500 + 400 + 1200 + 400
    m = _m(true_start + 40, true_end, true_start + 399, true_end - 399 + 40)
    p = lp.propose(m, [(f"t{i}", ltr, ltr) for i in range(3)], g, lp.Params())
    assert p.left == m.start and p.reason_left == "gap_in_added"

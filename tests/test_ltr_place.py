"""Placing a family model on one element: extend only, gates, anchor, obstacles."""

from __future__ import annotations

import dataclasses

import pytest

from ltrquest import ltr_model as lm
from ltrquest import ltr_place as lp

from reboundary_fixtures import FAMILY, build, members


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build(tmp_path_factory.mktemp("syn_place"))


@pytest.fixture(scope="module")
def g(fx):
    return lm.Genomes({fx.prefix: str(fx.genome)})


@pytest.fixture(scope="module")
def truth_model(fx):
    """The planted LTR itself: placement tests should not depend on MAFFT."""
    return lm.Model(f"{FAMILY}:truth", FAMILY, fx.ltr, 400.0, 13)


def member(fx, kind):
    e = fx.kind(kind)
    return [m for m in members(fx) if m.name == e.name][0], e


def test_local_reports_query_and_window_coordinates():
    q = "ACGTACGTTTGACCAGTTAG"
    w = "GGGGG" + q[2:] + "CCCC"
    h = lp.local(q, w)
    assert (h.q_start, h.q_end, h.w_start, h.w_end) == (2, len(q) - 1, 5, 5 + len(q) - 3)
    assert h.identity == 1.0


def test_intact_copy_is_not_a_candidate(fx, g, truth_model):
    m = [x for x in members(fx) if x.name == fx.elements[0].name][0]
    assert lp.propose(m, [truth_model], g, lp.Params()) is None


def test_deletion_near_the_left_end_is_recovered_by_the_core(fx, g, truth_model):
    m, e = member(fx, "del_left")
    p = lp.propose(m, [truth_model], g, lp.Params())
    assert p.gate_ok and (p.left, p.right) == (e.true_start, e.end)
    assert p.left_src == "core" and p.right_src == "." and p.ext5 == 50 and p.ext3 == 0


def test_mutation_patch_near_the_right_end_is_recovered(fx, g, truth_model):
    m, e = member(fx, "patch_right")
    p = lp.propose(m, [truth_model], g, lp.Params())
    assert p.gate_ok and (p.left, p.right) == (e.start, e.true_end) and p.ext3 == 70


def test_minus_strand_copy_extends_at_its_biological_five_prime_end(fx, g, truth_model):
    m, e = member(fx, "del_left_minus")
    p = lp.propose(m, [truth_model], g, lp.Params())
    assert p.orient == "-" and (p.left, p.right) == (e.start, e.true_end)
    assert p.ext5 == 50 and p.ext3 == 0


def test_large_insertion_is_crossed_by_the_terminal_anchor(fx, g, truth_model):
    m, e = member(fx, "ins_left")
    p = lp.propose(m, [truth_model], g, lp.Params())
    assert p.gate_ok and p.left == e.true_start and p.left_src == "anchor"


def test_without_the_anchor_the_insertion_fails_its_gate(fx, g, truth_model):
    m, e = member(fx, "ins_left")
    p = lp.propose(m, [truth_model], g, lp.Params(anchor=False))
    assert p is None or (not p.gate_ok and p.left == m.start)


def test_never_trims_an_over_extended_call(fx, g, truth_model):
    base = [x for x in members(fx) if x.family == FAMILY and x.len_left == 400][0]
    wide = dataclasses.replace(base, start=base.start - 30, end=base.end + 30)
    assert lp.propose(wide, [truth_model], g, lp.Params()) is None


def test_identity_gate_blocks_a_weak_end(fx, g, truth_model):
    m, e = member(fx, "del_left")
    p = lp.propose(m, [truth_model], g, lp.Params(min_identity=1.01))
    assert p is not None and not p.gate_ok and p.left == m.start and p.raw_left == e.true_start


def test_obstacle_names_the_deletion(fx, g, truth_model):
    m, _ = member(fx, "del_left")
    p = lp.propose(m, [truth_model], g, lp.Params())
    left, right = lp.obstacle(m, p, g)
    assert left.startswith("gap:") and int(left.split(":")[1]) >= 50 and right == "."


def test_tsd_probe_reads_the_planted_duplication_and_not_the_null(fx, g, k2l_api):
    m, e = member(fx, "del_left")
    found = lp.tsd_at(k2l_api, g, m, e.true_start, e.true_end)
    assert found not in (".", "NA") and len(found) >= 5
    assert lp.tsd_at(k2l_api, g, m, m.start, m.end) == "."

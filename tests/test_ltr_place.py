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


def test_nearest_templates_are_references_whose_ends_the_consensus_confirms(fx, g, truth_model):
    refs, _ = lm.select_references([m for m in members(fx) if m.family == FAMILY], "modal", FAMILY)
    tpl = lp.nearest_templates(truth_model, refs, g)
    assert len(tpl) >= 10 and all(t.model_id.startswith(f"{FAMILY}:tpl:") for t in tpl)
    assert all(len(t.seq) == 400 for t in tpl)          # the truncated ref (ins_left) is not one


def test_median_combining_agrees_with_best_on_a_clean_case(fx, g, truth_model):
    refs, _ = lm.select_references([m for m in members(fx) if m.family == FAMILY], "modal", FAMILY)
    tpl = lp.nearest_templates(truth_model, refs, g)
    m, e = member(fx, "del_left")
    p = lp.propose(m, tpl, g, lp.Params(combine="median"))
    assert p.gate_ok and p.left == e.true_start


def test_an_element_is_never_a_model_for_itself(fx, g, truth_model):
    """Leave-one-out is carried by the model, not by a caller that must remember it."""
    m, _ = member(fx, "del_left")
    mine = dataclasses.replace(truth_model, model_id=f"{FAMILY}:mine",
                               refs=frozenset({m.uid}))
    assert lp.candidate_models(m, [mine, truth_model], g, 3) == [truth_model]
    assert lp.propose(m, [mine], g, lp.Params()) is None
    assert lp.propose(m, [mine, truth_model], g, lp.Params()).model_id == truth_model.model_id


def test_the_gates_describe_the_coordinate_that_is_committed(fx, g, truth_model):
    """Two models agree on an end 20 bp too far out; the third, best-scoring one does not.

    The median end is the two's, so the identity that authorises it must be theirs
    as well -- with the best model's 1.0 the element would be committed to an end
    no model ever scored.
    """
    m, e = member(fx, "del_left")
    junk = ("GATTACA" * 3)[:20]
    over = [lm.Model(f"{FAMILY}:over{i}", FAMILY, junk + fx.ltr, 400.0, 5) for i in (1, 2)]
    c_over = lp.core(m, over[0].seq, g)
    assert c_over.left < lp.core(m, truth_model.seq, g).left == e.true_start
    assert c_over.id_left < lp.Params().min_identity        # only the outer end is wrong
    assert c_over.whole >= lp.Params().min_whole_identity   # the other gate does not save it
    p = lp.propose(m, [truth_model] + over, g, lp.Params(combine="median"))
    assert p.raw_left == c_over.left                      # the median end, as before
    assert p.id_left == pytest.approx(c_over.id_left)     # scored where it was committed
    assert not p.gate_ok and p.left == m.start


def test_a_reference_whose_inner_end_is_over_called_is_not_a_template(fx, g, truth_model):
    """An inner end nothing checked sets the target's opposite outer end (6 bp clears
    the outer-30 gate and quietly destroys the TSD), so it is checked like the outer
    ends: the consensus has to put it in the same place."""
    fam = [x for x in members(fx) if x.family == FAMILY]
    refs, _ = lm.select_references(fam, "modal", FAMILY)
    over = [dataclasses.replace(r, l1=r.l1 + 6, r0=r.r0 - 6) for r in refs]
    tpl = lp.nearest_templates(truth_model, over, g)
    assert [t.model_id for t in tpl] == [truth_model.model_id]     # none qualifies
    m, e = member(fx, "del_left")
    p = lp.propose(m, tpl, g, lp.Params(combine="median"))
    assert p.gate_ok and (p.left, p.right) == (e.true_start, e.end)


def test_a_placement_that_misses_the_call_is_refused(fx, g, truth_model):
    """The window is widened by either LTR's deficit, so it can reach a neighbouring
    copy. A placement that does not overlap the element's own call is not a placement
    of that element."""
    n = [x for x in members(fx) if x.family == FAMILY and x.len_left == 400
         and x.strand == "+"][0]
    ghost = dataclasses.replace(n, start=n.end + 60, l1=n.end + 99,
                                r0=n.end + 900, end=n.end + 1299)
    assert lp.core(ghost, truth_model.seq, g) is None
    assert lp.propose(ghost, [truth_model], g, lp.Params()) is None

"""Reconciliation: turning per-round calls into depth buckets.

This is the algorithmic heart of the pipeline. Round N cannot know whether the
element it just found sits inside something round N-1 found, so containment is
resolved once, over the pooled set, afterwards.
"""

from __future__ import annotations

import pytest

from ltrquest.reconcile import (
    IUPAC_DEPTH_SEQ,
    apply_depth_masking,
    build_direct_children,
    build_updated_nest_status,
    compute_chain_inward,
)


def rec(key: str, chrom: str, start: int, end: int) -> dict:
    """A reconciler record, cut down to the fields these functions read."""
    return {"key": key, "chrom": chrom, "s": start, "e": end}


class TestBuildDirectChildren:
    def test_a_lone_descendant_is_direct(self):
        assert build_direct_children({"outer": ["inner"]}) == {"outer": ["inner"]}

    def test_a_grandchild_is_not_a_direct_child(self):
        # outer contains both mid and inner; mid contains inner.
        all_in = {"outer": ["mid", "inner"], "mid": ["inner"]}
        children = build_direct_children(all_in)
        assert children["outer"] == ["mid"]
        assert children["mid"] == ["inner"]

    def test_siblings_are_both_direct(self):
        all_in = {"outer": ["a", "b"]}
        assert sorted(build_direct_children(all_in)["outer"]) == ["a", "b"]

    def test_no_containment_yields_nothing(self):
        assert build_direct_children({}) == {}

    def test_three_deep_chain_keeps_one_child_per_level(self):
        all_in = {
            "l3": ["l2", "l1", "l0"],
            "l2": ["l1", "l0"],
            "l1": ["l0"],
        }
        children = build_direct_children(all_in)
        assert children["l3"] == ["l2"]
        assert children["l2"] == ["l1"]
        assert children["l1"] == ["l0"]


class TestComputeChainInward:
    def test_an_element_with_nothing_inside_is_depth_zero(self):
        assert compute_chain_inward(["solo"], {}) == {"solo": 0}

    def test_one_layer_inside_is_depth_one(self):
        depths = compute_chain_inward(["outer", "inner"], {"outer": ["inner"]})
        assert depths == {"outer": 1, "inner": 0}

    def test_depth_counts_layers_not_elements(self):
        # Two elements side by side inside one host is still a single layer.
        depths = compute_chain_inward(["outer", "a", "b"], {"outer": ["a", "b"]})
        assert depths["outer"] == 1

    def test_depth_follows_the_longest_chain(self):
        children = {"outer": ["shallow", "mid"], "mid": ["deep"]}
        depths = compute_chain_inward(["outer", "shallow", "mid", "deep"], children)
        assert depths == {"outer": 2, "shallow": 0, "mid": 1, "deep": 0}

    def test_a_four_level_chain(self):
        children = {"l3": ["l2"], "l2": ["l1"], "l1": ["l0"]}
        depths = compute_chain_inward(["l3", "l2", "l1", "l0"], children)
        assert [depths[k] for k in ("l3", "l2", "l1", "l0")] == [3, 2, 1, 0]


class TestNestStatus:
    def test_unnested_elements_get_the_placeholder(self):
        assert build_updated_nest_status("solo", {}) == "."

    def test_a_host_records_its_tenant(self):
        status = build_updated_nest_status("outer", {"outer": ["inner"]})
        assert status == "nest-outer:inner"

    def test_a_tenant_records_its_host(self):
        status = build_updated_nest_status("inner", {"outer": ["inner"]})
        assert status == "nest-inner:outer"

    def test_an_element_can_be_both(self):
        all_in = {"outer": ["mid", "inner"], "mid": ["inner"]}
        status = build_updated_nest_status("mid", all_in)
        assert "nest-outer:inner" in status
        assert "nest-inner:outer" in status
        assert status.count(";") == 1

    def test_relations_are_deduplicated(self):
        status = build_updated_nest_status("outer", {"outer": ["inner", "inner"]})
        assert status == "nest-outer:inner"


class TestApplyDepthMasking:
    """Each nested element is masked with the IUPAC letter for its own depth."""

    def test_a_depth0_child_is_masked_with_n(self):
        outer = rec("chr1:1-100", "chr1", 1, 100)
        inner = rec("chr1:41-50", "chr1", 41, 50)
        out = apply_depth_masking(
            "A" * 100,
            outer,
            {outer["key"]: [inner["key"]]},
            {outer["key"]: outer, inner["key"]: inner},
            {outer["key"]: 1, inner["key"]: 0},
        )
        assert out[:40] == "A" * 40
        assert out[40:50] == "N" * 10
        assert out[50:] == "A" * 50

    def test_a_grandchild_overwrites_its_parents_mark(self):
        # depth2 host -> depth1 child (R) -> depth0 grandchild (N).
        outer = rec("chr1:1-100", "chr1", 1, 100)
        mid = rec("chr1:21-80", "chr1", 21, 80)
        deep = rec("chr1:41-60", "chr1", 41, 60)
        out = apply_depth_masking(
            "A" * 100,
            outer,
            {outer["key"]: [mid["key"]], mid["key"]: [deep["key"]]},
            {r["key"]: r for r in (outer, mid, deep)},
            {outer["key"]: 2, mid["key"]: 1, deep["key"]: 0},
        )
        assert out[:20] == "A" * 20        # host sequence
        assert out[20:40] == "R" * 20      # depth1 child
        assert out[40:60] == "N" * 20      # depth0 grandchild, painted over R
        assert out[60:80] == "R" * 20      # rest of the depth1 child
        assert out[80:] == "A" * 20

    def test_length_is_preserved(self):
        outer = rec("chr1:1-100", "chr1", 1, 100)
        inner = rec("chr1:41-50", "chr1", 41, 50)
        out = apply_depth_masking(
            "A" * 100, outer,
            {outer["key"]: [inner["key"]]},
            {outer["key"]: outer, inner["key"]: inner},
            {outer["key"]: 1, inner["key"]: 0},
        )
        assert len(out) == 100

    def test_a_child_on_another_contig_is_ignored(self):
        outer = rec("chr1:1-100", "chr1", 1, 100)
        stray = rec("chr2:41-50", "chr2", 41, 50)
        out = apply_depth_masking(
            "A" * 100, outer,
            {outer["key"]: [stray["key"]]},
            {outer["key"]: outer, stray["key"]: stray},
            {outer["key"]: 1, stray["key"]: 0},
        )
        assert out == "A" * 100

    def test_a_child_overhanging_the_host_is_clipped(self):
        outer = rec("chr1:1-100", "chr1", 1, 100)
        inner = rec("chr1:91-150", "chr1", 91, 150)
        out = apply_depth_masking(
            "A" * 100, outer,
            {outer["key"]: [inner["key"]]},
            {outer["key"]: outer, inner["key"]: inner},
            {outer["key"]: 1, inner["key"]: 0},
        )
        assert len(out) == 100
        assert out[90:] == "N" * 10

    def test_no_children_leaves_the_sequence_alone(self):
        outer = rec("chr1:1-100", "chr1", 1, 100)
        out = apply_depth_masking("ACGT" * 25, outer, {}, {outer["key"]: outer}, {})
        assert out == "ACGT" * 25

    @pytest.mark.parametrize("depth,char", list(enumerate(IUPAC_DEPTH_SEQ)))
    def test_every_depth_has_its_own_letter(self, depth, char):
        outer = rec("chr1:1-100", "chr1", 1, 100)
        inner = rec("chr1:41-50", "chr1", 41, 50)
        out = apply_depth_masking(
            "A" * 100, outer,
            {outer["key"]: [inner["key"]]},
            {outer["key"]: outer, inner["key"]: inner},
            {outer["key"]: depth + 1, inner["key"]: depth},
        )
        assert out[40:50] == char * 10


class TestIupacInvariants:
    """The depth letters are a contract between the driver and the reconciler."""

    def test_v_is_reserved_for_the_drivers_far_character(self):
        assert "V" not in IUPAC_DEPTH_SEQ

    def test_letters_are_unique(self):
        assert len(set(IUPAC_DEPTH_SEQ)) == len(IUPAC_DEPTH_SEQ)

    def test_letters_match_the_drivers_iupac_seq(self, driver_iupac_seq):
        # ltrquest.sh paints round N's hits with IUPAC_SEQ[N-1]; the reconciler
        # reads that letter back as a depth. If the two lists ever drift, nesting
        # is silently misattributed, so the coupling is asserted rather than
        # merely commented.
        assert list(IUPAC_DEPTH_SEQ) == driver_iupac_seq

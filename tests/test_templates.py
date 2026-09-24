"""Trusted templates for re-boundarying: which calls vouch for their own ends, and
which of them each target is placed with."""

from __future__ import annotations

import os

import pytest

from ltrquest import templates as tp
from ltrquest.ltr_model import Genomes, Member, rc
from reboundary_fixtures import build, members

TRUNCATED = ("del_left", "patch_right", "ins_left", "del_left_minus", "nested_del_left",
             "short3_left", "long3_right", "unstranded_del_left", "split_short")


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build(tmp_path_factory.mktemp("tpl"))


@pytest.fixture(scope="module")
def g(fx):
    return Genomes({fx.prefix: str(fx.genome)})


@pytest.fixture(scope="module")
def mem(fx):
    return members(fx)


def of_kind(fx, mem, kind):
    name = fx.kind(kind).name
    return next(m for m in mem if m.name == name)


def _m(**kw) -> Member:
    base = {"prefix": "p", "name": "x", "chrom": "chrT", "start": 1, "end": 2, "l1": 1,
            "r0": 2, "strand": "+", "orientation": "+", "family": "f", "depth": 0,
            "k2p": 0.01, "tsd": "ACGTA", "tsd_offset": "0,0", "nest_status": "."}
    base.update(kw)
    base["name"] = f"{base['chrom']}:{base['start']}-{base['end']}#LTR/x"
    return Member(**base)


# ---------------------------------------------------------------- who is trusted

def test_only_calls_with_an_exact_tsd_are_templates(fx, g, mem):
    kinds = sorted({next(e.kind for e in fx.elements if e.name == m.name)
                    for m in tp.trusted(mem, g)})
    assert kinds == ["host", "normal"]


def test_a_shifted_tsd_a_gap_or_a_nested_element_in_an_ltr_disqualifies_a_template(tmp_path):
    import random
    r = random.Random(3)
    seq = "".join(r.choice("ACGT") for _ in range(3000))
    seq = seq[:99] + "NNNNN" + seq[104:]                    # a gap at 100..104
    fa = tmp_path / "t.fa"
    fa.write_text(">chrT\n" + seq + "\n")
    g = Genomes({"p": str(fa)})
    gap = _m(start=90, l1=200, r0=800, end=910)
    nested = _m(start=1000, l1=1100, r0=1500, end=1600,
                nest_status="nest-outer:chrT:1050-1300")    # reaches into the left LTR
    clean = _m(start=1700, l1=1800, r0=2200, end=2300,
               nest_status="nest-outer:chrT:1850-2100")     # internal region only
    shifted = _m(start=2400, l1=2500, r0=2800, end=2900, tsd_offset="0,1")
    assert tp.trusted([gap, nested, clean, shifted], g) == [clean]


def test_targets_are_exactly_the_calls_without_a_tsd(fx, mem):
    kinds = sorted(next(e.kind for e in fx.elements if e.name == m.name)
                   for m in tp.targets(mem))
    assert kinds == sorted(["decayed", "conflict_left", "neighbour", "split_part",
                            *TRUNCATED])


def test_a_high_copy_family_is_capped_to_its_youngest_plus_a_seeded_draw():
    fam = [_m(start=100 * i, end=100 * i + 50, family="big", k2p=i / 1000) for i in range(1, 31)]
    small = [_m(start=5000 + 100 * i, end=5050 + 100 * i, family="small") for i in range(3)]
    loose = [_m(start=9000 + 100 * i, end=9050 + 100 * i, family=".") for i in range(40)]
    kept = tp.cap_per_family(fam + small + loose, cap=10)
    big = [m for m in kept if m.family == "big"]
    assert len(big) == 10
    assert [m.k2p for m in big[:5]] == [0.001, 0.002, 0.003, 0.004, 0.005]   # youngest half
    assert all(m.k2p > 0.005 for m in big[5:])                              # then the draw
    assert [m for m in kept if m.family == "small"] == small                 # under the cap
    assert len([m for m in kept if m.family == "."]) == 40                   # unlabelled: kept
    assert tp.cap_per_family(fam + small + loose, cap=10) == kept            # deterministic
    assert tp.cap_per_family(fam, cap=0) == fam                              # 0 = no cap


# ---------------------------------------------------------------- ranking (pure)

def test_each_template_counts_once_at_its_best_row_and_the_top_k_win():
    t = _m(start=10, end=20, tsd=".")
    a, b, c, d = (_m(start=100 * i, end=100 * i + 50) for i in range(1, 5))
    hits = [(t, a, 100.0, "+"), (t, a, 50.0, "-"), (t, c, 80.0, "+"), (t, b, 80.0, "+"),
            (t, d, 10.0, "+")]
    got = tp.rank(hits, k=3)[t.uid]
    assert [(m.start, rel) for m, rel in got] == [(100, "+"), (200, "+"), (300, "+")]


def test_a_template_whose_relation_contradicts_both_strands_is_dropped():
    t = _m(start=10, end=20, strand="+", tsd=".")
    wrong = _m(start=100, end=150, strand="-")             # opposite strand, same frame
    right = _m(start=200, end=250, strand="-")
    unknown = _m(start=300, end=350, strand=".")
    hits = [(t, wrong, 90.0, "+"), (t, right, 80.0, "-"), (t, unknown, 70.0, "+")]
    got = tp.rank(hits, k=5)[t.uid]
    assert [(m.start, rel) for m, rel in got] == [(200, "-"), (300, "+")]


def test_an_exact_tie_between_relations_is_broken_by_the_strands_else_plus():
    t = _m(start=10, end=20, strand="+", tsd=".")
    minus = _m(start=100, end=150, strand="-")
    unknown = _m(start=200, end=250, strand=".")
    hits = [(t, minus, 90.0, "+"), (t, minus, 90.0, "-"),
            (t, unknown, 80.0, "-"), (t, unknown, 80.0, "+")]
    got = tp.rank(hits, k=5)[t.uid]
    assert [(m.start, rel) for m, rel in got] == [(100, "-"), (200, "+")]


def test_excluded_templates_are_never_returned():
    t = _m(start=10, end=20, tsd=".")
    a, b = _m(start=100, end=150), _m(start=200, end=250)
    got = tp.rank([(t, a, 90.0, "+"), (t, b, 80.0, "+")], k=5, exclude=frozenset({a.uid}))
    assert [m.start for m, _ in got[t.uid]] == [200]


# ---------------------------------------------------------------- BLAST search

def test_every_truncated_call_gets_at_least_two_trusted_templates(fx, g, mem, blastn, tmp_path):
    trusted = tp.trusted(mem, g)
    got = tp.nearest(tp.targets(mem), trusted, g, k=5, workdir=str(tmp_path), blastn=blastn)
    uids = {m.uid for m in trusted}
    for kind in TRUNCATED:
        m = of_kind(fx, mem, kind)
        hits = got.get(m.uid, [])
        assert len(hits) >= 2, kind
        assert all(t.uid in uids and t.uid != m.uid for t, _ in hits), kind


def test_templates_of_both_strands_keep_their_own_relation(fx, g, mem, blastn, tmp_path):
    trusted = [m for m in tp.trusted(mem, g) if m.family == "syn_fam00001"]
    got = tp.nearest(tp.targets(mem), trusted, g, k=len(trusted), workdir=str(tmp_path),
                     blastn=blastn)
    for kind, frame in (("del_left", "+"), ("del_left_minus", "-"),
                        ("unstranded_del_left", "+")):
        hits = got[of_kind(fx, mem, kind).uid]
        assert {t.strand for t, _ in hits} == {"+", "-"}, kind
        for t, rel in hits:
            assert rel == ("+" if t.strand == frame else "-"), (kind, t.name, rel)


def test_makeblastdb_on_path_is_used_when_none_sits_beside_blastn(fx, g, mem, blastn, tmp_path,
                                                                   monkeypatch):
    """The start-up check accepts a makeblastdb found only on PATH; the search must then
    run that one, not a missing one beside blastn."""
    wrapper = tmp_path / "bin" / "blastn"
    wrapper.parent.mkdir()
    wrapper.write_text(f'#!/bin/sh\nexec "{blastn}" "$@"\n')
    wrapper.chmod(0o755)
    monkeypatch.setenv("PATH", os.path.dirname(blastn) + os.pathsep + os.environ["PATH"])
    got = tp.nearest(tp.targets(mem), tp.trusted(mem, g), g, k=5, workdir=str(tmp_path / "w"),
                     blastn=str(wrapper))
    assert got


def test_a_query_whose_blast_rows_are_not_contiguous_is_refused(fx, g, mem, blastn, tmp_path):
    """Each target is ranked as soon as BLAST moves past its rows, so only its top k is
    held; that needs every target's rows together, as blastn writes them."""
    stub = tmp_path / "bin" / "blastn"
    stub.parent.mkdir()
    stub.write_text("#!/usr/bin/env python3\nimport sys\n"
                    "out = sys.argv[sys.argv.index('-out') + 1]\n"
                    "open(out, 'w').write('q0L\\tt0\\t90.0\\tplus\\n'\n"
                    "                     'q1L\\tt1\\t80.0\\tplus\\n'\n"
                    "                     'q0R\\tt1\\t70.0\\tplus\\n')\n")
    stub.chmod(0o755)
    with pytest.raises(RuntimeError, match="not grouped"):
        tp.nearest(tp.targets(mem)[:2], tp.trusted(mem, g), g, k=5,
                   workdir=str(tmp_path / "w"), blastn=str(stub))


def test_the_search_is_deterministic(fx, g, mem, blastn, tmp_path):
    args = (tp.targets(mem), tp.trusted(mem, g), g)
    one = tp.nearest(*args, k=5, workdir=str(tmp_path / "a"), blastn=blastn)
    two = tp.nearest(*args, k=5, workdir=str(tmp_path / "b"), blastn=blastn)
    assert {u: [(t.uid, r) for t, r in v] for u, v in one.items()} == \
        {u: [(t.uid, r) for t, r in v] for u, v in two.items()}


def test_two_genomes_with_identical_keys_do_not_collide(tmp_path, blastn):
    fa, fb = build(tmp_path / "a", prefix="A_LTRs"), build(tmp_path / "b", prefix="B_LTRs")
    g = Genomes({fa.prefix: str(fa.genome), fb.prefix: str(fb.genome)})
    mem = members(fa) + members(fb)
    got = tp.nearest(tp.targets(mem), tp.trusted(mem, g), g, k=24, workdir=str(tmp_path / "w"),
                     blastn=blastn)
    target = of_kind(fa, members(fa), "del_left")
    prefixes = {t.prefix for t, _ in got[target.uid]}
    assert prefixes == {"A_LTRs", "B_LTRs"}
    assert len({t.uid for t, _ in got[target.uid]}) == len(got[target.uid])


# ---------------------------------------------------------------- side models

def test_each_side_is_placed_with_a_template_ltr_whose_outer_end_is_verified(fx, g, mem):
    t = next(m for m in tp.trusted(mem, g) if m.strand == "+" and m.family == "syn_fam00001")
    left = g.fetch(t.prefix, t.chrom, t.start, t.l1)[0]
    right = g.fetch(t.prefix, t.chrom, t.r0, t.end)[0]
    assert tp.models_for(t, "+", g) == (left, right)
    assert tp.models_for(t, "-", g) == (rc(right), rc(left))
    for rel in "+-":
        lm, rm = tp.models_for(t, rel, g)
        assert lm.startswith("TG") and rm.endswith("CA"), rel     # the TSD-verified ends

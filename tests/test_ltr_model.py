"""Genome access and the end-to-end alignment re-boundarying places templates with."""

from __future__ import annotations

import random

import pytest

from ltrquest import ltr_model as lm
from reboundary_fixtures import FAMILY, build, members


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build(tmp_path_factory.mktemp("syn_model"))


@pytest.fixture(scope="module")
def g(fx):
    return lm.Genomes({fx.prefix: str(fx.genome)})


@pytest.fixture(scope="module")
def fam(fx):
    return [m for m in members(fx) if m.family == FAMILY]


def test_genomes_fetch_clips_to_the_contig(fx, g):
    seq, first = g.fetch(fx.prefix, "chrS", -5, 10)
    assert first == 1 and len(seq) == 10
    n = g.length(fx.prefix, "chrS")
    assert g.fetch(fx.prefix, "chrS", n - 2, n + 50)[0] == g.fetch(
        fx.prefix, "chrS", n - 2, n)[0]


def test_genomes_survive_pickling(fx, g):
    import pickle
    g2 = pickle.loads(pickle.dumps(g))
    assert g2.fetch(fx.prefix, "chrS", 1, 20) == g.fetch(fx.prefix, "chrS",
                                                          1, 20)


def test_member_geometry(fam):
    m = [x for x in fam if x.name.startswith("chrS") and x.len_left == 400][0]
    assert m.len_right == 400 and m.key == m.name.split("#")[0] and m.uid.endswith(
        m.key)


def test_glocal_places_a_model_exactly():
    rng = random.Random(4)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    model = rand(300)
    window = rand(120) + model + rand(80)
    hit = lm.glocal(model, window)
    assert (hit.start, hit.end) == (120, 419)
    assert hit.id_first == hit.id_last == hit.id_all == 1.0

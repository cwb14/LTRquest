"""Family LTR models: termini from 5'/3' agreement, never from a TSD or a
motif."""

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


def test_modal_length_picks_the_dense_bin():
    assert abs(lm.modal_length([400] * 20 + [290, 330, 360, 1900]) - 400) <= 5


def test_select_references_takes_the_modal_class_and_honours_exclusion(fam):
    refs, modal = lm.select_references(fam, "modal", FAMILY)
    assert abs(modal - 400) <= 5
    assert all(abs(r.len_left - modal) <= 0.15 * modal for r in refs)
    assert len(refs) >= lm.MIN_REFS
    keys = {r.key for r in refs}
    refs2, _ = lm.select_references(fam, "modal", FAMILY,
                                    exclude=frozenset(list(keys)[:2]))
    assert not (set(list(keys)[:2]) & {r.key for r in refs2})


def test_select_references_by_tsd_uses_only_tsd_bearing_copies(fam):
    refs, _ = lm.select_references(fam, "tsd", FAMILY)
    assert refs and all(r.has_tsd for r in refs)


def test_consensus_from_alignment_finds_the_ltr_between_disagreeing_flanks():
    rng = random.Random(3)
    ltr = "TG" + "".join(rng.choice("ACGT") for _ in range(196)) + "CA"
    rows, sides = [], []
    for _ in range(12):   # pre-aligned rows, no MAFFT: 5' rows = flank A +
        rows.append("A" * 60 + ltr + "G" * 60)  # LTR + internal G, 3' rows =
        sides.append(5)                          # internal C + LTR + flank T, so
        rows.append("C" * 60 + ltr + "T" * 60)  # the two sides disagree
        sides.append(3)                          # everywhere but the LTR
    assert lm.consensus_from_alignment(rows, sides) == ltr


def test_glocal_places_a_model_exactly():
    rng = random.Random(4)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    model = rand(300)
    window = rand(120) + model + rand(80)
    hit = lm.glocal(model, window)
    assert (hit.start, hit.end) == (120, 419)
    assert hit.id_first == hit.id_last == hit.id_all == 1.0


def test_consensus_model_recovers_the_planted_ltr(fx, g, fam, mafft):
    refs, modal = lm.select_references(fam, "modal", FAMILY)
    model = lm.consensus_model(FAMILY, refs, g, modal, mafft)
    assert model is not None and abs(len(model.seq) - 400) <= 2
    assert model.seq[:20] == fx.ltr[:20] and model.seq[-20:] == fx.ltr[-20:]
    assert lm.glocal(model.seq, fx.ltr).id_all >= 0.99


def test_polish_ends_restores_clipped_termini(fx, g, fam):
    refs, _ = lm.select_references(fam, "modal", FAMILY)
    seq, (ds, de) = lm.polish_ends(fx.ltr[3:-3], refs, g)
    assert (ds, de) == (3, 3) and seq[:10] == fx.ltr[:10] and seq[-10:] == fx.ltr[-10:]


def test_polish_ends_trims_overhangs(fx, g, fam):
    refs, _ = lm.select_references(fam, "modal", FAMILY)
    seq, (ds, de) = lm.polish_ends("ACG" + fx.ltr + "TTA", refs, g)
    assert (ds, de) == (-3, -3) and seq[:10] == fx.ltr[:10] and seq[-10:] == fx.ltr[-10:]


def test_a_model_does_not_carry_its_kmer_index_through_a_pickle():
    """`nearest` makes up to 80 models per family; every one used to ship its k-mer
    set back to the parent and into --cache."""
    import pickle
    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(1600))
    mo = lm.Model("f:consensus", "f", seq, 400.0, 12)
    assert len(mo.kmers) > 1000                     # built on demand
    blob = pickle.dumps(mo)
    assert len(blob) < 2 * len(seq)
    back = pickle.loads(blob)
    assert back == mo and back.kmers == mo.kmers


def test_a_model_records_the_copies_it_was_built_from(fx, g, fam, mafft):
    refs, modal = lm.select_references(fam, "modal", FAMILY)
    model = lm.consensus_model(FAMILY, refs, g, modal, mafft)
    assert model.refs == {r.uid for r in refs}


def test_ratio_ok():
    m = lm.Model("f:consensus", "f", "A" * 459, 400.0, 20)
    assert lm.ratio_ok(m, 1.15) and not lm.ratio_ok(m, 1.10)


def test_binomial_tail():
    assert lm.binom_sf(0, 5, 0.3) == 1.0
    assert abs(lm.binom_sf(5, 5, 0.5) - 1 / 32) < 1e-12


def test_binomial_tail_survives_a_family_with_thousands_of_candidates():
    # The real failure: math.comb(5000, i) is an int with thousands of digits.
    assert 0.0 < lm.binom_sf(120, 5000, 0.02) < 1.0
    assert lm.binom_sf(3000, 5000, 0.02) >= 0.0          # deep tail: underflows to 0, never raises
    assert lm.binom_sf(120, 5000, 0.02) > lm.binom_sf(160, 5000, 0.02)


def test_binomial_tail_still_matches_the_exact_sum_for_small_n():
    import math as _math
    for k, n, p in ((3, 20, 0.1), (5, 12, 0.5), (1, 8, 0.01)):
        exact = sum(_math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))
        assert abs(lm.binom_sf(k, n, p) - exact) < 1e-12


def test_tsd_enrichment_has_an_explicit_untested_branch():
    assert lm.tsd_enrichment(3, 3, 0.01, min_n=5) == "untested"
    assert lm.tsd_enrichment(6, 8, 0.01, min_n=5) == "pass"
    assert lm.tsd_enrichment(0, 8, 0.01, min_n=5) == "fail"


def test_subfamily_models_on_a_one_type_family_give_the_family_ltr(fx, g, fam, mafft):
    refs, modal = lm.select_references(fam, "modal", FAMILY, n_young=150, n_random=150)
    models = lm.subfamily_models(FAMILY, refs, g, modal, 0.5, mafft)
    assert models and all(abs(len(mo.seq) - 400) <= 2 for mo in models)
    assert models[0].seq[:20] == fx.ltr[:20]

"""The driver: family QC, Kmer2LTR arbitration, and what a real rewrite leaves behind."""

from __future__ import annotations

import json
import os
import shutil

import pytest

from ltrquest import reboundary as rb
from ltrquest.detect import revcomp as revcomp_record
from ltrquest.kmer2ltr import COLUMNS
from ltrquest.ltr_model import Member, Model
from ltrquest.ltr_place import Proposal
from ltrquest.reboundary_io import fetch_records

from reboundary_fixtures import build

TOOLS = os.environ.get("LTRQUEST_TOOLS_DIR", ".")


def fake_member(**kw):
    base = dict(prefix="p", name="c:100-900#LTR/Gypsy/X", chrom="c", start=100, end=900,
                l1=199, r0=801, strand="+", orientation="+", family="f", depth=0, k2p=0.01,
                tsd=".")
    base.update(kw)
    return Member(**base)


def fake_proposal(m, left, right, gate=True, id_left=0.9, id_right=0.9):
    return Proposal(m, "f:consensus", "+", left, right, "core", "core", id_left, id_right, 0.9,
                    left, right, m.l1, m.r0, gate)


def test_rebase_restates_coordinates_against_the_cut_record():
    f = ["NA"] * len(COLUMNS)
    i = {c: k for k, c in enumerate(COLUMNS)}
    f[i["ltr5_start"]], f[i["ltr5_end"]] = "11", "210"
    f[i["ltr3_start"]], f[i["ltr3_end"]] = "811", "1010"
    f[i["flank5_len"]], f[i["flank3_len"]] = "10", "5"
    out = rb.rebase(f, 10, "c:110-1109#x")
    assert out[i["seq_id"]] == "c:110-1109#x"
    assert [out[i[c]] for c in ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end", "seq_len")] \
        == ["1", "200", "801", "1000", "1000"]
    assert out[i["flank5_len"]] == out[i["flank3_len"]] == "0"


def test_length_gate_matches_detection():
    assert rb.length_ok(400, 390, 405, 1, 3000)
    assert not rb.length_ok(400, 150, 405, 1, 3000)      # ratio < 0.65
    assert not rb.length_ok(90, 90, 95, 1, 3000)         # LTR < 100
    assert not rb.length_ok(400, 400, 400, 1, 200)       # element < 300


def test_credit_for_constant_and_model():
    m = fake_member()
    p = fake_proposal(m, 50, 900, id_left=0.9)
    assert rb.credit_for(p, rb.Settings(credit="200")) == 200.0
    s = rb.Settings(credit="model")
    assert rb.credit_for(p, s) == pytest.approx(2 * s.place.anchor_len * 0.9)


def test_qc_passes_families_whose_proposals_gain_tsds():
    fam = rb.FamilyResult("f", 12, [Model("f:consensus", "f", "A" * 400, 400.0, 12)], 400.0, "ok")
    for k in range(6):
        m = fake_member(name=f"c:{1000 * k + 100}-{1000 * k + 900}#x", start=1000 * k + 100,
                        end=1000 * k + 900, l1=1000 * k + 199, r0=1000 * k + 801)
        fam.proposals.append(fake_proposal(m, m.start - 50, m.end))
        fam.tsd_new[m.uid], fam.tsd_null[m.uid] = "ACGTA", "."
    status, p0 = rb.qc([fam], rb.Settings())
    assert status["f"] == "pass" and p0 == 0.0
    assert rb.family_accepts("pass", rb.Settings())
    assert rb.family_accepts("untested", rb.Settings(untested="accept"))
    assert not rb.family_accepts("untested", rb.Settings(untested="skip"))
    assert not rb.family_accepts("fail", rb.Settings())


def test_mutation_rate_comes_from_the_run_record(tmp_path):
    (tmp_path / "p.detect.json").write_text(json.dumps({"settings": {"mutation_rate": "1.3e-8"}}))
    assert rb.resolve_mutation_rate(str(tmp_path), ["p"], None) == pytest.approx(1.3e-8)
    assert rb.resolve_mutation_rate(str(tmp_path), ["p"], 5e-9) == 5e-9
    assert rb.resolve_mutation_rate(str(tmp_path / "none"), ["p"], None) == 3e-8


def test_a_non_numeric_credit_is_refused():
    with pytest.raises(SystemExit):
        rb.main(["--prefix", "p", "--genome", "g.fa", "--credit", "lots"])


def sidecar(fx):
    lines = (fx.indir / f"{fx.prefix}_reboundary.tsv").read_text().splitlines()
    head = lines[0].lstrip("#").split("\t")
    rows = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    return {r["old_seq_id"].split("#")[0]: r for r in rows}


@pytest.fixture(scope="module")
def ran(tmp_path_factory, k2l_api, mafft):
    fx = build(tmp_path_factory.mktemp("syn_run"))
    raw = {p.name: p.read_bytes() for p in fx.indir.glob(f"{fx.prefix}_depth*_ltr.*")
           if "_clean_" not in p.name}
    counts = rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)],
                    rb.Settings(threads=2, mafft=mafft), TOOLS)
    return fx, counts, raw


@pytest.mark.parametrize("kind", ["del_left", "patch_right", "del_left_minus", "nested_del_left"])
def test_truncated_calls_are_extended_to_their_true_ends(ran, kind):
    fx, _, _ = ran
    e = fx.kind(kind)
    row = sidecar(fx)[e.key]
    assert row["decision"] == "extended", row
    assert row["new_seq_id"].split("#")[0] == f"chrS:{e.true_start}-{e.true_end}"
    assert row["tsd_new"] not in (".", "NA")


def test_a_clash_with_a_neighbour_is_rejected(ran):
    fx, _, _ = ran
    assert sidecar(fx)[fx.kind("conflict_left").key]["reason"] == "overlaps_element"


def test_an_insertion_beyond_ltrquests_pair_rules_is_rejected_with_a_reason(ran):
    fx, _, _ = ran
    row = sidecar(fx)[fx.kind("ins_left").key]
    assert row["decision"] == "rejected" and row["end_source5"] == "anchor"
    assert row["reason"] in ("kmer2ltr_reverted", "length_filter", "kmer2ltr_not_pass")


def test_intact_copies_are_not_candidates(ran):
    fx, _, _ = ran
    side = sidecar(fx)
    assert not [e for e in fx.elements if e.kind in ("normal", "decayed") and e.key in side]


def test_rewritten_rows_come_from_kmer2ltr_and_keep_the_annotation(ran):
    fx, _, _ = ran
    e = fx.kind("del_left")
    new = f"chrS:{e.true_start}-{e.true_end}#LTR/Gypsy/Synth"
    lines = (fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv").read_text().splitlines()
    head = lines[0].lstrip("#").split("\t")
    row = dict(zip(head, [l for l in lines if l.startswith(new + "\t")][0].split("\t")))
    assert row["ltr5_start"] == "1" and row["seq_len"] == str(e.true_end - e.true_start + 1)
    assert row["flank5_len"] == row["flank3_len"] == "0" and row["status"] == "pass"
    assert row["orientation"] == "+" and row["strand"] == "+" and row["family"] == e.family
    assert row["tsd"] not in (".", "NA")


def test_a_minus_record_is_stored_reverse_complemented(ran):
    fx, _, _ = ran
    e = fx.kind("del_left_minus")
    new = f"chrS:{e.true_start}-{e.true_end}#LTR/Gypsy/Synth"
    genome = "".join(fx.genome.read_text().splitlines()[1:])
    rec = fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], {new})[new]
    assert rec == revcomp_record(genome[e.true_start - 1:e.true_end])


def test_the_host_masks_the_recovered_bases_and_points_at_the_new_name(ran):
    fx, _, _ = ran
    e, host = fx.kind("nested_del_left"), fx.kind("host")
    rec = fetch_records([str(fx.indir / f"{fx.prefix}_depth1_clean_ltr.fa")],
                        {host.name})[host.name]
    assert set(rec[e.true_start - host.start:e.start - host.start]) == {"N"}
    t1 = (fx.indir / f"{fx.prefix}_depth1_clean_ltr.tsv").read_text()
    assert f"nest-outer:chrS:{e.true_start}-{e.true_end}" in t1


def test_raw_tables_are_untouched(ran):
    fx, _, raw = ran
    for name, data in raw.items():
        assert (fx.indir / name).read_bytes() == data, name


def test_genomes_from_different_family_namespaces_are_refused(tmp_path, k2l_api, mafft):
    fx = build(tmp_path / "syn_ns")
    path = fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv"
    lines = path.read_text().splitlines()
    lines[1] = lines[1].replace("\tsyn_fam00001\t", "\tother_fam00001\t")
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(SystemExit, match="namespaces"):
        rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=1, mafft=mafft),
               TOOLS)


def test_a_prefix_with_no_clean_tables_is_skipped_not_fatal(tmp_path, k2l_api, mafft):
    fx = build(tmp_path / "syn_partial")
    counts = rb.run(str(fx.indir), [fx.prefix, "missing_p"], [str(fx.genome), str(fx.genome)],
                    rb.Settings(threads=2, mafft=mafft), TOOLS)
    assert counts.get("extended", 0) == 4
    e = fx.kind("del_left")
    assert sidecar(fx)[e.key]["decision"] == "extended"
    assert not (fx.indir / "missing_p_reboundary.tsv").exists()


def test_all_prefixes_empty_raises_systemexit(tmp_path, k2l_api, mafft):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    with pytest.raises(SystemExit, match="no _clean_ depth tables"):
        rb.run(str(tmp_path), ["missing"], [str(genome)], rb.Settings(threads=1, mafft=mafft),
               TOOLS)


def test_a_second_run_extends_nothing_and_the_cli_works(tmp_path, k2l_api, mafft):
    fx = build(tmp_path / "syn_twice")
    args = ["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
            "-t", "2", "--mafft", mafft, "--tools-dir", TOOLS]
    assert rb.main(args) == 0
    first = sidecar(fx)
    assert sum(r["decision"] == "extended" for r in first.values()) == 4
    counts = rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)],
                    rb.Settings(threads=2, mafft=mafft), TOOLS)
    assert counts.get("extended", 0) == 0


def test_posthoc_backs_up_once_reruns_from_originals_and_restores(tmp_path, k2l_api, mafft):
    fx = build(tmp_path / "syn_posthoc")
    originals = {p.name: p.read_bytes() for p in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.*")}
    (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").write_text("##gff-version 3\n#old\n")
    (fx.indir / f"{fx.prefix}_plots").mkdir()
    (fx.indir / f"{fx.prefix}_plots" / "x.pdf").write_text("old plot")
    base = ["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
            "-t", "2", "--mafft", mafft, "--tools-dir", TOOLS, "--posthoc", "--no-plots"]

    assert rb.main(base) == 0
    backup = fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}"
    assert {p.name for p in backup.iterdir()} >= set(originals) | {
        f"{fx.prefix}_all_depth_LTR_cleaned.gff3", f"{fx.prefix}_plots"}
    gff = (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").read_text()
    assert "boundary_source=family_model" in gff and "#old" not in gff
    assert not (fx.indir / f"{fx.prefix}_plots").exists()

    assert rb.main(base + ["--no-anchor"]) == 0      # starts again from the backup
    e = fx.kind("del_left")
    assert sidecar(fx)[e.key]["decision"] == "extended"

    assert rb.main(["--indir", str(fx.indir), "--prefix", fx.prefix, "--restore"]) == 0
    assert not backup.exists() and not (fx.indir / f"{fx.prefix}_reboundary.tsv").exists()
    for name, data in originals.items():
        assert (fx.indir / name).read_bytes() == data, name
    assert "#old" in (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").read_text()
    assert (fx.indir / f"{fx.prefix}_plots" / "x.pdf").read_text() == "old plot"


def test_restore_without_a_backup_says_so(tmp_path):
    with pytest.raises(SystemExit, match="nothing to restore"):
        rb.main(["--indir", str(tmp_path), "--prefix", "p", "--restore"])


def test_an_interrupted_backup_stops_the_run(tmp_path):
    (tmp_path / f"p{rb.BACKUP_SUFFIX}.partial").mkdir()
    with pytest.raises(SystemExit, match="interrupted"):
        rb.prepare_posthoc(str(tmp_path), "p")


def clean_bytes(fx):
    return {p.name: p.read_bytes() for p in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.*")}


def test_posthoc_with_an_empty_backup_keeps_the_live_files(tmp_path):
    fx = build(tmp_path / "syn_empty_backup")
    before = clean_bytes(fx)
    (fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}").mkdir()
    with pytest.raises(SystemExit, match=rb.BACKUP_SUFFIX):
        rb.prepare_posthoc(str(fx.indir), fx.prefix)
    assert clean_bytes(fx) == before


def test_posthoc_with_a_partial_backup_refuses_before_removing_anything(tmp_path):
    fx = build(tmp_path / "syn_partial_backup")
    before = clean_bytes(fx)
    bdir = fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}"
    bdir.mkdir()
    for p in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.*"):
        shutil.copy2(p, bdir / p.name)
    (bdir / f"{fx.prefix}_depth1_clean_ltr.fa").unlink()      # a table with no FASTA beside it
    with pytest.raises(SystemExit, match="depth1_clean_ltr"):
        rb.prepare_posthoc(str(fx.indir), fx.prefix)
    assert clean_bytes(fx) == before


def test_an_interrupted_restore_can_be_re_run(tmp_path, monkeypatch):
    fx = build(tmp_path / "syn_restore")
    originals = clean_bytes(fx)
    (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").write_text("##gff-version 3\n#old\n")
    (fx.indir / f"{fx.prefix}_plots").mkdir()
    (fx.indir / f"{fx.prefix}_plots" / "x.pdf").write_text("old plot")
    rb.prepare_posthoc(str(fx.indir), fx.prefix)
    bdir = fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}"
    held = {p.name for p in bdir.iterdir()}
    for name in originals:                                   # as a finished run leaves them
        (fx.indir / name).write_text("rewritten\n")
    (fx.indir / f"{fx.prefix}{rb.SIDECAR_SUFFIX}").write_text("#rewritten\n")

    real, seen = shutil.copy2, []

    def flaky(src, dst, **kw):
        seen.append(src)
        if len(seen) == 2:
            raise KeyboardInterrupt                          # killed part way through
        return real(src, dst, **kw)

    monkeypatch.setattr(rb.shutil, "copy2", flaky)
    with pytest.raises(KeyboardInterrupt):
        rb.restore(str(fx.indir), fx.prefix)
    monkeypatch.undo()
    assert {p.name for p in bdir.iterdir()} == held           # the backup is still whole

    rb.restore(str(fx.indir), fx.prefix)                      # the re-run finishes the job
    assert not bdir.exists()
    for name, data in originals.items():
        assert (fx.indir / name).read_bytes() == data, name
    assert "#old" in (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").read_text()
    assert (fx.indir / f"{fx.prefix}_plots" / "x.pdf").read_text() == "old plot"
    assert not (fx.indir / f"{fx.prefix}{rb.SIDECAR_SUFFIX}").exists()


def test_build_models_rejects_an_unknown_method():
    from ltrquest.ltr_model import Genomes
    m = fake_member()
    with pytest.raises(ValueError, match="unknown method"):
        rb.build_models("f", [m], Genomes({}), rb.Settings(method="bogus"))


@pytest.mark.parametrize("method", ["consensus", "subfamily", "nearest"])
def test_every_method_recovers_the_deletion(tmp_path, mafft, method):
    from ltrquest.ltr_model import Genomes
    from ltrquest.ltr_place import propose
    from reboundary_fixtures import FAMILY, members
    fx = build(tmp_path / "syn_methods")
    g = Genomes({fx.prefix: str(fx.genome)})
    fam = [m for m in members(fx) if m.family == FAMILY]
    s = rb.Settings(method=method, mafft=mafft)
    models, _, status = rb.build_models(FAMILY, fam, g, s)
    assert status == "ok" and models
    e = fx.kind("del_left")
    m = [x for x in fam if x.name == e.name][0]
    p = propose(m, models, g, rb.place_params(s))
    assert p.gate_ok and p.left == e.true_start


@pytest.mark.parametrize("method", ["consensus", "subfamily", "nearest"])
def test_no_element_is_re_boundaried_by_a_model_it_helped_build(tmp_path, mafft, method):
    """Leave-one-out in the production path: `build_models` is called the way the
    driver calls it, with no exclusion, and no element may still be placed by a model
    its own sequence went into."""
    from ltrquest.ltr_model import Genomes
    from ltrquest.ltr_place import candidate_models, propose
    from reboundary_fixtures import FAMILY, members
    fx = build(tmp_path / "syn_loo")
    g = Genomes({fx.prefix: str(fx.genome)})
    fam = [m for m in members(fx) if m.family == FAMILY]
    s = rb.Settings(method=method, mafft=mafft)
    models, _, status = rb.build_models(FAMILY, fam, g, s)
    assert status == "ok"
    built = {u for mo in models for u in mo.refs}
    assert len(built) >= 10                       # every model says who built it
    for m in (x for x in fam if x.uid in built):
        for mo in candidate_models(m, models, g, len(models)):
            assert m.uid not in mo.refs
        p = propose(m, models, g, rb.place_params(s))
        assert p is None or m.key not in p.model_id


def test_a_family_below_the_reference_floor_says_so_rather_than_failing(tmp_path):
    from ltrquest.ltr_model import Genomes
    from reboundary_fixtures import FAMILY, members
    fx = build(tmp_path / "syn_floor")
    g = Genomes({fx.prefix: str(fx.genome)})
    fam = [m for m in members(fx) if m.family == FAMILY]
    models, _, status = rb.build_models(FAMILY, fam, g, rb.Settings(), frozenset(
        m.key for m in fam[:-5]))
    assert (models, status) == ([], "too_few_references")


def test_nearest_uses_median_combining():
    assert rb.place_params(rb.Settings(method="nearest")).combine == "median"
    assert rb.place_params(rb.Settings(method="consensus")).combine == "best"

"""The driver: templates, placement, merges, Kmer2LTR re-measurement, and what a real
rewrite leaves behind."""

from __future__ import annotations

import json
import os
import shutil

import pytest

from ltrquest import reboundary as rb
from ltrquest.detect import revcomp as revcomp_record
from ltrquest.kmer2ltr import COLUMNS
from ltrquest.reboundary_io import NEW_SUFFIX, OLD_SUFFIX, fetch_records
from reboundary_fixtures import FAMILY, Element, _row, build

TOOLS = os.environ.get("LTRQUEST_TOOLS_DIR", ".")
MOVED_KINDS = ("del_left", "patch_right", "del_left_minus", "nested_del_left", "ins_left",
               "short3_left", "long3_right", "unstranded_del_left", "split_short")


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


def test_the_length_gate_keeps_detections_floors_but_not_its_ratio():
    assert rb.length_ok(400, 390, 405, 1, 3000)
    assert rb.length_ok(1900, 360, 1905, 1, 3000)        # a partner that lost part of its LTR
    assert not rb.length_ok(90, 90, 95, 1, 3000)         # LTR < 100
    assert not rb.length_ok(400, 400, 80, 1, 3000)       # alignment < 90
    assert not rb.length_ok(400, 400, 400, 1, 200)       # element < 300


@pytest.mark.parametrize("cigar, cut5, cut3, want", [
    ("100=", 3, 0, (0, 99, 503, 599)),                  # 3 paired bases leave the partner too
    ("2I98=", 3, 0, (0, 99, 501, 599)),                 # 2 of them were unpaired
    ("3=2D95=", 3, 0, (0, 99, 503, 599)),               # the partner's run after the cut stays
    ("100=", 0, 4, (0, 95, 500, 599)),                  # a right trim mirrors it
    ("97=3D", 0, 4, (0, 98, 500, 599)),                 # 3 of the 4 had no partner
    ("100=", 0, 0, (0, 99, 500, 599)),
    (".", 3, 0, (0, 99, 500, 599)),                     # no alignment: nothing to follow
])
def test_trim_spans_retract_the_partner_by_what_the_cut_bases_paired_with(cigar, cut5, cut3, want):
    assert rb.trim_spans(cigar, (0, 99, 500, 599), cut5, cut3) == want


def test_mutation_rate_comes_from_the_run_record(tmp_path):
    (tmp_path / "p.detect.json").write_text(json.dumps({"settings": {"mutation_rate": "1.3e-8"}}))
    assert rb.resolve_mutation_rate(str(tmp_path), ["p"], None) == pytest.approx(1.3e-8)
    assert rb.resolve_mutation_rate(str(tmp_path), ["p"], 5e-9) == 5e-9
    assert rb.resolve_mutation_rate(str(tmp_path / "none"), ["p"], None) == 3e-8


def test_the_genome_is_indexed_before_any_worker_starts(tmp_path, monkeypatch, blastn):
    # Workers open the genome lazily, all at once, and pyfaidx locks per process: an
    # unindexed genome had its .fai built by every worker, and one read it while
    # another had truncated it (KeyError on a real chromosome).
    fx = build(tmp_path / "syn_index")
    fai = f"{fx.genome}.fai"
    assert not os.path.exists(fai)
    seen = {}

    class Started(Exception):
        pass

    def pool(*args, **kwargs):
        seen["fai"] = os.path.isfile(fai) and os.path.getsize(fai) > 0
        raise Started

    monkeypatch.setattr(rb.k2l, "api", lambda *args, **kwargs: None)
    monkeypatch.setattr(rb, "Pool", pool)
    with pytest.raises(Started):
        rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(), TOOLS)
    assert seen == {"fai": True}


@pytest.mark.parametrize("flag", [
    ["--min-ext", "0"], ["--templates", "0"], ["--min-support", "0"],
    ["--templates", "3", "--min-support", "4"], ["--max-trim", "-1"], ["--anchor-len", "0"],
    ["--max-indel", "-5"], ["--min-identity", "0"], ["--min-identity", "1.5"], ["-t", "0"],
    ["--mutation-rate", "0"],
])
def test_out_of_range_settings_are_refused_before_any_work(tmp_path, flag, capsys):
    """A setting outside its range would run and quietly write nonsense (--min-ext 0
    makes every confirmed call a rejected row with no reason; --min-support above
    --templates can never move anything): the CLI refuses it up front."""
    with pytest.raises(SystemExit) as e:
        rb.main(["--indir", str(tmp_path), "--prefix", "p", "--genome", "g.fa"] + flag)
    assert e.value.code == 2
    assert flag[-2] in capsys.readouterr().err


@pytest.mark.parametrize("flag", [["--method", "nearest"], ["--credit", "5000"],
                                  ["--mafft", "mafft"], ["--min-copies", "10"]])
def test_the_family_model_flags_are_gone(flag):
    with pytest.raises(SystemExit):
        rb.main(["--prefix", "p", "--genome", "g.fa"] + flag)


def sidecar(fx):
    lines = (fx.indir / f"{fx.prefix}_reboundary.tsv").read_text().splitlines()
    head = lines[0].lstrip("#").split("\t")
    rows = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    return {r["old_seq_id"].split("#")[0]: r for r in rows}


def table_row(fx, depth, name):
    lines = (fx.indir / f"{fx.prefix}_depth{depth}_clean_ltr.tsv").read_text().splitlines()
    head = lines[0].lstrip("#").split("\t")
    hits = [l for l in lines if l.startswith(name + "\t")]
    return dict(zip(head, hits[0].split("\t"))) if hits else None


@pytest.fixture(scope="module")
def ran(tmp_path_factory, k2l_api, blastn):
    fx = build(tmp_path_factory.mktemp("syn_run"))
    raw = {p.name: p.read_bytes() for p in fx.indir.glob(f"{fx.prefix}_depth*_ltr.*")
           if "_clean_" not in p.name}
    before = {p.name: p.read_text() for p in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.fa")}
    counts = rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=2), TOOLS)
    return fx, counts, raw, before


@pytest.mark.parametrize("kind", MOVED_KINDS)
def test_every_mis_called_kind_is_moved_to_its_true_ends(ran, kind):
    fx, _, _, _ = ran
    e = fx.kind(kind)
    row = sidecar(fx)[e.key]
    assert row["decision"] == "moved", row
    assert row["new_seq_id"].split("#")[0] == f"chrS:{e.true_start}-{e.true_end}"


def test_most_moved_ends_land_on_their_tsd(ran):
    fx, _, _, _ = ran
    side = sidecar(fx)
    assert sum(side[fx.kind(k).key]["tsd_new"] not in (".", "NA") for k in MOVED_KINDS) >= 8


def test_the_second_call_of_one_element_is_merged_into_the_first(ran):
    fx, _, _, _ = ran
    a, b = fx.kind("split_short"), fx.kind("split_part")
    row = sidecar(fx)[b.key]
    new = sidecar(fx)[a.key]["new_seq_id"]
    assert (row["decision"], row["reason"], row["merged_into"]) == ("merged", "split_call", new)
    assert table_row(fx, 0, b.name) is None
    assert not fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], {b.name})


def test_a_clash_with_a_neighbour_is_rejected(ran):
    fx, _, _, _ = ran
    assert sidecar(fx)[fx.kind("conflict_left").key]["reason"] == "overlaps_element"


def _add_call(fx, e):
    """Append `e` to the depth-0 clean table and FASTA, as detection would have left it."""
    contig = "".join(fx.genome.read_text().splitlines()[1:])
    seq = contig[e.start - 1:e.end]
    base = fx.indir / f"{fx.prefix}_depth0_clean_ltr"
    with open(f"{base}.tsv", "a") as fh:
        fh.write("\t".join(_row(e)) + "\n")
    with open(f"{base}.fa", "a") as fh:
        fh.write(f">{e.name}\n" + "".join(seq[i:i + 60] + "\n"
                                          for i in range(0, len(seq), 60)))


def test_two_calls_moved_to_the_same_ends_do_not_abort_the_run(tmp_path, k2l_api, blastn):
    """One element called twice, one call inside the other: 3 bp short on the left, and
    2 bp long on the right. Extending the one and trimming the other gives both the
    element's true ends, and so one key; the rewrite refused that and the whole run wrote
    nothing. The first claim in genome order moves; the other keeps its call."""
    fx = build(tmp_path / "syn_converge")
    inner = fx.kind("short3_left")
    outer = Element("dup_long2_right", FAMILY, "+", inner.true_start, inner.true_end,
                    inner.true_start, inner.true_end + 2, inner.l1 + 2, inner.r0 - 3)
    _add_call(fx, outer)
    rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=2), TOOLS)
    side = sidecar(fx)
    assert side[outer.key]["decision"] == "moved", side[outer.key]
    assert side[outer.key]["new_seq_id"].split("#")[0] == \
        f"chrS:{inner.true_start}-{inner.true_end}"
    assert (side[inner.key]["decision"], side[inner.key]["reason"]) == \
        ("rejected", "duplicates_move")
    assert table_row(fx, 0, inner.name) is not None                # it keeps its call
    assert side[fx.kind("del_left").key]["decision"] == "moved"     # and the rest went ahead


def test_intact_and_tsd_bearing_copies_are_left_alone(ran):
    fx, _, _, _ = ran
    side = sidecar(fx)
    assert not [e.kind for e in fx.elements
                if e.kind in ("normal", "decayed", "host") and e.key in side]


def test_unpaired_bases_an_insertion_brings_in_do_not_age_the_element(ran, k2l_api):
    """ins_left's 40 bp remainder and 1.5 kb insertion have no partner in the other LTR:
    they are gap columns, so the K2P is the one the called pair measures."""
    fx, _, _, before = ran
    e = fx.kind("ins_left")
    name = e.name
    rec = [l for l in before[f"{fx.prefix}_depth0_clean_ltr.fa"].split(">") if l.startswith(name)]
    seq = "".join(rec[0].splitlines()[1:])
    spans = (0, e.l1 - e.start, e.r0 - e.start, e.end - e.start)
    ref = k2l_api.classify("x", seq, period_rule="outermost", spans=spans)
    new = table_row(fx, 0, f"chrS:{e.true_start}-{e.true_end}#LTR/Gypsy/Synth")
    assert abs(float(new["k2p"]) - ref.k2p) <= 0.005
    assert sidecar(fx)[e.key]["end_source5"] == "anchor"
    assert int(sidecar(fx)[e.key]["unpaired5"]) > 0


def test_rewritten_rows_come_from_kmer2ltr_and_keep_the_annotation(ran):
    fx, _, _, _ = ran
    e = fx.kind("del_left")
    row = table_row(fx, 0, f"chrS:{e.true_start}-{e.true_end}#LTR/Gypsy/Synth")
    assert row["ltr5_start"] == "1" and row["seq_len"] == str(e.true_end - e.true_start + 1)
    assert row["flank5_len"] == row["flank3_len"] == "0" and row["status"] == "pass"
    assert row["orientation"] == "+" and row["strand"] == "+" and row["family"] == e.family
    assert row["tsd"] not in (".", "NA") and row["cigar"] not in (".", "")


def test_a_minus_record_is_stored_reverse_complemented(ran):
    fx, _, _, _ = ran
    e = fx.kind("del_left_minus")
    new = f"chrS:{e.true_start}-{e.true_end}#LTR/Gypsy/Synth"
    genome = "".join(fx.genome.read_text().splitlines()[1:])
    rec = fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], {new})[new]
    assert rec == revcomp_record(genome[e.true_start - 1:e.true_end])


def test_a_trimmed_record_is_cut_to_its_new_ends(ran):
    fx, _, _, _ = ran
    e = fx.kind("long3_right")
    new = f"chrS:{e.true_start}-{e.true_end}#LTR/Gypsy/Synth"
    genome = "".join(fx.genome.read_text().splitlines()[1:])
    rec = fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], {new})[new]
    assert rec == genome[e.true_start - 1:e.true_end]
    assert sidecar(fx)[e.key]["ext3"] == "-3"


def test_the_host_masks_the_recovered_bases_and_points_at_the_new_name(ran):
    fx, _, _, _ = ran
    e, host = fx.kind("nested_del_left"), fx.kind("host")
    rec = fetch_records([str(fx.indir / f"{fx.prefix}_depth1_clean_ltr.fa")],
                        {host.name})[host.name]
    assert set(rec[e.true_start - host.start:e.start - host.start]) == {"N"}
    t1 = (fx.indir / f"{fx.prefix}_depth1_clean_ltr.tsv").read_text()
    assert f"nest-outer:chrS:{e.true_start}-{e.true_end}" in t1


def test_raw_tables_are_untouched(ran):
    fx, _, raw, _ = ran
    for name, data in raw.items():
        assert (fx.indir / name).read_bytes() == data, name


def test_a_merge_is_undone_when_the_survivor_fails_kmer2ltr(tmp_path, k2l_api, blastn,
                                                            monkeypatch):
    fx = build(tmp_path / "syn_merge_fail")
    a, b = fx.kind("split_short"), fx.kind("split_part")
    real = rb.arbitrate

    def failing(args):
        v = real(args)
        return (rb.Verdict(v.uid, "kmer2ltr_not_pass", None, "weak_pair", "NA", "NA")
                if args[0].member.name == a.name else v)

    monkeypatch.setattr(rb, "arbitrate", failing)
    monkeypatch.setattr(rb, "Pool", _SerialPool)
    rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=1), TOOLS)
    assert sidecar(fx)[a.key]["reason"] == "kmer2ltr_not_pass"
    assert table_row(fx, 0, b.name) is not None                   # the partner stays


def test_a_call_absorbed_by_a_survivor_that_drops_out_before_kmer2ltr_says_why(
        tmp_path, k2l_api, blastn, monkeypatch):
    """An absorbed merge partner skips its own move. When its survivor then drops out
    before Kmer2LTR (here: the survivor's record is missing), the partner keeps its call
    and its row says why -- never a rejected row with no reason."""
    fx = build(tmp_path / "syn_absorbed")
    a, b = sorted((fx.kind("del_left"), fx.kind("patch_right")), key=lambda e: e.start)
    real_mp, real_fetch = rb.merge_partner, rb.fetch_records

    def forced(index, by_uid, m, left, right, tol=5):
        return b.key if m.name == a.name else real_mp(index, by_uid, m, left, right, tol)

    def without_a(paths, names):
        return {k: v for k, v in real_fetch(paths, names).items() if k != a.name}

    monkeypatch.setattr(rb, "merge_partner", forced)
    monkeypatch.setattr(rb, "fetch_records", without_a)
    monkeypatch.setattr(rb, "Pool", _SerialPool)
    rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=1), TOOLS)
    side = sidecar(fx)
    assert side[a.key]["reason"] == "record_missing"
    assert (side[b.key]["decision"], side[b.key]["reason"]) == ("rejected", "merge_partner_failed")
    assert all(r["reason"] != "." for r in side.values() if r["decision"] == "rejected")
    assert table_row(fx, 0, b.name) is not None                   # the partner stays


def test_a_survivor_is_never_itself_absorbed(tmp_path, k2l_api, blastn, monkeypatch):
    """A call that absorbs a partner must not be absorbed in turn: its partner would be
    neither retired nor released, left rejected with no reason, and its own move lost."""
    fx = build(tmp_path / "syn_chain")
    y, x, s = sorted((fx.kind("del_left"), fx.kind("patch_right"), fx.kind("short3_left")),
                     key=lambda e: e.start)
    real = rb.merge_partner

    def chained(index, by_uid, m, left, right, tol=5):
        return {x.name: y.key, s.name: x.key}.get(m.name) or real(index, by_uid, m, left,
                                                                    right, tol)

    monkeypatch.setattr(rb, "merge_partner", chained)
    monkeypatch.setattr(rb, "Pool", _SerialPool)
    rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=1), TOOLS)
    side = sidecar(fx)
    assert all(r["reason"] != "." for r in side.values() if r["decision"] == "rejected")
    assert side[y.key]["decision"] == "merged"
    assert side[y.key]["merged_into"] == side[x.key]["new_seq_id"]
    assert side[s.key]["decision"] == "moved"


class _SerialPool:
    """Pool stand-in running in-process, so monkeypatched workers are the ones used."""

    def __init__(self, n, initializer=None, initargs=()):
        if initializer:
            initializer(*initargs)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def imap_unordered(self, f, jobs, chunksize=1):
        return map(f, jobs)

    def map(self, f, jobs, chunksize=1):
        return list(map(f, jobs))


def test_kmer2ltr_is_handed_the_called_pair(tmp_path, k2l_api, blastn, monkeypatch):
    fx = build(tmp_path / "syn_spans")
    e = fx.kind("patch_right")                       # a pure extension on the right
    calls = []

    def spy(name, seq, **kw):
        calls.append((name, kw.get("spans"), kw.get("tsd_credit")))
        return k2l_api.classify(name, seq, **kw)

    monkeypatch.setattr(rb.k2l, "api", lambda *a, **k: k2l_api._replace(classify=spy))
    monkeypatch.setattr(rb, "Pool", _SerialPool)
    rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=1), TOOLS)
    got = [c for c in calls if c[0].startswith(f"chrS:{e.start}-{e.true_end}#")]
    assert got and got[0][1] == (0, e.l1 - e.start, e.r0 - e.start, e.end - e.start)
    assert got[0][2] == rb.CREDIT_BITS


def test_a_second_run_moves_nothing_and_the_cli_works(tmp_path, k2l_api, blastn):
    fx = build(tmp_path / "syn_twice")
    args = ["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
            "-t", "2", "--tools-dir", TOOLS]
    assert rb.main(args) == 0
    assert sum(r["decision"] == "moved" for r in sidecar(fx).values()) == len(MOVED_KINDS)
    counts = rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=2), TOOLS)
    assert counts.get("moved", 0) == 0 and counts.get("merged", 0) == 0


def test_a_prefix_with_no_clean_tables_is_skipped_not_fatal(tmp_path, k2l_api, blastn):
    fx = build(tmp_path / "syn_partial")
    counts = rb.run(str(fx.indir), [fx.prefix, "missing_p"], [str(fx.genome), str(fx.genome)],
                    rb.Settings(threads=2), TOOLS)
    assert counts.get("moved", 0) == len(MOVED_KINDS)
    assert not (fx.indir / "missing_p_reboundary.tsv").exists()


def test_all_prefixes_empty_raises_systemexit(tmp_path, k2l_api, blastn):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    with pytest.raises(SystemExit, match="no _clean_ depth tables"):
        rb.run(str(tmp_path), ["missing"], [str(genome)], rb.Settings(threads=1), TOOLS)


def test_a_run_with_no_templates_changes_nothing_and_says_so(tmp_path, k2l_api, blastn, capsys):
    fx = build(tmp_path / "syn_none")
    for path in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.tsv"):     # no call keeps a TSD
        lines = path.read_text().splitlines()
        head = lines[0].lstrip("#").split("\t")
        i = head.index("tsd")
        out = [lines[0]] + ["\t".join(c if k != i else "." for k, c in enumerate(l.split("\t")))
                            for l in lines[1:]]
        path.write_text("\n".join(out) + "\n")
    before = clean_bytes(fx)
    counts = rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=1), TOOLS)
    assert counts == {} and "no templates" in capsys.readouterr().err
    after = clean_bytes(fx)
    assert {k: v for k, v in after.items() if "reboundary" not in k} == before
    assert sidecar(fx) == {}


def test_a_missing_blastn_fails_before_any_work(tmp_path):
    fx = build(tmp_path / "syn_noblast")
    before = clean_bytes(fx)
    with pytest.raises(SystemExit, match="blastn not found"):
        rb.main(["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
                 "--blastn", str(tmp_path / "no_such_blastn"), "--posthoc"])
    assert clean_bytes(fx) == before
    assert not (fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}").exists()


def test_an_old_kmer2ltr_is_refused_with_its_own_exit_code(tmp_path, blastn, monkeypatch):
    fx = build(tmp_path / "syn_oldk2l")

    def old(*args, **kwargs):
        raise rb.k2l.IncompatibleKmer2LTR("force-pairs credited flanks")

    monkeypatch.setattr(rb.k2l, "api", old)
    before = clean_bytes(fx)
    assert rb.main(["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome",
                    str(fx.genome), "--posthoc"]) == rb.EXIT_INCOMPATIBLE_KMER2LTR
    assert clean_bytes(fx) == before
    assert not (fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}").exists()


def test_posthoc_backs_up_once_reruns_from_originals_and_restores(tmp_path, k2l_api, blastn):
    fx = build(tmp_path / "syn_posthoc")
    originals = {p.name: p.read_bytes() for p in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.*")}
    (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").write_text("##gff-version 3\n#old\n")
    (fx.indir / f"{fx.prefix}_plots").mkdir()
    (fx.indir / f"{fx.prefix}_plots" / "x.pdf").write_text("old plot")
    base = ["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
            "-t", "2", "--tools-dir", TOOLS, "--posthoc", "--no-plots"]

    assert rb.main(base) == 0
    backup = fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}"
    assert {p.name for p in backup.iterdir()} >= set(originals) | {
        f"{fx.prefix}_all_depth_LTR_cleaned.gff3", f"{fx.prefix}_plots"}
    gff = (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").read_text()
    assert "boundary_source=templates" in gff and "#old" not in gff
    assert not (fx.indir / f"{fx.prefix}_plots").exists()

    assert rb.main(base + ["--no-anchor"]) == 0      # starts again from the backup
    e = fx.kind("del_left")
    assert sidecar(fx)[e.key]["decision"] == "moved"

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


def test_posthoc_refuses_tables_a_pipeline_run_already_re_boundaried(tmp_path):
    """A pipeline run re-boundaries in place: its sidecar is there, no backup is. A
    post-hoc run would back up the moved calls as "originals" and move them again,
    which is less precise than the first pass."""
    fx = build(tmp_path / "syn_already")
    (fx.indir / f"{fx.prefix}_reboundary.tsv").write_text("#old_seq_id\tnew_seq_id\n")
    before = clean_bytes(fx)
    with pytest.raises(SystemExit, match="already re-boundaried"):
        rb.prepare_posthoc(str(fx.indir), fx.prefix)
    assert clean_bytes(fx) == before
    assert not (fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}").exists()


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


def test_a_repeated_prefix_is_refused_before_anything_is_written(tmp_path, blastn):
    """One prefix twice would stage every one of its files twice and lose the originals."""
    fx = build(tmp_path / "syn_dupe")
    before = clean_bytes(fx)
    with pytest.raises(SystemExit, match="more than once") as got:
        rb.run(str(fx.indir), [fx.prefix, fx.prefix], [str(fx.genome)] * 2,
               rb.Settings(threads=1), TOOLS)
    assert fx.prefix in str(got.value)
    assert clean_bytes(fx) == before
    assert not (fx.indir / f"{fx.prefix}{rb.SIDECAR_SUFFIX}").exists()


def test_a_repeated_prefix_is_refused_by_the_cli_before_the_backup(tmp_path):
    fx = build(tmp_path / "syn_dupe_cli")
    before = clean_bytes(fx)
    with pytest.raises(SystemExit, match="more than once"):
        rb.main(["--indir", str(fx.indir), "--prefix", fx.prefix, fx.prefix,
                 "--genome", str(fx.genome), str(fx.genome), "--posthoc"])
    assert clean_bytes(fx) == before
    assert not (fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}").exists()


def test_a_leftover_from_an_interrupted_swap_stops_the_next_run(tmp_path, blastn):
    """The killed swap took the target away, so only the directory itself can reveal it."""
    fx = build(tmp_path / "syn_leftover")
    target = fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv"
    (fx.indir / (target.name + OLD_SUFFIX)).write_bytes(target.read_bytes())
    target.unlink()                                  # killed between the two renames
    with pytest.raises(SystemExit, match="interrupted") as got:
        rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)],
               rb.Settings(threads=1), TOOLS)
    assert target.name + OLD_SUFFIX in str(got.value)


def test_a_leftover_new_file_stops_the_cli_before_it_backs_anything_up(tmp_path):
    fx = build(tmp_path / "syn_leftover_new")
    left = fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa{NEW_SUFFIX}"
    left.write_text(">x\nACGT\n")
    before = clean_bytes(fx)
    with pytest.raises(SystemExit, match="interrupted") as got:
        rb.main(["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
                 "--posthoc"])
    assert left.name in str(got.value)
    assert clean_bytes(fx) == before
    assert not (fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}").exists()


def test_restore_names_the_derived_output_it_cannot_replace(tmp_path, capsys):
    """A GFF3 built from extended coordinates must not sit silently beside restored tables."""
    fx = build(tmp_path / "syn_restore_warn")
    rb.prepare_posthoc(str(fx.indir), fx.prefix)          # no GFF3 yet, so none is backed up
    gff = fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3"
    gff.write_text("##gff-version 3\n")                   # as regenerate() would leave it
    capsys.readouterr()
    rb.restore(str(fx.indir), fx.prefix)
    err = capsys.readouterr().err
    assert gff.exists()                                   # the backup cannot replace it
    assert gff.name in err and "WARNING" in err


def test_sidecar_rows_are_ordered_without_ties():
    rows = [{"old_seq_id": "c:100-900#x"}, {"old_seq_id": "c:100-500#x"}]
    assert sorted(rows, key=rb._row_order) == sorted(rows[::-1], key=rb._row_order)


def test_verbose_adds_per_step_progress(tmp_path, k2l_api, blastn, capsys):
    out = {}
    for flag in ([], ["-v"]):
        fx = build(tmp_path / f"syn_v{len(flag)}")
        assert rb.main(["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome",
                        str(fx.genome), "-t", "2", "--tools-dir", TOOLS] + flag) == 0
        out[bool(flag)] = capsys.readouterr().out
    for marker in ("templates (after the per-family cap)", "side outcomes", "Kmer2LTR status"):
        assert marker in out[True] and marker not in out[False]


def test_a_soft_masked_record_is_checked_against_its_row_too():
    """Detection keeps the genome's case, so a soft-masked genome gives lowercase
    records: the check must compare them, not skip every lowercase base."""
    genomic = "ACGTTGCAAGGCTTAACCGTAGGATCCA"
    flipped = genomic[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert rb.record_matches(genomic.lower(), genomic)
    assert rb.record_matches(genomic.lower()[:5] + "RRRR" + genomic.lower()[9:], genomic)
    assert not rb.record_matches(flipped.lower(), genomic)


def test_a_record_that_does_not_match_its_row_is_not_spliced(tmp_path, k2l_api, blastn):
    """A record stored reverse-complemented under an orientation `+` row would get genome
    segments spliced onto a flipped core. It is refused, never rewritten."""
    fx = build(tmp_path / "syn_flipped")
    e = fx.kind("del_left")
    path = fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa"
    out, name = [], None
    recs = {}
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            recs[name] = ""
        else:
            recs[name] += line
    recs[e.name] = revcomp_record(recs[e.name])
    path.write_text("".join(f">{n}\n{s}\n" for n, s in recs.items()))
    rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=2), TOOLS)
    row = sidecar(fx)[e.key]
    assert (row["decision"], row["reason"]) == ("rejected", "record_mismatch")
    assert fetch_records([str(path)], {e.name})[e.name] == recs[e.name]     # left as it was

"""Reading and rewriting a run's clean side: nothing raw, nothing half-written."""

from __future__ import annotations

import dataclasses
import os

import pytest

from ltrquest import reboundary_io as rio
from ltrquest.detect import revcomp as revcomp_record
from ltrquest.kmer2ltr import COLUMNS
from ltrquest.reconcile import IUPAC_DEPTH_SEQ
from reboundary_fixtures import build, members


@pytest.fixture
def fx(tmp_path):
    return build(tmp_path / "syn")


def test_load_clean_tables_finds_both_depths_and_their_fastas(fx):
    tables = rio.load_clean_tables(str(fx.indir), fx.prefix)
    assert [t.depth for t in tables] == [0, 1]
    assert all(os.path.isfile(t.fasta) for t in tables)


def test_load_refuses_a_table_that_was_never_annotated(fx):
    path = fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv"
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    cut = head.index("family")
    path.write_text("\n".join("\t".join(l.split("\t")[:cut] + l.split("\t")[cut + 1:])
                              for l in lines) + "\n")
    with pytest.raises(ValueError, match="no family column.*needs the strand, family, nesting"):
        rio.load_clean_tables(str(fx.indir), fx.prefix)


def test_members_from_reads_called_coordinates(fx):
    got = {m.name: m for m in rio.members_from(fx.prefix, rio.load_clean_tables(str(fx.indir),
                                                                                 fx.prefix))}
    e = fx.kind("del_left")
    m = got[e.name]
    assert (m.start, m.end, m.l1, m.r0, m.strand, m.depth) == (e.start, e.end, e.l1, e.r0, "+", 0)
    assert got[fx.kind("host").name].depth == 1 and len(got) == len(fx.elements)


def test_fetch_records_returns_only_the_named(fx):
    names = {fx.kind("del_left").name}
    got = rio.fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], names)
    assert set(got) == names


def test_forward_and_stored_leave_depth_letters_alone():
    rec = "ACGTNRDYACGT"
    assert rio.forward(rec, "-") == revcomp_record(rec)
    assert rio.stored(rio.forward(rec, "-"), "-") == rec
    assert rio.forward(rec, "+") == rec


def test_sanitize_turns_every_ambiguity_into_n():
    assert rio.sanitize("acgtNRDY") == "ACGTNNNN"


def test_paint_mirrors_on_reverse_stored_records():
    rec = "A" * 20
    assert rio.paint(rec, "+", 101, [(103, 105)], "R") == "AA" + "RRR" + "A" * 15
    assert rio.paint(rec, "-", 101, [(103, 105)], "R") == "A" * 15 + "RRR" + "AA"


def test_rekey_nest_and_hosts():
    value = "nest-outer:c:1-9;nest-inner:c:50-900"
    assert rio.rekey_nest(value, {"c:1-9": "c:1-12"}) == "nest-outer:c:1-12;nest-inner:c:50-900"
    assert rio.hosts_of(value) == ["c:50-900"]
    assert rio.rekey_nest(".", {"c:1-9": "c:1-12"}) == "."


def test_conflict_rules(fx):
    ms = {m.name: m for m in members(fx)}
    index = rio.SpanIndex(ms.values())
    inner = ms[fx.kind("nested_del_left").name]
    host = ms[fx.kind("host").name]
    assert rio.conflict(index, inner, inner.start - 50, inner.end) is None
    assert rio.conflict(index, inner, host.start - 10, inner.end) == "host_exceeded"
    c = ms[fx.kind("conflict_left").name]
    assert rio.conflict(index, c, fx.kind("conflict_left").true_start, c.end) == "overlaps_element"
    x = fx.kind("neighbour")
    assert rio.conflict(index, c, x.start - 5, c.end) == "engulfs_element"


def test_mutual_conflicts_keep_the_first_claim(fx):
    a = members(fx)[0]
    b = dataclasses.replace(a, name="chrS:900000-900500#x", start=a.end + 100, end=a.end + 600,
                            l1=a.end + 200, r0=a.end + 500)
    lost = rio.mutual_conflicts([(a, a.start, a.end + 60), (b, b.start - 60, b.end)])
    assert lost == {b.uid}


def test_sidecar_round_trips_through_read_map(tmp_path):
    row = dict.fromkeys(rio.SIDECAR_COLUMNS, ".")
    row.update(old_seq_id="c:150-900#LTR/Gypsy/X", new_seq_id="c:100-900#LTR/Gypsy/X",
               decision="extended", ext5="50", ext3="0")
    rejected = dict(row, new_seq_id=".", decision="rejected", reason="gate_identity",
                    old_seq_id="c:2000-2500#LTR/Gypsy/X")
    path = tmp_path / "p_reboundary.tsv"
    path.write_text(rio.sidecar_text([row, rejected]))
    got = rio.read_map(str(path))
    assert set(got) == {"c:100-900"}
    assert got["c:100-900"] == rio.Rebound("c:150-900#LTR/Gypsy/X", "c:100-900#LTR/Gypsy/X", 50, 0)


def _leftovers(tmp_path):
    return sorted(p.name for p in tmp_path.iterdir()
                  if p.name.endswith(rio.NEW_SUFFIX) or p.name.endswith(rio.OLD_SUFFIX))


def test_commit_is_all_or_nothing(tmp_path):
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    a.write_text("old-a")
    c = rio.Commit()
    with c.open(str(a)) as fh:
        fh.write("new-a")
    with c.open(str(b)) as fh:
        fh.write("new-b")
    assert a.read_text() == "old-a" and not b.exists()
    c.abort()
    assert a.read_text() == "old-a" and not b.exists() and not _leftovers(tmp_path)
    c2 = rio.Commit()
    with c2.open(str(a)) as fh:
        fh.write("new-a")
    c2.commit()
    assert a.read_text() == "new-a" and not _leftovers(tmp_path)


def test_a_commit_that_fails_part_way_puts_every_original_back(tmp_path, monkeypatch):
    """A kill between two renames must not leave a table rewritten and its FASTA not."""
    a, b, c = tmp_path / "a.txt", tmp_path / "b.txt", tmp_path / "c.txt"
    a.write_text("old-a")
    c.write_text("old-c")                      # b does not exist yet
    commit = rio.Commit()
    for path, text in ((a, "new-a"), (b, "new-b"), (c, "new-c")):
        with commit.open(str(path)) as fh:
            fh.write(text)

    real, seen = os.replace, []

    def flaky(src, dst):
        seen.append(src)
        if len(seen) == 4:                     # part way through the rename phase
            raise OSError(28, "No space left on device")
        return real(src, dst)

    monkeypatch.setattr(rio.os, "replace", flaky)
    with pytest.raises(OSError):
        commit.commit()
    monkeypatch.undo()

    assert a.read_text() == "old-a"
    assert not b.exists()                      # created by the commit, so taken away again
    assert c.read_text() == "old-c"
    assert not _leftovers(tmp_path)


def test_a_rollback_that_cannot_finish_says_exactly_what_is_where(tmp_path, monkeypatch):
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    a.write_text("old-a")
    b.write_text("old-b")
    commit = rio.Commit()
    for path, text in ((a, "new-a"), (b, "new-b")):
        with commit.open(str(path)) as fh:
            fh.write(text)

    real, seen = os.replace, []

    def flaky(src, dst):
        seen.append(src)
        if len(seen) >= 3:                     # the swap fails and so does putting it back
            raise OSError(5, "Input/output error")
        return real(src, dst)

    monkeypatch.setattr(rio.os, "replace", flaky)
    with pytest.raises(rio.CommitError) as got:
        commit.commit()
    monkeypatch.undo()
    assert str(a) + rio.OLD_SUFFIX in str(got.value)
    assert os.path.isfile(str(a) + rio.OLD_SUFFIX)     # nothing was thrown away
    assert (tmp_path / "b.txt").read_text() == "old-b"


def test_an_original_left_by_a_killed_swap_is_never_overwritten(tmp_path):
    a = tmp_path / "a.txt"
    a.write_text("half-committed")
    (tmp_path / ("a.txt" + rio.OLD_SUFFIX)).write_text("the original")
    commit = rio.Commit()
    with commit.open(str(a)) as fh:
        fh.write("newer still")
    with pytest.raises(rio.CommitError, match="interrupted mid-swap"):
        commit.commit()
    commit.abort()
    assert (tmp_path / ("a.txt" + rio.OLD_SUFFIX)).read_text() == "the original"
    assert a.read_text() == "half-committed"


def test_staging_one_path_twice_is_refused(tmp_path):
    """The second swap would move the first swap's rewrite into the '.old' holding the original."""
    a = tmp_path / "a.txt"
    a.write_text("the original")
    commit = rio.Commit()
    with commit.open(str(a)) as fh:
        fh.write("rewrite-1")
    with pytest.raises(ValueError, match="staged twice"):
        commit.open(str(a))
    commit.abort()
    assert a.read_text() == "the original" and not _leftovers(tmp_path)


def test_a_leftover_is_found_even_when_its_target_is_gone(tmp_path):
    """A kill between the two renames leaves no target to discover, only the '.old'."""
    (tmp_path / ("p_depth0_clean_ltr.tsv" + rio.OLD_SUFFIX)).write_text("the original")
    (tmp_path / ("p_depth1_clean_ltr.fa" + rio.NEW_SUFFIX)).write_text(">x\nACGT\n")
    (tmp_path / ("p" + rio.SIDECAR_SUFFIX + rio.NEW_SUFFIX)).write_text("#\n")
    (tmp_path / "genome.fa.old").write_text("not ours")
    got = [os.path.basename(p) for p in rio.leftover_staging(str(tmp_path))]
    assert got == ["p_depth0_clean_ltr.tsv" + rio.OLD_SUFFIX,
                   "p_depth1_clean_ltr.fa" + rio.NEW_SUFFIX,
                   "p" + rio.SIDECAR_SUFFIX + rio.NEW_SUFFIX]
    assert rio.leftover_staging(str(tmp_path / "gone")) == []


def test_a_proposal_that_fills_its_host_exactly_is_rejected(fx):
    """Extending to the host's own span would re-key the inner onto the host's key."""
    ms = {m.name: m for m in members(fx)}
    index = rio.SpanIndex(ms.values())
    inner = ms[fx.kind("nested_del_left").name]
    host = ms[fx.kind("host").name]
    assert rio.conflict(index, inner, host.start, host.end) == "duplicates_host"
    assert rio.conflict(index, inner, host.start + 1, host.end) is None


def test_rewrite_refuses_two_records_that_would_share_a_key(fx):
    tables = rio.load_clean_tables(str(fx.indir), fx.prefix)
    ms = {m.name: m for m in rio.members_from(fx.prefix, tables)}
    e, host = fx.kind("nested_del_left"), fx.kind("host")
    m = ms[e.name]
    new_name = f"chrS:{host.start}-{host.end}#LTR/Gypsy/Synth"      # the host's own key
    old_row = [r for r in tables[0].rows if r[0] == e.name][0]
    fields = [new_name] + old_row[1:len(COLUMNS)]
    acc = rio.Accepted(m, new_name, fields, "A" * (host.end - host.start + 1),
                       host.start, host.end)
    before = {p.name: p.read_bytes() for p in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.*")}
    c = rio.Commit()
    with pytest.raises(ValueError, match="same key"):
        rio.rewrite(tables, {m.key: acc}, IUPAC_DEPTH_SEQ, c)
    c.abort()
    for name, data in before.items():
        assert (fx.indir / name).read_bytes() == data, name


def test_rewrite_renames_rekeys_and_paints_the_host(fx):
    tables = rio.load_clean_tables(str(fx.indir), fx.prefix)
    ms = {m.name: m for m in rio.members_from(fx.prefix, tables)}
    e = fx.kind("nested_del_left")
    m = ms[e.name]
    new_name = f"chrS:{e.true_start}-{e.end}#LTR/Gypsy/Synth"
    old_row = [r for r in tables[0].rows if r[0] == e.name][0]
    fields = [new_name] + old_row[1:len(COLUMNS)]
    genome = "".join(fx.genome.read_text().splitlines()[1:])
    record = genome[e.true_start - 1:e.end]
    acc = rio.Accepted(m, new_name, fields, record, e.true_start, e.end)
    raw_before = (fx.indir / f"{fx.prefix}_depth0_ltr.tsv").read_bytes()
    c = rio.Commit()
    rio.rewrite(tables, {m.key: acc}, IUPAC_DEPTH_SEQ, c)
    c.commit()

    t0 = (fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv").read_text()
    assert new_name in t0 and e.name not in t0
    t1 = (fx.indir / f"{fx.prefix}_depth1_clean_ltr.tsv").read_text()
    assert f"nest-outer:chrS:{e.true_start}-{e.end}" in t1
    fa0 = rio.fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], {new_name})
    assert fa0[new_name] == record
    host = fx.kind("host")
    fa1 = rio.fetch_records([str(fx.indir / f"{fx.prefix}_depth1_clean_ltr.fa")], {host.name})
    painted = fa1[host.name][e.true_start - host.start:e.end - host.start + 1]
    assert set(painted) == {"N"}
    assert (fx.indir / f"{fx.prefix}_depth0_ltr.tsv").read_bytes() == raw_before


# ---------------------------------------------------------------- trims, merges (v2)

def _mini_chain(tmp_path):
    """A 3-level nest written as the pipeline writes it: H2 (depth 2, stored +) holds
    H1 (depth 1, stored -), which holds E (depth 0, stored +)."""
    import random
    from reboundary_fixtures import TAIL, Element, _row
    r = random.Random(5)
    seq = "".join(r.choice("ACGT") for _ in range(6000))
    (tmp_path / "g.fa").write_text(">chrS\n" + seq + "\n")
    h2 = Element("h2", "fam_a", "+", 101, 5900, 101, 5900, 400, 5601, depth=2)
    h1 = Element("h1", "fam_b", "-", 1001, 4900, 1001, 4900, 1300, 4601, depth=1)
    e = Element("e", "fam_c", "+", 2001, 3900, 2001, 3900, 2300, 3601, depth=0)
    h2.nest_status = f"nest-outer:{h1.key};nest-outer:{e.key}"
    h1.nest_status = f"nest-outer:{e.key};nest-inner:{h2.key}"
    e.nest_status = f"nest-inner:{h1.key};nest-inner:{h2.key}"

    def painted(x, kids):
        rec = list(seq[x.start - 1:x.end])
        for k in kids:                                     # parent first, deeper after
            for p in range(k.start, k.end + 1):
                rec[p - x.start] = IUPAC_DEPTH_SEQ[k.depth]
        rec = "".join(rec)
        return revcomp_record(rec) if x.strand == "-" else rec

    records = {e.key: painted(e, []), h1.key: painted(h1, [e]), h2.key: painted(h2, [h1, e])}
    header = "#" + "\t".join(COLUMNS + TAIL)
    for x in (e, h1, h2):
        base = tmp_path / f"p_depth{x.depth}_clean_ltr"
        (base.with_suffix(".tsv")).write_text(header + "\n" + "\t".join(_row(x)) + "\n")
        (base.with_suffix(".fa")).write_text(f">{x.name}\n{records[x.key]}\n")
    return seq, e, h1, h2


def test_a_trim_gives_each_host_back_what_reconcile_would_have_painted(tmp_path):
    seq, e, h1, h2 = _mini_chain(tmp_path)
    tables = rio.load_clean_tables(str(tmp_path), "p")
    ms = {m.uid: m for m in rio.members_from("p", tables)}
    m = next(x for x in ms.values() if x.key == e.key)
    restore = rio.trim_restore(m, [(3898, 3900)], ms, IUPAC_DEPTH_SEQ,
                               lambda prefix, chrom, a, b: seq[a - 1:b])
    new_name = f"chrS:2001-3897#LTR/Gypsy/Synth"
    old_row = tables[0].rows[0]
    acc = rio.Accepted(m, new_name, [new_name] + old_row[1:len(COLUMNS)], seq[2000:3897],
                       2001, 3897, restore=restore)
    c = rio.Commit()
    rio.rewrite(tables, {m.key: acc}, IUPAC_DEPTH_SEQ, c)
    c.commit()
    fa1 = rio.fetch_records([str(tmp_path / "p_depth1_clean_ltr.fa")], {h1.name})[h1.name]
    fa2 = rio.fetch_records([str(tmp_path / "p_depth2_clean_ltr.fa")], {h2.name})[h2.name]
    fwd1 = rio.forward(fa1, "-")
    assert fwd1[3898 - 1001:3900 - 1001 + 1] == seq[3897:3900]         # parent: the genome
    assert set(fwd1[2001 - 1001:3897 - 1001 + 1]) == {"N"}              # the call itself: still N
    assert fa2[3898 - 101:3900 - 101 + 1] == "RRR"                       # grandparent: H1's R
    assert set(fa2[2001 - 101:3897 - 101 + 1]) == {"N"}


def test_a_prebuilt_host_index_restores_exactly_what_a_member_scan_does(tmp_path):
    """At 10^5 members a scan of every member per trimmed nested element and host costs
    minutes; the index built once must give the same restore."""
    seq, e, h1, h2 = _mini_chain(tmp_path)
    tables = rio.load_clean_tables(str(tmp_path), "p")
    ms = {m.uid: m for m in rio.members_from("p", tables)}
    m = next(x for x in ms.values() if x.key == e.key)

    def fetch(prefix, chrom, a, b):
        return seq[a - 1:b]
    scan = rio.trim_restore(m, [(3898, 3900)], ms, IUPAC_DEPTH_SEQ, fetch)
    index = rio.nested_index(ms.values())
    assert rio.trim_restore(m, [(3898, 3900)], ms, IUPAC_DEPTH_SEQ, fetch, nested=index) == scan


def _split(fx):
    ms = {m.uid: m for m in members(fx)}
    by_name = {m.name: m for m in ms.values()}
    a, b = by_name[fx.kind("split_short").name], by_name[fx.kind("split_part").name]
    return ms, rio.SpanIndex(ms.values()), a, b, fx.kind("split_short")


def test_the_second_call_of_one_element_is_found_as_its_merge_partner(fx):
    ms, index, a, b, e = _split(fx)
    assert rio.merge_partner(index, ms, a, a.start, e.true_end) == b.key
    assert rio.merge_partner(index, ms, a, a.start, e.true_end - 20) is None    # end not shared


@pytest.mark.parametrize("change", [
    {"strand": "-"}, {"depth": 1},
    {"nest_status": "nest-inner:chrS:1-99999999"},
    {"nest_status": "nest-outer:chrS:1-2"}, {"tsd": "ACGTA", "tsd_offset": "0,0"},
])
def test_a_second_call_that_is_not_the_same_element_is_no_merge_partner(fx, change):
    ms, _, a, b, e = _split(fx)
    other = dataclasses.replace(b, **change)
    ms = dict(ms)
    ms[b.uid] = other
    assert rio.merge_partner(rio.SpanIndex(ms.values()), ms, a, a.start, e.true_end) is None


def test_a_second_call_in_another_family_is_still_a_merge_partner(fx):
    """A split call's family is clustered from its own, partly wrong, LTR pair, so it
    says nothing about which element the call belongs to."""
    ms, _, a, b, e = _split(fx)
    ms = dict(ms)
    ms[b.uid] = dataclasses.replace(b, family="other_fam")
    assert rio.merge_partner(rio.SpanIndex(ms.values()), ms, a, a.start, e.true_end) == b.key


def test_a_split_pair_inside_a_common_host_is_still_a_merge_partner():
    """The host spans the added bases too; it is the element's host, not a second
    call."""
    from ltrquest.ltr_model import Member

    def mm(s, e, name, nest, depth):
        return Member(prefix="p", name=f"c:{s}-{e}#{name}", chrom="c", start=s, end=e,
                      l1=s + 399, r0=e - 399, strand="+", orientation="+", family="f",
                      depth=depth, k2p=0.01, tsd=".", nest_status=nest)
    host = mm(1, 30000, "h", "nest-outer:c:5000-15000;nest-outer:c:8000-18000", 0)
    a = mm(5000, 15000, "a", "nest-inner:c:1-30000", 1)
    b = mm(8000, 18000, "b", "nest-inner:c:1-30000", 1)
    ms = {m.uid: m for m in (host, a, b)}
    assert rio.merge_partner(rio.SpanIndex(ms.values()), ms, a, 5000, 18000) == b.key


def test_a_retired_call_leaves_no_row_record_or_reference(fx):
    tables = rio.load_clean_tables(str(fx.indir), fx.prefix)
    ms = {m.name: m for m in rio.members_from(fx.prefix, tables)}
    e, x = fx.kind("split_short"), fx.kind("split_part")
    m = ms[e.name]
    new_name = f"chrS:{e.start}-{e.true_end}#LTR/Gypsy/Synth"
    old_row = [r for r in tables[0].rows if r[0] == e.name][0]
    genome = "".join(fx.genome.read_text().splitlines()[1:])
    acc = rio.Accepted(m, new_name, [new_name] + old_row[1:len(COLUMNS)],
                       genome[e.start - 1:e.true_end], e.start, e.true_end)
    c = rio.Commit()
    rio.rewrite(tables, {m.key: acc}, IUPAC_DEPTH_SEQ, c, retired=frozenset({ms[x.name].key}))
    c.commit()
    t0 = (fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv").read_text()
    assert new_name in t0 and x.name not in t0 and e.name not in t0
    fa0 = rio.fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")],
                            {x.name, new_name})
    assert set(fa0) == {new_name}


def test_a_retired_calls_nest_tokens_are_dropped_not_rekeyed():
    value = "nest-outer:c:1-9;nest-outer:c:5-20;nest-inner:c:1-900"
    got = rio.rekey_nest(rio.drop_nest(value, frozenset({"c:5-20"})), {"c:1-9": "c:1-20"})
    assert got == "nest-outer:c:1-20;nest-inner:c:1-900"
    assert rio.drop_nest("nest-outer:c:5-20", frozenset({"c:5-20"})) == "."


def test_a_survivor_may_take_its_retired_partners_key(fx):
    """The whole element can be one of the two calls: the survivor then takes that key,
    which is only legal because the partner leaves."""
    tables = rio.load_clean_tables(str(fx.indir), fx.prefix)
    ms = {m.name: m for m in rio.members_from(fx.prefix, tables)}
    e, x = fx.kind("split_short"), fx.kind("split_part")
    m, xm = ms[e.name], ms[x.name]
    new_name = x.name                                         # the partner's own key
    old_row = [r for r in tables[0].rows if r[0] == e.name][0]
    acc = rio.Accepted(m, new_name, [new_name] + old_row[1:len(COLUMNS)], "A" * 10,
                       xm.start, xm.end)
    c = rio.Commit()
    with pytest.raises(ValueError, match="same key"):
        rio.rewrite(tables, {m.key: acc}, IUPAC_DEPTH_SEQ, c)
    c.abort()
    c = rio.Commit()
    rio.rewrite(tables, {m.key: acc}, IUPAC_DEPTH_SEQ, c, retired=frozenset({xm.key}))
    c.abort()


def test_mutual_conflicts_do_not_depend_on_input_order_when_starts_tie():
    """Workers finish in any order; two calls sharing a start must not let that order
    pick which one keeps the bases both claim."""
    from ltrquest.ltr_model import Member

    def mm(s, e):
        return Member(prefix="p", name=f"c:{s}-{e}#x", chrom="c", start=s, end=e, l1=s + 99,
                      r0=e - 99, strand="+", orientation="+", family="f", depth=0, k2p=0.01,
                      tsd=".")
    items = [(mm(1000, 5000), 900, 5000), (mm(1000, 6000), 900, 6000)]
    assert rio.mutual_conflicts(items) == rio.mutual_conflicts(items[::-1])


def test_a_trim_that_would_leave_a_nested_element_outside_is_refused():
    from ltrquest.ltr_model import Member

    def mm(s, e, name):
        return Member(prefix="p", name=f"c:{s}-{e}#{name}", chrom="c", start=s, end=e,
                      l1=s + 99, r0=e - 99, strand="+", orientation="+", family="f", depth=0,
                      k2p=0.01, tsd=".")
    host, inner = mm(100, 5000, "h"), mm(4995, 5000, "i")
    index = rio.SpanIndex([host, inner])
    assert rio.conflict(index, host, 100, 4997) == "nest_broken"
    assert rio.conflict(index, host, 100, 5000) is None


def _mm(s, e, name="x"):
    from ltrquest.ltr_model import Member
    return Member(prefix="p", name=f"c:{s}-{e}#{name}", chrom="c", start=s, end=e,
                  l1=s + 99, r0=e - 99, strand="+", orientation="+", family="f", depth=0,
                  k2p=0.01, tsd=".")


def test_a_trim_onto_the_exact_span_of_a_call_inside_is_refused():
    """The trimmed call would take the key of the call it contains."""
    outer, inner = _mm(1000, 5002, "o"), _mm(1000, 5000, "i")
    index = rio.SpanIndex([outer, inner])
    assert rio.conflict(index, outer, 1000, 5000) == "duplicates_element"
    assert rio.conflict(index, outer, 1000, 5001) is None


def test_two_moves_to_one_span_keep_the_first_claim_in_any_input_order():
    """A call and a second call of it around it: one extends and one trims to the same
    ends. Neither move conflicts with the other call as it stands, so only the moves
    themselves show that both would take one key."""
    outer, inner = _mm(1000, 5002, "o"), _mm(1003, 5000, "i")
    items = [(inner, 1000, 5000), (outer, 1000, 5000)]
    assert rio.conflict(rio.SpanIndex([outer, inner]), inner, 1000, 5000) is None
    assert rio.conflict(rio.SpanIndex([outer, inner]), outer, 1000, 5000) is None
    assert rio.duplicate_moves(items) == {inner.uid}
    assert rio.duplicate_moves(items[::-1]) == {inner.uid}
    assert rio.duplicate_moves([(inner, 1000, 5000), (outer, 1000, 5001)]) == set()


def test_sidecar_v2_round_trips_and_legacy_extended_rows_still_map(tmp_path):
    row = dict.fromkeys(rio.SIDECAR_COLUMNS, ".")
    row.update(old_seq_id="c:150-900#LTR/Gypsy/X", new_seq_id="c:100-897#LTR/Gypsy/X",
               decision="moved", ext5="50", ext3="-3")
    merged = dict(row, old_seq_id="c:600-897#LTR/Gypsy/X", new_seq_id="c:100-897#LTR/Gypsy/X",
                  decision="merged", reason="split_call", merged_into="c:100-897#LTR/Gypsy/X")
    path = tmp_path / "p_reboundary.tsv"
    path.write_text(rio.sidecar_text([row, merged]))
    got = rio.read_map(str(path))
    assert got == {"c:100-897": rio.Rebound("c:150-900#LTR/Gypsy/X",
                                             "c:100-897#LTR/Gypsy/X", 50, -3)}
    legacy = tmp_path / "old_reboundary.tsv"
    legacy.write_text("#old_seq_id\tnew_seq_id\tdecision\text5\text3\n"
                      "c:150-900#X\tc:100-900#X\textended\t50\t0\n")
    assert set(rio.read_map(str(legacy))) == {"c:100-900"}

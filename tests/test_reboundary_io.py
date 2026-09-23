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
    with pytest.raises(ValueError, match="family"):
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
    row = {c: "." for c in rio.SIDECAR_COLUMNS}
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

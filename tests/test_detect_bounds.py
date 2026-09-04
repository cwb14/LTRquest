import pytest

from ltrquest import detect, kmer2ltr

HDR = "\t".join(kmer2ltr.COLUMNS)


def _row(name, seq_len, l5s, l3e, f5, f3):
    v = dict.fromkeys(kmer2ltr.COLUMNS, "NA")
    v.update(seq_id=name, seq_len=str(seq_len), status="pass",
             ltr5_start=str(l5s), ltr3_end=str(l3e),
             flank5_len=str(f5), flank3_len=str(f3))
    return "\t".join(v[c] for c in kmer2ltr.COLUMNS)


def test_assert_bounded_accepts_a_clean_table(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n" + _row("chr1:1-100", 100, 1, 100, 0, 0) + "\n")
    detect.assert_bounded(str(p))


def test_assert_bounded_rejects_residual_flank(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n" + _row("chr1:1-100", 100, 4, 97, 3, 3) + "\n")
    with pytest.raises(AssertionError, match="chr1:1-100"):
        detect.assert_bounded(str(p))


def test_rebase_shifts_coordinates_and_renames(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n"
                 + _row("chr1:100-1100#LTR/x", 1001, 21, 981, 20, 20) + "\n"
                 + _row("chr1:5000-5200#LTR/x", 201, 1, 201, 0, 0) + "\n")
    detect.rebase_to_trimmed(str(p), {"chr1:100-1100#LTR/x": "chr1:120-1080#LTR/x"})
    rows = list(kmer2ltr.read_rows(str(p)))
    assert len(rows) == 1
    r = rows[0]
    assert r["seq_id"] == "chr1:120-1080#LTR/x"
    assert (r["ltr5_start"], r["ltr3_end"], r["seq_len"]) == ("1", "961", "961")
    assert (r["flank5_len"], r["flank3_len"]) == ("0", "0")
    # _row() leaves ltr5_end/ltr3_start unset (kmer2ltr's own "no value" marker);
    # shifting must pass a missing field through unchanged, not arithmetic it
    # into a number or raise.
    assert (r["ltr5_end"], r["ltr3_start"]) == ("NA", "NA")
    detect.assert_bounded(str(p))


def test_rebase_shifts_ltr5_end_and_ltr3_start_too(tmp_path):
    """Every element-relative coordinate moves, not just the bounds that name
    the record -- ltr5_end/ltr3_start feed the GFF3 sub-features."""
    v = dict.fromkeys(kmer2ltr.COLUMNS, "NA")
    v.update(seq_id="chr1:100-1100#LTR/x", seq_len="1001", status="pass",
             ltr5_start="21", ltr5_end="320", ltr3_start="700", ltr3_end="981",
             flank5_len="20", flank3_len="20")
    row = "\t".join(v[c] for c in kmer2ltr.COLUMNS)
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n" + row + "\n")

    detect.rebase_to_trimmed(str(p), {"chr1:100-1100#LTR/x": "chr1:120-1080#LTR/x"})
    r = next(kmer2ltr.read_rows(str(p)))
    assert (r["ltr5_end"], r["ltr3_start"]) == ("300", "680")


def test_rebase_skips_a_truncated_row(tmp_path):
    """A short row can't supply every coordinate column rebase_to_trimmed
    reads; matching relabel_kmer2ltr_tsv's own guard, it is dropped rather
    than raising IndexError, even though its name is one rebase_to_trimmed
    would otherwise carry forward."""
    good = _row("chr1:100-1100#LTR/x", 1001, 21, 981, 20, 20)
    victim = _row("chr1:5000-5200#LTR/x", 201, 1, 201, 0, 0)
    truncated = "\t".join(victim.split("\t")[:5])  # ends before ltr3_start
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n" + good + "\n" + truncated + "\n")

    rename = {"chr1:100-1100#LTR/x": "chr1:120-1080#LTR/x",
              "chr1:5000-5200#LTR/x": "chr1:5001-5200#LTR/x"}
    detect.rebase_to_trimmed(str(p), rename)

    rows = list(kmer2ltr.read_rows(str(p)))
    assert [r["seq_id"] for r in rows] == ["chr1:120-1080#LTR/x"]


def test_classification_relabel_matches_the_table_once_rebased(tmp_path):
    """The seam a classifier crosses after rebase_to_trimmed: cls.tsv classifies
    the bounded record, so its map is keyed the way the rebased table is -- on
    the bounded locus, not the candidate one. Restating cls_names through
    rename's *values* keys it on exactly the names the rebased table carries,
    and relabel_kmer2ltr_tsv succeeds against it without tripping its own
    guard."""
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n" + _row("chr1:100-1100", 1001, 21, 981, 20, 20) + "\n")
    rename = {"chr1:100-1100": "chr1:120-1080"}
    detect.rebase_to_trimmed(str(p), rename)

    cls_names = {"chr1:120-1080": "chr1:120-1080#LTR/Copia/Ale"}
    element_names = detect.rekey_through(
        cls_names, {v: v for v in rename.values()}, "TEBinSorter classifications")
    assert element_names == cls_names

    n_kept, n_dropped = detect.relabel_kmer2ltr_tsv(str(p), element_names)
    assert (n_kept, n_dropped) == (1, 0)
    assert next(kmer2ltr.read_rows(str(p)))["seq_id"] == "chr1:120-1080#LTR/Copia/Ale"


def test_classification_relabel_through_the_candidate_side_mismatches(tmp_path):
    """cls_names is keyed on the bounded locus in both directions of this
    crossing, so restating it through rename's inverse is self-consistent and
    rekey_through raises nothing -- it just lands back on the candidate locus,
    which does not match the rebased table's seq_id. relabel_kmer2ltr_tsv is
    the guard that catches this shape of mismatch."""
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n" + _row("chr1:100-1100", 1001, 21, 981, 20, 20) + "\n")
    rename = {"chr1:100-1100": "chr1:120-1080"}
    detect.rebase_to_trimmed(str(p), rename)

    cls_names = {"chr1:120-1080": "chr1:120-1080#LTR/Copia/Ale"}
    mismatched = detect.rekey_through(
        cls_names, {new: old for old, new in rename.items()},
        "TEBinSorter classifications")
    assert mismatched == {"chr1:100-1100": "chr1:120-1080#LTR/Copia/Ale"}

    with pytest.raises(RuntimeError, match="key elements differently"):
        detect.relabel_kmer2ltr_tsv(str(p), mismatched)


def _bounds_row(seq_id, l5s, l5e, l3s, l3e):
    """A row carrying all four LTR-boundary columns, unlike _row() above
    (which leaves ltr5_end/ltr3_start at kmer2ltr's own "no value" marker)."""
    v = dict.fromkeys(kmer2ltr.COLUMNS, "NA")
    v.update(seq_id=seq_id, seq_len=str(l3e), status="pass",
             ltr5_start=str(l5s), ltr5_end=str(l5e),
             ltr3_start=str(l3s), ltr3_end=str(l3e))
    return "\t".join(v[c] for c in kmer2ltr.COLUMNS)


def test_ltr_bounds_from_table_computes_absolute_coordinates(tmp_path):
    # S=5001; s_lLTR=5001+1-1=5001, e_lLTR=5001+300-1=5300,
    # s_rLTR=5001+702-1=5702, e_rLTR=5001+1001-1=6001.
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n"
                 + _bounds_row("chr3:5001-6001", l5s=1, l5e=300, l3s=702, l3e=1001) + "\n")

    bounds = detect.ltr_bounds_from_table(str(p))
    assert bounds == {"chr3:5001-6001": (5001, 5300, 5702, 6001)}


def test_ltr_bounds_from_table_strips_the_classification_suffix(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n"
                 + _bounds_row("chr3:5001-6001#LTR/Copia/Ale",
                               l5s=1, l5e=300, l3s=702, l3e=1001) + "\n")

    bounds = detect.ltr_bounds_from_table(str(p))
    assert bounds == {"chr3:5001-6001": (5001, 5300, 5702, 6001)}


def test_ltr_bounds_from_table_skips_a_row_with_NA_coordinates(tmp_path):
    """_row() leaves ltr5_end/ltr3_start at kmer2ltr's own "no value" marker;
    a row that cannot supply all four LTR coordinates contributes nothing
    rather than raising."""
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n" + _row("chr1:1-100", 100, 1, 100, 0, 0) + "\n")

    assert detect.ltr_bounds_from_table(str(p)) == {}


def test_ltr_bounds_from_table_keys_match_the_tables_own_row_loci(tmp_path):
    """The regression this guards against: a boundary map built from a
    separate merged SCN file is keyed on each round's pre-trim locus, and
    Kmer2LTR renames every element to its post-trim locus before this table
    is written -- so a map from that source silently stops matching the rows
    it is meant to describe. Deriving the key from the same row as the
    bounds forecloses that: whatever locus a row carries is the key its own
    bounds land under, so this equality holds no matter what the round
    renamed the element to. A derivation that fell back to an SCN-shaped
    source would not touch these rows' own seq_id at all, and would return
    the wrong keys (or none)."""
    rows = [
        _bounds_row("chr2:2000-2500#LTR/Gypsy/Tekay", l5s=1, l5e=200, l3s=301, l3e=500),
        _bounds_row("chr5:100000-100800", l5s=1, l5e=250, l3s=551, l3e=800),
    ]
    p = tmp_path / "t.tsv"
    p.write_text("#" + HDR + "\n" + "\n".join(rows) + "\n")

    own_loci = {r["seq_id"].split("#", 1)[0] for r in kmer2ltr.read_rows(str(p))}
    assert detect.ltr_bounds_from_table(str(p)).keys() == own_loci

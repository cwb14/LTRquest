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

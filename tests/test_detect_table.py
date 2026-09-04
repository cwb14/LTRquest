import textwrap

import pytest

from ltrquest import detect, kmer2ltr

K2L = textwrap.dedent("""\
    {header}
    chr1:100-2100#LTR/unknown\t2001\tpass\t1\t300\t1702\t2001\t300\t300\t0\t0\t305\t300\t4\t2\t5\t0.98\t0.02\t0.0203\t0.008\t250.0\t3.1\t300=\ttg...ca\t338333\t+\tAAGCT\t0,0\tAAGCT
    chr1:5000-5200#LTR/unknown\t201\tno_pair\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t+\tNA\tNA\tNA
    chr1:9000-9400#LTR/unknown\t401\tpass\t1\t80\t322\t401\t80\t80\t0\t0\t82\t80\t1\t0\t2\t0.99\t0.0125\t0.0126\t0.004\t60.0\t2.0\t80=\ttg...ca\t210000\t+\t.\tNA\t.
    """).format(header="\t".join(kmer2ltr.COLUMNS))


def test_element_header_is_thirtythree_columns():
    assert len(detect.ELEMENT_COLUMNS) == 33
    assert detect.ELEMENT_COLUMNS[:29] == kmer2ltr.COLUMNS
    assert detect.ELEMENT_COLUMNS[29:] == ["strand", "family", "domains", "nest_status"]
    assert detect.ELEMENT_HEADER.startswith("#seq_id\t")
    assert detect.ELEMENT_HEADER.rstrip("\n").endswith("nest_status")


def test_filter_drops_non_pass_and_short_ltrs(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text(K2L)
    kept, dropped, malformed = detect.filter_kmer2ltr_in_place(str(p))
    assert (kept, dropped, malformed) == (1, 2, 0)
    text = p.read_text()
    assert text.startswith("seq_id\t") or text.startswith("#seq_id\t")
    assert "chr1:100-2100" in text
    assert "chr1:5000-5200" not in text   # no_pair
    assert "chr1:9000-9400" not in text   # ltr5_len 80 < 100


def test_filter_keeps_the_header(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text(K2L)
    detect.filter_kmer2ltr_in_place(str(p))
    first = p.read_text().splitlines()[0]
    assert "seq_id" in first and "tsd_input" in first


def test_tsd_names_uses_bare_coordinates_and_skips_dot(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text(K2L)
    names = detect.tsd_names_from_kmer2ltr(str(p))
    assert names == {"chr1:100-2100": "AAGCT"}


def test_detect_header_is_thirtyone_columns():
    assert len(detect.DETECT_COLUMNS) == 31
    assert detect.DETECT_COLUMNS[:29] == kmer2ltr.COLUMNS
    assert detect.DETECT_COLUMNS[29:] == ["domains", "nest_status"]


CANDS = (">chr1:100-2100#LTR/unknown\n" + "AC" * 1000 + "A\n"
         ">chr1:5000-5200#LTR/unknown\n" + "GT" * 100 + "G\n"
         ">chr1:9000-9400#LTR/unknown\n" + "TA" * 200 + "T\n")

CLS = ("#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\n"
       "chr1:100-2100\tLTR\tCopia\tAle\tyes\t+\tGAG RT\n"
       "chr1:9000-9400\tLINE\tunknown\tunknown\tyes\t+\tRT\n")


def _k2l_with_flanks(tmp_path):
    """The row for chr1:100-2100 with 20 bp of 5' flank and 30 bp of 3' flank."""
    lines = K2L.splitlines()
    cols = lines[1].split("\t")
    i = kmer2ltr.COLUMNS.index
    cols[i("ltr5_start")], cols[i("ltr3_end")] = "21", "1971"
    cols[i("flank5_len")], cols[i("flank3_len")] = "20", "30"
    p = tmp_path / "k.tsv"
    p.write_text("\n".join([lines[0], "\t".join(cols)] + lines[2:]) + "\n")
    return p


def test_bounded_fasta_cuts_the_flanks_off_and_moves_the_locus(tmp_path):
    tsv = _k2l_with_flanks(tmp_path)
    fa = tmp_path / "c.fa"
    fa.write_text(CANDS)
    out = tmp_path / "b.fa"
    rename = detect.bounded_fasta(str(fa), str(tsv), str(out))

    assert rename["chr1:100-2100#LTR/unknown"] == "chr1:120-2070#LTR/unknown"
    records = dict(detect.iter_fasta(str(out)))
    seq = records["chr1:120-2070#LTR/unknown"]
    assert len(seq) == 1971 - 21 + 1
    assert seq == ("AC" * 1000 + "A")[20:1971]


def test_bounded_fasta_skips_excluded_loci(tmp_path):
    tsv = _k2l_with_flanks(tmp_path)
    fa = tmp_path / "c.fa"
    fa.write_text(CANDS)
    out = tmp_path / "b.fa"
    rename = detect.bounded_fasta(str(fa), str(tsv), str(out),
                                  exclude={"chr1:100-2100"})
    assert "chr1:100-2100#LTR/unknown" not in rename
    assert "chr1:120-2070#LTR/unknown" not in dict(detect.iter_fasta(str(out)))


def test_ltr_names_keeps_only_the_ltr_calls(tmp_path):
    p = tmp_path / "cls.tsv"
    p.write_text(CLS)
    assert detect.ltr_names_from_cls_tsv(str(p)) == {
        "chr1:100-2100": "chr1:100-2100#LTR/Copia/Ale"
    }


def test_relabel_renames_and_drops_the_rest(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text(K2L)
    kept, dropped = detect.relabel_kmer2ltr_tsv(
        str(p), {"chr1:100-2100#LTR/unknown": "chr1:120-2070#LTR/Copia/Ale"})
    assert (kept, dropped) == (1, 2)
    rows = list(kmer2ltr.read_rows(str(p)))
    assert [r["seq_id"] for r in rows] == ["chr1:120-2070#LTR/Copia/Ale"]


def test_bounded_fasta_skips_records_kmer2ltr_could_not_bound(tmp_path):
    tsv = _k2l_with_flanks(tmp_path)
    fa = tmp_path / "c.fa"
    fa.write_text(CANDS)
    out = tmp_path / "b.fa"
    rename = detect.bounded_fasta(str(fa), str(tsv), str(out))

    assert "chr1:5000-5200#LTR/unknown" not in rename      # no_pair
    assert not any(n.startswith("chr1:5000")
                   for n in dict(detect.iter_fasta(str(out))))


def test_bounded_fasta_keeps_an_unparseable_name_but_still_cuts(tmp_path):
    lines = K2L.splitlines()
    cols = lines[1].split("\t")
    i = kmer2ltr.COLUMNS.index
    cols[i("seq_id")] = "contig_with_no_locus"
    cols[i("ltr5_start")], cols[i("ltr3_end")] = "21", "1971"
    cols[i("flank5_len")], cols[i("flank3_len")] = "20", "30"
    tsv = tmp_path / "k.tsv"
    tsv.write_text(lines[0] + "\n" + "\t".join(cols) + "\n")

    fa = tmp_path / "c.fa"
    fa.write_text(">contig_with_no_locus\n" + "AC" * 1000 + "A\n")
    out = tmp_path / "b.fa"
    rename = detect.bounded_fasta(str(fa), str(tsv), str(out))

    assert rename == {"contig_with_no_locus": "contig_with_no_locus"}
    assert dict(detect.iter_fasta(str(out)))["contig_with_no_locus"] == \
        ("AC" * 1000 + "A")[20:1971]


def test_bounded_fasta_rejects_flanks_that_contradict_the_bounds(tmp_path):
    lines = K2L.splitlines()
    cols = lines[1].split("\t")
    i = kmer2ltr.COLUMNS.index
    cols[i("ltr5_start")], cols[i("ltr3_end")] = "21", "1971"
    cols[i("flank5_len")], cols[i("flank3_len")] = "5", "30"   # 5 != 21 - 1
    tsv = tmp_path / "k.tsv"
    tsv.write_text(lines[0] + "\n" + "\t".join(cols) + "\n")

    fa = tmp_path / "c.fa"
    fa.write_text(CANDS)
    with pytest.raises(RuntimeError, match="chr1:100-2100"):
        detect.bounded_fasta(str(fa), str(tsv), str(tmp_path / "b.fa"))


RENAME = {"chr1:100-2100#LTR/unknown": "chr1:120-2070#LTR/unknown",
          "chr1:9000-9400#LTR/unknown": "chr1:9000-9400#LTR/unknown"}


def test_rekey_through_moves_every_key_onto_the_trimmed_locus():
    out = detect.rekey_through({"chr1:100-2100": "AAGCT",
                                "chr1:9000-9400": "TTTT"}, RENAME, "TSDs")
    assert out == {"chr1:120-2070": "AAGCT", "chr1:9000-9400": "TTTT"}


def test_rekey_through_drops_keys_no_element_carried_forward():
    out = detect.rekey_through({"chr1:100-2100": "AAGCT",
                                "chr1:5000-5200": "GGGG"}, RENAME, "TSDs")
    assert out == {"chr1:120-2070": "AAGCT"}


def test_rekey_through_raises_when_nothing_matches():
    with pytest.raises(RuntimeError, match="TSDs"):
        detect.rekey_through({"chr9:1-2": "AAGCT"}, RENAME, "TSDs")


def test_rekey_through_accepts_an_empty_input():
    assert detect.rekey_through({}, RENAME, "TSDs") == {}
    assert detect.rekey_through({"chr1:100-2100": "AAGCT"}, {}, "TSDs") == {}


def test_relabel_raises_when_the_names_key_on_something_else(tmp_path):
    p = tmp_path / "k.tsv"
    p.write_text(K2L)
    with pytest.raises(RuntimeError, match="key elements differently"):
        detect.relabel_kmer2ltr_tsv(str(p), {"chr9:1-2": "chr9:1-2#LTR/Copia/Ale"})


def test_bounded_fasta_keeps_one_record_per_bounded_span(tmp_path):
    """Two candidates that bound to the same span are one element found twice."""
    hdr = K2L.splitlines()[0]
    i = kmer2ltr.COLUMNS.index

    def row(name, seq_len, l5s, l3e, f5, f3):
        cols = K2L.splitlines()[1].split("\t")
        cols[i("seq_id")], cols[i("seq_len")] = name, str(seq_len)
        cols[i("ltr5_start")], cols[i("ltr3_end")] = str(l5s), str(l3e)
        cols[i("flank5_len")], cols[i("flank3_len")] = str(f5), str(f3)
        return "\t".join(cols)

    tsv = tmp_path / "k.tsv"
    tsv.write_text("\n".join([hdr,
                              row("chr1:100-2100", 2001, 21, 1971, 20, 30),
                              row("chr1:110-2080", 1971, 11, 1961, 10, 10)]) + "\n")
    fa = tmp_path / "c.fa"
    fa.write_text(">chr1:100-2100\n" + "AC" * 1000 + "A\n"
                  ">chr1:110-2080\n" + "AC" * 985 + "A\n")
    out = tmp_path / "b.fa"
    rename = detect.bounded_fasta(str(fa), str(tsv), str(out))

    assert rename == {"chr1:100-2100": "chr1:120-2070"}
    assert list(dict(detect.iter_fasta(str(out)))) == ["chr1:120-2070"]
    assert len(set(rename.values())) == len(rename)


# A good element and a short TSD+ neighbour that overlaps it almost exactly and
# shares both LTRs, so the pre-purge has nothing to protect it as a nested TE.
SCN = ("1000 3000 2001 1000 1200 201 2800 3000 201 95.0 0 chr1\n"
       "1050 2950 1901 1050 1200 151 2800 2950 151 95.0 0 chr1\n")


def _pair_table(tmp_path):
    i = kmer2ltr.COLUMNS.index

    def row(name, seq_len, ltr_len, tsd):
        cols = ["NA"] * len(kmer2ltr.COLUMNS)
        cols[i("seq_id")], cols[i("seq_len")], cols[i("status")] = name, str(seq_len), "pass"
        cols[i("ltr5_len")] = cols[i("ltr3_len")] = cols[i("aln_len")] = str(ltr_len)
        cols[i("tsd")] = tsd
        return "\t".join(cols)

    p = tmp_path / "k.tsv"
    p.write_text("\n".join([K2L.splitlines()[0],
                            row("chr1:1000-3000", 2001, 201, "."),      # keeps
                            row("chr1:1050-2950", 1901, 80, "AAGCT")])  # too short
                 + "\n")
    return p


def test_a_purger_the_length_filter_drops_cannot_take_a_neighbour_with_it(tmp_path):
    scn = tmp_path / "s.scn"
    scn.write_text(SCN)
    bounds = detect.load_scn_ltr_boundaries(str(scn))
    tsv = _pair_table(tmp_path)

    # Unfiltered, the short TSD+ element dominates its good neighbour.
    unfiltered = detect.tsd_names_from_kmer2ltr(str(tsv))
    assert detect.pre_purge_tsd_dominated(
        str(scn), set(unfiltered), threshold=0.80, ltr_bounds=bounds
    ) == {"chr1:1000-3000"}

    # Gating on the retained elements leaves nothing able to purge it.
    detect.filter_kmer2ltr_in_place(str(tsv))
    retained = detect.tsd_names_from_kmer2ltr(str(tsv))
    assert retained == {}
    assert detect.pre_purge_tsd_dominated(
        str(scn), set(retained), threshold=0.80, ltr_bounds=bounds
    ) == set()


BOUNDED_FA = (">chr1:120-2070\n" + "AC" * 1000 + "A\n"
              ">chr1:9000-9400\n" + "TA" * 200 + "T\n")

CLS_WITH_STRAND = ("#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\n"
                   "chr1:120-2070\tLTR\tCopia\tAle\tyes\t-\tGAG RT\n"
                   "chr1:9000-9400\tLTR\tGypsy\tTat\tyes\t+\tGAG RT\n")

LIB_RENAME = {"chr1:120-2070": "chr1:120-2070#LTR/Copia/Ale",
             "chr1:9000-9400": "chr1:9000-9400#LTR/Gypsy/Tat"}


def test_bounded_fasta_oriented_flips_only_the_minus_strand_record(tmp_path):
    fa = tmp_path / "bounded.fa"
    fa.write_text(BOUNDED_FA)
    cls = tmp_path / "cls.tsv"
    cls.write_text(CLS_WITH_STRAND)
    out = tmp_path / "lib.fa"

    flipped = detect.bounded_fasta_oriented(
        str(fa), str(cls), str(out), set(LIB_RENAME.values()), LIB_RENAME)

    assert flipped == {"chr1:120-2070#LTR/Copia/Ale"}
    records = dict(detect.iter_fasta(str(out)))
    assert records["chr1:120-2070#LTR/Copia/Ale"] == detect.revcomp("AC" * 1000 + "A")
    assert records["chr1:9000-9400#LTR/Gypsy/Tat"] == "TA" * 200 + "T"


def test_bounded_fasta_oriented_respects_keep_names(tmp_path):
    fa = tmp_path / "bounded.fa"
    fa.write_text(BOUNDED_FA)
    cls = tmp_path / "cls.tsv"
    cls.write_text(CLS_WITH_STRAND)
    out = tmp_path / "lib.fa"

    flipped = detect.bounded_fasta_oriented(
        str(fa), str(cls), str(out), {"chr1:9000-9400#LTR/Gypsy/Tat"}, LIB_RENAME)

    assert flipped == set()
    assert list(dict(detect.iter_fasta(str(out)))) == ["chr1:9000-9400#LTR/Gypsy/Tat"]


def test_bounded_fasta_oriented_defaults_to_forward_when_cls_tsv_has_no_row(tmp_path):
    fa = tmp_path / "bounded.fa"
    fa.write_text(">chr1:1-50\n" + "A" * 30 + "T" * 20 + "\n")
    cls = tmp_path / "cls.tsv"
    cls.write_text("#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\n")
    out = tmp_path / "lib.fa"

    flipped = detect.bounded_fasta_oriented(str(fa), str(cls), str(out), {"chr1:1-50"}, {})

    assert flipped == set()
    assert dict(detect.iter_fasta(str(out)))["chr1:1-50"] == "A" * 30 + "T" * 20

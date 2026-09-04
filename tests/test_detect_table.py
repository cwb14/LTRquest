import textwrap

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


def test_annotate_inserts_between_kmer2ltr_and_domains():
    i = detect.ELEMENT_COLUMNS.index("strand")
    assert detect.ELEMENT_COLUMNS[i - 1] == "tsd_input"
    assert detect.ELEMENT_COLUMNS[i + 1] == "family"
    assert detect.ELEMENT_COLUMNS[i + 2] == "domains"


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

from ltrquest import flag_fp

HEADER = "#seq_id\tseq_len\tstatus\tltr5_len\tk2p\ttsd\tdomains\tnest_status\n"
ROWS = (
    "chr1:100-2100#LTR/Gypsy/Tekay\t2001\tpass\t300\t0.02\tAAGCT"
    "\tRT|Tekay@500-700\tnest-outer:chr1:900-1500\n"
    "chr1:900-1500#LTR/Copia/Ale\t601\tpass\t120\t0.01\t.\t.\tnest-inner:chr1:100-2100\n"
)


def test_load_depth_tsvs_reads_the_documented_fields(tmp_path):
    p = tmp_path / "s_depth0_ltr.tsv"
    p.write_text(HEADER + ROWS)
    recs = flag_fp.load_depth_tsvs([str(p)])
    a = recs["chr1:100-2100#LTR/Gypsy/Tekay"]
    assert a.ltr_len == 300
    assert abs(a.k2p - 0.02) < 1e-9
    assert a.superfamily == "Gypsy"
    assert a.raw_domains == [("RT", 500, 700)]
    assert a.insertions == [(900, 1500)]
    b = recs["chr1:900-1500#LTR/Copia/Ale"]
    assert b.raw_domains == []
    assert b.insertions == []


def test_load_depth_tsvs_survives_a_reordered_schema(tmp_path):
    """The load-bearing test: `domains` is no longer second-from-last and
    `ltr5_len` is no longer column 2, so a positional reader gets this wrong."""
    p = tmp_path / "s_depth0_ltr.tsv"
    p.write_text("#seq_id\tk2p\tdomains\tnest_status\tltr5_len\tstatus\n"
                 "chr1:100-2100#LTR/Gypsy/Tekay\t0.02\tRT|Tekay@500-700\t.\t300\tpass\n")
    recs = flag_fp.load_depth_tsvs([str(p)])
    e = recs["chr1:100-2100#LTR/Gypsy/Tekay"]
    assert e.ltr_len == 300
    assert abs(e.k2p - 0.02) < 1e-9
    assert e.raw_domains == [("RT", 500, 700)]


def test_clean_depth_tsv_drops_fp_rows_and_scrubs_dangling_nest(tmp_path):
    p = tmp_path / "s_depth0_ltr.tsv"
    p.write_text(HEADER + ROWS)
    out = tmp_path / "s_depth0_clean_ltr.tsv"
    removed, scrubbed = flag_fp.clean_depth_tsv(str(p), str(out), {"chr1:100-2100"})
    text = out.read_text()
    assert removed == 1
    assert scrubbed >= 1
    assert "chr1:100-2100" not in text
    assert "chr1:900-1500" in text
    assert "nest-inner" not in text
    assert text.startswith("#seq_id")


def _metrics(dominance, recon):
    return flag_fp.FamilyMetrics(rep="r", n=10, reconstitution=recon,
                                 dominance=dominance, entropy=0.0,
                                 cross_superfamily=False, n_unknown=0, n_mixture=0)


def test_classify_family_gates():
    assert flag_fp.classify_family(_metrics(dominance=0.9, recon=0.0)) == "safe"
    assert flag_fp.classify_family(_metrics(dominance=0.1, recon=0.9)) == "recovered"
    assert flag_fp.classify_family(_metrics(dominance=0.1, recon=0.1)) == "false_positive"

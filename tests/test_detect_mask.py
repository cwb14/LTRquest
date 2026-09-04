from ltrquest import detect

NAME = "chr1:101-140#LTR/Gypsy"
HDR = "#seq_id\tdomains\tnest_status\n"


def _write(tmp_path, seq):
    fa = tmp_path / "lib.fa"
    fa.write_text(f">{NAME}\n{seq}\n")
    tsv = tmp_path / "t.tsv"
    tsv.write_text(HDR + f"{NAME}\t.\tnest-outer:chr1:131-140\n")
    return str(fa), str(tsv)


def _seq(fa_path):
    return "".join(line.strip() for line in open(fa_path) if not line.startswith(">"))


def test_forward_record_masks_the_last_ten_bases(tmp_path):
    fa, tsv = _write(tmp_path, "A" * 30 + "C" * 10)
    detect.mask_same_round_inners_in_fa(fa, tsv, "R")
    assert _seq(fa) == "A" * 30 + "R" * 10


def test_revcomped_record_masks_the_first_ten_bases(tmp_path):
    """The outer's genomic interval (its last ten bases) is unchanged, but
    this record is one `bounded_fasta_oriented` reverse-complemented, so
    those same ten bases are now this record's first ten."""
    fa, tsv = _write(tmp_path, "G" * 10 + "T" * 30)
    detect.mask_same_round_inners_in_fa(fa, tsv, "R", revcomped={NAME})
    assert _seq(fa) == "R" * 10 + "T" * 30
